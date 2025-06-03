#include "ConfigurationGenerator.h"

ConfigurationGenerator::ConfigurationGenerator(const SimulationParameters &sim_params){

    init_config_index = sim_params.init_config_index;

    init_M = sim_params.init_M;

    // Will be set during equilibration
    average_cumulative_loop_size = 0;
    N_l = 0;
    max_loop_size = 0;
    n = 0;
    M = init_M;

    beta = sim_params.beta;

    N = sim_params.N;

    a_parameter = sim_params.new_M_multiplier;

    equilibration_steps = sim_params.equilibration_steps;
    mc_steps = sim_params.simulation_steps;

    initial_config_generator.seed(init_config_seed);
    diagonal_update_generator.seed(diagonal_update_seed);
    off_diagonal_update_generator.seed(off_diagonal_update_seed);
    disconnected_spin_flip_generator.seed(disconnected_spin_flip_seed);
    loop_start_pos_generator.seed(loop_start_pos_seed);
    metropolis_generator.seed(metropolis_seed);

    spin_labels = {-1,1};

    ConfigurationGenerator::setInitialConfiguration();
    ConfigurationGenerator::populateOperatorLocations(M);
}

void ConfigurationGenerator::populateOperatorLocations(const int &num_fill_zeros) {

    for (int i=0; i<num_fill_zeros; i++){
        operator_locations.push_back(0);
    }
}

void ConfigurationGenerator::setInitialConfiguration() {

    if (init_config_index == 0){
        for (int i=0; i<N; i++){
            if (i % 2 == 0){
                spin_configuration.push_back(1);
            }
            else {
                spin_configuration.push_back(-1);
            }
        }
    }
    else if (init_config_index == 2){
        std::uniform_int_distribution<> initial_config_dist(0, 1);
        for (int i=0; i<N; i++){
            spin_configuration.push_back(spin_labels.at(initial_config_dist(initial_config_generator)));
        }
    }
    else if (init_config_index == 1){
        for (int i=0; i<N; i++){
            spin_configuration.push_back(1);
        }
    }
    else if (init_config_index == -1){
        for (int i=0; i<N; i++){
            spin_configuration.push_back(-1);
        }
    }
    else{
        throw std::runtime_error("Invalid initial configuration index");
    }
}

void ConfigurationGenerator::diagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types)
{

    p_list.clear();
    b_list.clear();

    std::uniform_real_distribution<double> metropolis_acceptance_distribution(0.0,1.0);
    std::discrete_distribution<int> diag_update_distribution(prob_tables.max_norm_probabilities.begin(),
                                                             prob_tables.max_norm_probabilities.end());

    for (int t=0; t<M; t++){
        double M_minus_n = M - n;
        if (operator_locations.at(t) == 0) {
            double u_1 = metropolis_acceptance_distribution(metropolis_generator);
            double acceptance_prob = beta * prob_tables.max_diagonal_norm/M_minus_n;
            if (u_1 < acceptance_prob) {
                int b_prime = diag_update_distribution(diagonal_update_generator);
                int i_b = prob_tables.lattice_sites.at(b_prime).at(0);
                int j_b = prob_tables.lattice_sites.at(b_prime).at(1);
                std::vector<int> diag_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                spin_configuration.at(i_b), spin_configuration.at(j_b)};
                int vertex_type = vertex_types.getVertexTypeIndex(diag_config);
                double u_2 = metropolis_acceptance_distribution(metropolis_generator);
                double acceptance_prob_2 = prob_tables.diagonal_probabilities.at(vertex_type).at(b_prime);
                if (u_2 < acceptance_prob_2) {
                    int bond_num = b_prime + 1;
                    operator_locations.at(t) = 2*bond_num;
                    n++;
                    p_list.push_back(t);
                    b_list.push_back(bond_num);
                }
            }
        }
        else if (operator_locations.at(t) % 2 == 1){
            p_list.push_back(t);
            int b = (operator_locations.at(t) - 1)/2;
            b_list.push_back(b);
            int bond_index = b - 1;
            int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
            int j_b = prob_tables.lattice_sites.at(bond_index).at(1);
            spin_configuration.at(i_b) = -spin_configuration.at(i_b);
            spin_configuration.at(j_b) = -spin_configuration.at(j_b);
        }
        else {
            double acceptance_prob = (M_minus_n + 1.0)/(beta * prob_tables.max_diagonal_norm);
            double u_1 = metropolis_acceptance_distribution(metropolis_generator);
            if (u_1 < acceptance_prob) {
                operator_locations.at(t) = 0;
                n--;
            }
            else {
                p_list.push_back(t);
                int b = operator_locations.at(t)/2;
                b_list.push_back(b);
            }
        }
    }
}

void ConfigurationGenerator::initializeOffDiagonalUpdates() {

    first_vertex_leg.clear();
    last_vertex_leg.clear();
    linked_list.clear();
    vertex_configuration.clear();

    num_total_legs = num_legs_per_vertex * n;

    for (int i=0; i<N; i++){
        first_vertex_leg.push_back(-1);
        last_vertex_leg.push_back(-1);
    }

    for (int i=0; i<n; i++){
        vertex_configuration.push_back(0);
        for (int j=0; j<num_legs_per_vertex; j++){
            linked_list.push_back(0);
        }
    }
}

void ConfigurationGenerator::offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types) {

    if (n == 0) {
        std::uniform_int_distribution<> spin_flip_dist(0, 1);
        for (int i = 0; i < N; i++) {
            int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
            spin_configuration.at(i) = multiplier * spin_configuration.at(i);
        }
    }
    else {
        ConfigurationGenerator::initializeOffDiagonalUpdates();
        for (int p=0; p<n; p++){
            int t = p_list.at(p);
            int bond_num = b_list.at(p);
            int bond_index = bond_num - 1;
            int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
            int j_b = prob_tables.lattice_sites.at(bond_index).at(1);

            if (last_vertex_leg.at(i_b) == -1) {
                first_vertex_leg.at(i_b) = num_legs_per_vertex*p;
                last_vertex_leg.at(i_b) = num_legs_per_vertex*p + 2;
            }
            else{
                linked_list.at(num_legs_per_vertex*p) = last_vertex_leg.at(i_b);
                linked_list.at(last_vertex_leg.at(i_b)) = num_legs_per_vertex*p;
                last_vertex_leg.at(i_b) = num_legs_per_vertex*p + 2;
            }

            if (last_vertex_leg.at(j_b) == -1) {
                first_vertex_leg.at(j_b) = num_legs_per_vertex*p + 1;
                last_vertex_leg.at(j_b) = num_legs_per_vertex*p + 3;
            }
            else {
                linked_list.at(num_legs_per_vertex*p + 1) = last_vertex_leg.at(j_b);
                linked_list.at(last_vertex_leg.at(j_b)) = num_legs_per_vertex*p + 1;
                last_vertex_leg.at(j_b) = num_legs_per_vertex*p + 3;
            }

            if (operator_locations.at(t) % 2 == 1) {
                std::vector<int> current_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                -spin_configuration.at(i_b), -spin_configuration.at(j_b)};
                vertex_configuration.at(p) = vertex_types.getVertexTypeIndex(current_config);
                spin_configuration.at(i_b) = -spin_configuration.at(i_b);
                spin_configuration.at(j_b) = -spin_configuration.at(j_b);
            }
            else {
                std::vector<int> current_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                spin_configuration.at(i_b), spin_configuration.at(j_b)};
                vertex_configuration.at(p) = vertex_types.getVertexTypeIndex(current_config);
            }
        }

        for (int i=0; i<N; i++){
            if (last_vertex_leg.at(i) != -1) {
                linked_list.at(last_vertex_leg.at(i)) = first_vertex_leg.at(i);
                linked_list.at(first_vertex_leg.at(i)) = last_vertex_leg.at(i);
            }
        }

        bool skip_loop_update = false;
        int loop_size = 0;
        cumulative_loop_size = 0;

        for (int loop_num=0; loop_num<N_l; loop_num++){

            if (loop_size >= max_loop_size) {
                skip_loop_update = true;
                break;
            }

            cumulative_loop_size += loop_size;
            loop_size = 0;

            std::uniform_int_distribution<> loop_start_pos_dist(0, num_total_legs - 1);
            int j_0 = loop_start_pos_dist(loop_start_pos_generator);
            int j_current = j_0;

            while (loop_size < max_loop_size) {

                int l_e = j_current % num_legs_per_vertex;
                int p = (j_current - l_e)/num_legs_per_vertex;
                int vertex_type = vertex_configuration.at(p);

                int composite_leg_index = num_legs_per_vertex*l_e;
                int row_index = num_composite_leg_indices*vertex_type + composite_leg_index;
                int loc_index = num_legs_per_vertex*(num_legs_per_vertex-1)*vertex_type + (num_legs_per_vertex-1)*l_e;

                std::vector<double> out_leg_probs;
                for (int k=0; k<num_legs_per_vertex-1; k++){
                    int l_k = vertex_types.allowed_exit_legs.at(loc_index + k);
                    double prob_l_k = prob_tables.loop_update_probabilities.at(row_index + l_k).at(b_list.at(p)-1);
                    out_leg_probs.push_back(prob_l_k);
                }

                std::discrete_distribution<int> exit_leg_dist(out_leg_probs.begin(),out_leg_probs.end());
                int l_x_index = exit_leg_dist(off_diagonal_update_generator);
                int l_x = vertex_types.allowed_exit_legs.at(loc_index + l_x_index);
                vertex_configuration.at(p) = vertex_types.getFlippedSpinsVertexIndex(l_e, l_x, vertex_type);

                j_current = num_legs_per_vertex*p + l_x;
                if (l_e != l_x)
                    loop_size++;
                if (j_current == j_0){
                    break;
                }
                else {
                    j_current = linked_list.at(j_current);
                    if (j_current == j_0) {
                        break;
                    }
                }
            }
        }

        bool free_spins = false;
        std::vector<int> free_spin_indices;
        std::uniform_int_distribution<> spin_flip_dist(0, 1);

        if (!skip_loop_update) {
            for (int p=0; p<n; p++){
                int t = p_list.at(p);
                int b = b_list.at(p);
                operator_locations.at(t) = 2*b + vertex_types.is_off_diag.at(vertex_configuration.at(p));
            }
            for (int i=0; i<N; i++){
                if (first_vertex_leg.at(i) == -1){
                    int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                    spin_configuration.at(i) = multiplier * spin_configuration.at(i);
                    free_spins = true;
                    free_spin_indices.push_back(i);
                }
                else{
                    int l = first_vertex_leg.at(i) % num_legs_per_vertex;
                    int p = (first_vertex_leg.at(i) - l) / num_legs_per_vertex;
                    int vertex_type = vertex_configuration.at(p);
                    std::vector<int> config = vertex_types.getVertexConfig(vertex_type);
                    spin_configuration.at(i) = config.at(l);
                }
            }
            if (free_spins) {
                for (int i=0; i<num_free_flips; i++){
                    for (int j=0; j<free_spin_indices.size(); j++){
                        int s = free_spin_indices.at(j);
                        int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                        spin_configuration.at(s) = multiplier * spin_configuration.at(s);
                    }
                }
            }
        }
        for (int i=0; i<N; i++) {
            if (first_vertex_leg.at(i) == -1) {
                int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                spin_configuration.at(i) = multiplier * spin_configuration.at(i);
                free_spins = true;
                free_spin_indices.push_back(i);
            }
        }
        if (free_spins) {
            for (int i=0; i<num_free_flips; i++){
                for (int j=0; j<free_spin_indices.size(); j++){
                    int s = free_spin_indices.at(j);
                    int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                    spin_configuration.at(s) = multiplier * spin_configuration.at(s);
                }
            }
        }
    }
}

void ConfigurationGenerator::simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables,
                                                            Estimators &estimators,
                                                            const SimulationParameters &sim_params,
                                                            const VertexTypes &vertex_types) {

    int max_n = n;
    int n_sum = 0;
    int cumulative_loop_size_sum = 0;
    for (int step=0; step<equilibration_steps; step++){
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        if (n > max_n) {
            max_n = n;
        }
        n_sum += n;
        int M_new = std::max((int)(a_parameter * max_n),M);
        ConfigurationGenerator::populateOperatorLocations(int (M_new - M));
        M = M_new;
        int avg_n = n_sum / step;
        max_loop_size = 100*avg_n;
        N_l = 2*M;
        num_free_flips = N_l;
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        cumulative_loop_size_sum += cumulative_loop_size;
    }

    average_cumulative_loop_size =  cumulative_loop_size_sum / equilibration_steps;
    N_l = std::max(2*M/average_cumulative_loop_size, init_M);
    num_free_flips = N_l;

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        estimators.updateAllProperties(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset);
    }

    estimators.outputStepData(sim_params);
}


