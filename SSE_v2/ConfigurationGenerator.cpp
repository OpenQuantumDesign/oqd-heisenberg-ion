#include "ConfigurationGenerator.h"

ConfigurationGenerator::ConfigurationGenerator(const SimulationParameters &sim_params, const ProbabilityTables &prob_tables){

    init_config_index = sim_params.init_config_index;

    init_M = (int)(2 * sim_params.beta * prob_tables.max_diagonal_norm);

    // Will be set during equilibration
    average_cumulative_loop_size = 1;
    max_loop_size = 10;
    n = 0;
    M = init_M;
    N_l = 1;
    num_clusters = 100;
    count_non_skipped_loop_updates = 0;
    skip_loop_update = true;

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
    metropolis_generator_1.seed(metropolis_seed_1);
    metropolis_generator_2.seed(metropolis_seed_2);
    metropolis_generator_3.seed(metropolis_seed_3);

    spin_labels = {-1,1};

    ConfigurationGenerator::setInitialConfiguration();
    ConfigurationGenerator::populateOperatorLocations(M);

    num_winding = 0;
}

void ConfigurationGenerator::populateOperatorLocations(const int &num_fill_zeros) {

    std::vector<int> operator_locations_new;
    for (int i=0;i<n;i++){
        operator_locations_new.push_back(operator_locations.at(p_list.at(i)));
        p_list[i] = i;
    }
    operator_locations = operator_locations_new;
    for (int i=0; i<num_fill_zeros; i++){
        operator_locations.push_back(0);
    }
}

void ConfigurationGenerator::setInitialConfiguration() {

    if (init_config_index == 0){
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
    else if (init_config_index > 1){
        for (int i=0; i<N; i++){
            if (i % init_config_index == 0){
                spin_configuration.push_back(1);
            }
            else {
                spin_configuration.push_back(-1);
            }
        }
    }
    else {
        throw std::runtime_error("Invalid initial configuration index provided\n");
    }
    std::cout << init_config_index << "\n";
    for (int spin_index = 0; spin_index < N; spin_index++) {
        std::cout << spin_configuration.at(spin_index) << ",";
    }
    std::cout << "\n";
}

void ConfigurationGenerator::diagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types)
{

    p_list.clear();
    b_list.clear();
    id_list.clear();

    //std::uniform_real_distribution<double> metropolis_acceptance_distribution(0.0,1.0);
    std::discrete_distribution<int> diag_update_distribution(prob_tables.max_norm_probabilities.begin(),
                                                             prob_tables.max_norm_probabilities.end());

    for (int t=0; t<M; t++){
        double M_minus_n = M - n;
        if (operator_locations.at(t) == 0) {
            double u_1 = metropolis_acceptance_distribution(metropolis_generator_1);
            double acceptance_prob = beta * prob_tables.max_diagonal_norm/M_minus_n;
            if (u_1 < acceptance_prob) {
                int b_prime = diag_update_distribution(diagonal_update_generator);
                int i_b = prob_tables.lattice_sites.at(b_prime).at(0);
                int j_b = prob_tables.lattice_sites.at(b_prime).at(1);
                std::vector<int> diag_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                spin_configuration.at(i_b), spin_configuration.at(j_b)};
                int vertex_type = vertex_types.getVertexTypeIndex(diag_config);
                double u_2 = metropolis_acceptance_distribution(metropolis_generator_2);
                double acceptance_prob_2 = prob_tables.diagonal_probabilities.at(vertex_type).at(b_prime);
                if (u_2 < acceptance_prob_2) {
                    int bond_num = b_prime + 1;
                    operator_locations.at(t) = 2*bond_num;
                    n++;
                    p_list.push_back(t);
                    b_list.push_back(bond_num);
                    if (vertex_type == 1 or vertex_type == 2){
                        disorder_vertex_locations.push_back(t);
                    }
                }
                else{
                    id_list.push_back(t);
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
            double u_1 = metropolis_acceptance_distribution(metropolis_generator_3);
            if (u_1 < acceptance_prob) {
                operator_locations.at(t) = 0;
                n--;
                id_list.push_back(t);
            }
            else {
                p_list.push_back(t);
                int b = operator_locations.at(t)/2;
                b_list.push_back(b);
                int i_b = prob_tables.lattice_sites.at(b-1).at(0);
                int j_b = prob_tables.lattice_sites.at(b-1).at(1);
                std::vector<int> diag_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                spin_configuration.at(i_b), spin_configuration.at(j_b)};
                int vertex_type = vertex_types.getVertexTypeIndex(diag_config);
                if (vertex_type == 1 or vertex_type == 2){
                    disorder_vertex_locations.push_back(t);
                }
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
    loop_start_pos_dist = std::uniform_int_distribution<>(0, num_total_legs - 1);

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

void ConfigurationGenerator::multiBranchClusterUpdate(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types){

    if (n == 0) {
        for (int i = 0; i < N; i++) {
            int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
            spin_configuration.at(i) = multiplier * spin_configuration.at(i);
        }
        num_free_spins = N;
    }
    else {
        int num_disorder_vertex_indices = disorder_vertex_locations.size()-1;
        std::uniform_int_distribution<> distrib(0, num_disorder_vertex_indices);
        ConfigurationGenerator::populateLinkedList(prob_tables, vertex_types);
        //for (int i = 0; i<num_clusters; i++) {

        //}

    }
}

void ConfigurationGenerator::populateLinkedList(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types) {

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
}

void ConfigurationGenerator::offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types) {

    if (n == 0) {
        for (int i = 0; i < N; i++) {
            int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
            spin_configuration.at(i) = multiplier * spin_configuration.at(i);
        }
        num_free_spins = N;
    }
    else {
        populateLinkedList(prob_tables, vertex_types);
        int loop_size = 0;
        cumulative_loop_size = 0;

        for (int loop_num=0; loop_num<N_l; loop_num++){

            skip_loop_update = true;
            loop_size = 0;

            int j_0 = loop_start_pos_dist(loop_start_pos_generator);
            int j_current = j_0;

            for (int iter=0;iter<max_loop_size;iter++){

                int l_e = j_current % num_legs_per_vertex;
                int p = (j_current - l_e)/num_legs_per_vertex;
                int vertex_type = vertex_configuration.at(p);

                int composite_leg_index = num_legs_per_vertex*l_e;
                int row_index = num_composite_leg_indices*vertex_type + composite_leg_index;
                int loc_index = num_legs_per_vertex*(num_legs_per_vertex-1)*vertex_type + (num_legs_per_vertex-1)*l_e;

                out_leg_probs.clear();
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
                if (l_e != l_x) {
                    loop_size++;
                }
                if (j_current == j_0){
                    skip_loop_update = false;
                    cumulative_loop_size += loop_size;
                    break;
                }
                else {
                    j_current = linked_list.at(j_current);
                    if (j_current == j_0) {
                        skip_loop_update = false;
                        cumulative_loop_size += loop_size;
                        break;
                    }
                }
            }

            if (skip_loop_update) {
                cumulative_loop_size = 0;
                break;
            }
        }

        bool free_spins = false;
        num_free_spins = 0;
        std::vector<int> free_spin_indices;

        if (!skip_loop_update) {
            count_non_skipped_loop_updates++;
            num_winding = 0;
            for (int p = 0; p < n; p++) {
                int t = p_list.at(p);
                int b = b_list.at(p);
                operator_locations.at(t) = 2 * b + vertex_types.is_off_diag.at(vertex_configuration.at(p));
                num_winding += vertex_types.twist_mapping.at(vertex_configuration.at(p));
            }
            for (int i = 0; i < N; i++) {
                if (first_vertex_leg.at(i) == -1) {
                    int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                    spin_configuration.at(i) = multiplier * spin_configuration.at(i);
                    free_spins = true;
                    free_spin_indices.push_back(i);
                    num_free_spins++;
                } else {
                    int l = first_vertex_leg.at(i) % num_legs_per_vertex;
                    int p = (first_vertex_leg.at(i) - l) / num_legs_per_vertex;
                    int vertex_type = vertex_configuration.at(p);
                    std::vector<int> config = vertex_types.getVertexConfig(vertex_type);
                    spin_configuration.at(i) = config.at(l);
                }
            }
        }
        else {
            for (int i=0; i<N; i++) {
                if (first_vertex_leg.at(i) == -1) {
                    int multiplier = spin_labels.at(spin_flip_dist(disconnected_spin_flip_generator));
                    spin_configuration.at(i) = multiplier * spin_configuration.at(i);
                    free_spins = true;
                    free_spin_indices.push_back(i);
                    num_free_spins++;
                }
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

void ConfigurationGenerator::randomSpinFlipsXXZ() {

    int flip_spins_bool = spin_flip_dist(disconnected_spin_flip_generator);
    if (flip_spins_bool == 1){
        for (int j=0; j<N; j++){
            spin_configuration.at(j) = -1 * spin_configuration.at(j);
        }
    }
}

void ConfigurationGenerator::flipAllSpins() {
    for (int j=0; j<N; j++){
        spin_configuration.at(j) = -1 * spin_configuration.at(j);
    }
}

void ConfigurationGenerator::simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables,
                                                            Estimators &estimators,
                                                            const SimulationParameters &sim_params,
                                                            const VertexTypes &vertex_types) {

    int max_n = n;
    int n_sum = 0;
    int cumulative_loop_size_sum = 0;
    int avg_n;
    for (int step=0; step<equilibration_steps; step++){
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        if (n > max_n) {
            max_n = n;
        }
        n_sum += n;
        int M_new = std::max((int)(a_parameter * max_n),M);
        //int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),M);
        ConfigurationGenerator::populateOperatorLocations(int (M_new - n));
        M = M_new;
        avg_n = n_sum / step;
        average_cumulative_loop_size = cumulative_loop_size_sum / count_non_skipped_loop_updates;
        max_loop_size = std::max(100*avg_n,1);
        N_l = std::max(2*avg_n/average_cumulative_loop_size, 1);
        num_free_flips = N_l;
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        cumulative_loop_size_sum += cumulative_loop_size;
    }

    average_cumulative_loop_size =  cumulative_loop_size_sum / count_non_skipped_loop_updates;
    N_l = std::max(2*avg_n/average_cumulative_loop_size, 1);
    num_free_flips = N_l;
    max_loop_size = std::max(100*avg_n,1);

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        estimators.updateAllProperties(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset, num_winding);
    }

    estimators.outputStepData(sim_params);
}

void ConfigurationGenerator::simulateProbabilisticLoopsXXZ(const ProbabilityTables &prob_tables,
                                                            Estimators &estimators,
                                                            const SimulationParameters &sim_params,
                                                            const VertexTypes &vertex_types) {

    int max_n = n;
    int avg_n = 0;
    double avg_num_samples = 1000.0;
    int n_sum = 0.0;
    int cumulative_loop_size_sum = 1;

    double avg_n_d = 0.0;
    double avg_cumul_loop_size_d = 1.0;
    double prob_non_skip_update = 1.0;
    double N_l_d = 1.0;

    for (int step=0; step<equilibration_steps; step++) {
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        if (n > max_n) {
            max_n = n;
        }
        n_list.push_back(n);
        n_sum += n;

        avg_n_d = n_sum / ((double)step+1.0);
        avg_n = std::round(avg_n_d);

        //int M_new = std::max((int) (a_parameter * max_n), M);
        /*
        int M_new = (int)(a_parameter * avg_n_d);
        if (M_new - n <= 0){
            M_new = n+1;
        }
         */
        //int M_new = std::max(2*n,M);
        int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M = M_new;

        max_loop_size = 10000 * std::max(avg_n, 1);

        /*
        cumulative_loop_size_sum += cumulative_loop_size;
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        avg_cumul_loop_size_d = cumulative_loop_size_sum/((double)step+1.0);
        average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
         */

        //N_l = 2 * std::max((int)(avg_n / (average_cumulative_loop_size * prob_non_skip_update)), avg_n/2);
        //num_free_flips = N_l;

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        //ConfigurationGenerator::flipAllSpins();
        ConfigurationGenerator::randomSpinFlipsXXZ();

        /*
        if (!skip_loop_update){
            cumulative_loop_size_sum += cumulative_loop_size;
            cumulative_loop_size_list.push_back(cumulative_loop_size);
            avg_cumul_loop_size_d = ((double)cumulative_loop_size_sum)/((double)step+1);
            prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
            N_l_d = 10.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
            N_l = std::round(N_l_d);
            N_l = std::max(N_l,1);
            num_free_flips = N_l;
        }
         */

        cumulative_loop_size_sum += cumulative_loop_size;
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        avg_cumul_loop_size_d = ((double)cumulative_loop_size_sum)/((double)step+1);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
        N_l = std::round(N_l_d);
        N_l = std::max(N_l,1);
        num_free_flips = N_l;

        if (cumulative_loop_size_list.size() == avg_num_samples){
            break;
        }

        /*
        std::cout << "step = " << step << "\n";
        std::cout << "average n = " << avg_n_d << "\n";
        std::cout << "M = " << M << "\n";
        std::cout << "average cumulative loop size = " << avg_cumul_loop_size_d << "\n";
        std::cout << "N_l = " << N_l_d << "\n";
        std::cout << "Max loop size = " << max_loop_size << "\n";
         */

    }

    if (cumulative_loop_size_list.size() == 0){
        throw std::runtime_error("cumulative loop size list is empty\n");
    }

    avg_cumul_loop_size_d = ConfigurationGenerator::computeAverage(cumulative_loop_size_list,
                                                                   avg_num_samples);
    avg_n_d = ConfigurationGenerator::computeAverage(n_list, avg_num_samples);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/(avg_num_samples);

    avg_n = std::round(avg_n_d);
    average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
    //int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),M);
    //int M_new = std::max((int) (a_parameter * max_n), M);
    /*
    int M_new = (int)(a_parameter * avg_n_d);
    if (M_new - n <= 0){
        M_new = n+1;
    }
     */
    //int M_new = std::max(2*n,M);
    int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
    ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
    M = M_new;
    count_non_skipped_loop_updates = 0;

    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::round(N_l_d);
    N_l = std::max(N_l, 1);
    num_free_flips = N_l;

    for (int step=0; step<equilibration_steps; step++){
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);

        if (n > max_n) {
            max_n = n;
        }

        avg_n_d = ((avg_n_d * avg_num_samples) - n_list.at(n_list.size() - (int)avg_num_samples) + n)/avg_num_samples;
        avg_n = std::round(avg_n_d);
        n_list.push_back(n);

        //int M_new = std::max((int) (a_parameter * max_n), M);
        /*
        int M_new = (int)(a_parameter * avg_n_d);
        if (M_new - n <= 0){
            M_new = n+1;
        }
        */
        //int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),M);
        //int M_new = std::max(2*n,M);
        int M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M = M_new;

        max_loop_size = std::max(10000*avg_n,1);

        /*
        avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                 - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                                 + cumulative_loop_size)/avg_num_samples;
        average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        if (average_cumulative_loop_size < 0){
            std::cout << "Check" << "\n";
        }
         */

        //N_l = 2 * std::max((int)(avg_n / (average_cumulative_loop_size * prob_non_skip_update)), avg_n/2);
        //num_free_flips = N_l;

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        estimators.updateAllProperties(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset, num_winding);
        //ConfigurationGenerator::flipAllSpins();

        //estimators.updateAllProperties(n, spin_configuration, sim_params,
        //                               prob_tables.spectrum_offset, num_winding);
        ConfigurationGenerator::randomSpinFlipsXXZ();

        /*
        if (!skip_loop_update){
            avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                     - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                                     + cumulative_loop_size)/avg_num_samples;
            average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
            cumulative_loop_size_list.push_back(cumulative_loop_size);
            prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
            N_l_d = 10.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
            N_l = std::round(N_l_d);
            N_l = std::max(N_l, 1);
            num_free_flips = N_l;
        }
         */

        avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                 - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                                 + cumulative_loop_size)/avg_num_samples;
        average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
        N_l = std::round(N_l_d);
        N_l = std::max(N_l, 1);
        num_free_flips = N_l;

        if (step % 10000 == 0) {
            std::cout << "step = " << step << "\n";
            std::cout << "average n = " << avg_n_d << "\n";
            std::cout << "M = " << M << "\n";
            std::cout << "average cumulative loop size = " << avg_cumul_loop_size_d << "\n";
            std::cout << "Prob non-skipped update = " << prob_non_skip_update << "\n";
            std::cout << "N_l = " << N_l_d << "\n";
            std::cout << "max loop size = " << max_loop_size << "\n";
            for (int spin_index = 0; spin_index < N; spin_index++) {
                std::cout << spin_configuration.at(spin_index) << ",";
            }
            std::cout << "\n";
        }
    }

    max_loop_size = std::max(10000*avg_n,1);


    /*
    if (!skip_loop_update){
        avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                 - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                                 + cumulative_loop_size)/avg_num_samples;
        average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)equilibration_steps);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
        N_l = std::round(N_l_d);
        N_l = std::max(N_l, 1);
        num_free_flips = N_l;
    }
     */

    /*
    avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                             - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                             + cumulative_loop_size)/avg_num_samples;
    average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
    cumulative_loop_size_list.push_back(cumulative_loop_size);
     */

    avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                             - cumulative_loop_size_list.at(cumulative_loop_size_list.size() - (int)avg_num_samples)
                             + cumulative_loop_size)/avg_num_samples;
    average_cumulative_loop_size = std::round(avg_cumul_loop_size_d);
    cumulative_loop_size_list.push_back(cumulative_loop_size);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)equilibration_steps);
    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::round(N_l_d);
    N_l = std::max(N_l, 1);
    num_free_flips = N_l;

    //N_l = 2 * std::max((int)(avg_n / (average_cumulative_loop_size * prob_non_skip_update)), avg_n/2);
    //num_free_flips = N_l;

    std::cout << "Equilibration Finished" << "\n";
    std::cout << "average n = " << avg_n_d << "\n";
    std::cout << "M = " << M << "\n";
    std::cout << "average cumulative loop size = " << avg_cumul_loop_size_d << "\n";
    std::cout << "Prob non-skipped update = " << prob_non_skip_update << "\n";
    std::cout << "N_l = " << N_l_d << "\n";
    std::cout << "max loop size = " << max_loop_size << "\n";

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        /*
        int M_new = std::max(2*n,M);
        ConfigurationGenerator::populateOperatorLocations(int(M_new - M));
        M = M_new;
        */

        estimators.updateAllProperties(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset, num_winding);

        //ConfigurationGenerator::flipAllSpins();

        //estimators.updateAllProperties(n, spin_configuration, sim_params,
        //                               prob_tables.spectrum_offset, num_winding);

        ConfigurationGenerator::randomSpinFlipsXXZ();
    }

    estimators.outputStepData(sim_params);
}

double ConfigurationGenerator::computeAverage(std::vector<int> &vector_entries, const double &num_samples_in){

    double sum_entries = 0.0;
    int num_samples = (int)num_samples_in;
    for (int i=0; i<num_samples; i++){
        sum_entries += vector_entries.at(vector_entries.size()-i-1);
    }
    sum_entries = sum_entries/num_samples_in;
    double avg_val = sum_entries;

    return avg_val;
}


