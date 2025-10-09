#include "SimulationParameters.h"

SimulationParameters::SimulationParameters(const std::string &num_spins, const std::string &Delta_in,
                                           const std::string &h_in, const std::string &alpha_in,
                                           const std::string &gamma_in, const bool &dist_dep_offset,
                                           const std::string &ksi_in,const std::string &J_in,
                                           const std::string &temperature_in, const std::string &loop_type_in,
                                           const int &simulation_steps_in, const int &eq_steps_in,
                                           const double &new_M_multiplier_in, const int &init_config_in,
                                           const std::string &root_folder_in, const std::vector<int> &input_spin_config,
                                           const int &init_M_in, const int &init_n_in, const std::vector<int> &init_op_locs,
                                           const double &winding_in, const int &initial_start_config_index,
                                           const int &hamiltonian_type_in, const std::string &boundary_conditions_in) {

    N = std::stoi(num_spins);
    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    Delta = std::stod(Delta_in);
    h = std::stod(h_in);
    h_B = h/(J * ((double)N - 1.0));

    gamma = std::stod(gamma_in);
    ksi = std::stod(ksi_in);

    setJAndAlpha(alpha_in, J_in);

    hamiltonian_type = hamiltonian_type_in;

    boundary_conditions = std::stoi(boundary_conditions_in);

    if (hamiltonian_type != 2) {
        throw std::runtime_error("This constructor is for probabilistic loops (hamiltonian type = 2)\n");
    }

    loop_type = loop_type_in;
    simulation_steps = simulation_steps_in;
    equilibration_steps = eq_steps_in;
    //init_M = init_M_in;

    new_M_multiplier = new_M_multiplier_in;

    T = std::stod(temperature_in);
    beta = J/T;

    if (dist_dep_offset){
        distance_dep_offset = "1";
    }
    else {
        distance_dep_offset = "0";
    }

    init_config_index = init_config_in;
    if (init_config_index != -2) {
        throw std::runtime_error("init_config_index should be -2 for this constructor\n");
    }
    SimulationParameters::setInitialConfigurationsFromInput(input_spin_config, init_M_in, init_n_in, init_op_locs,
                                                            winding_in, initial_start_config_index);

    file_prefix = "N_" + num_spins + "_hamiltonian_type_" + std::to_string(hamiltonian_type_in)
                  + "_Delta_" + Delta_in + "_h_" + h_in + "_alpha_" + alpha_in + "_gamma_" +
                  gamma_in + "_ksi_" + ksi_in + "_J_" + J_in + "_dist_dep_offset_" + distance_dep_offset +
                  "_boundary_" + boundary_conditions_in;

    out_dir_subfolder = "/Results/SSE/" + file_prefix + "_" + "T" + "_" + temperature_in + "_" + loop_type +
            "_input_config_" + std::to_string(init_config_in) +
            "_initial_input_config_" + std::to_string(initial_start_config_index);

    root_folder = root_folder_in;
}

SimulationParameters::SimulationParameters(const std::string &num_spins, const std::string &Delta_in,
                                           const std::string &h_in, const std::string &alpha_in,
                                           const std::string &gamma_in, const bool &dist_dep_offset,
                                           const std::string &ksi_in,const std::string &J_in,
                                           const std::string &temperature_in, const std::string &loop_type_in,
                                           const int &simulation_steps_in, const int &eq_steps_in,
                                           const double &new_M_multiplier_in, const int &init_config_in,
                                           const std::string &root_folder_in, const int &hamiltonian_type_in,
                                           const std::string &boundary_conditions_in) {

    N = std::stoi(num_spins);
    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    Delta = std::stod(Delta_in);
    h = std::stod(h_in);
    h_B = h/(J * ((double)N - 1.0));
    gamma = std::stod(gamma_in);
    ksi = std::stod(ksi_in);
    hamiltonian_type = hamiltonian_type_in;

    setJAndAlpha(alpha_in, J_in);

    boundary_conditions = std::stoi(boundary_conditions_in);

    if (hamiltonian_type != 2) {
        throw std::runtime_error("This constructor is for probabilistic loops (hamiltonian type = 2)\n");
    }

    loop_type = loop_type_in;
    simulation_steps = simulation_steps_in;
    equilibration_steps = eq_steps_in;
    //init_M = init_M_in;

    new_M_multiplier = new_M_multiplier_in;

    T = std::stod(temperature_in);
    beta = J/T;

    if (dist_dep_offset){
        distance_dep_offset = "1";
    }
    else {
        distance_dep_offset = "0";
    }

    init_config_index = init_config_in;
    if (init_config_index == -2) {
        throw std::runtime_error("init_config_index should not be -2 for this constructor\n");
    }
    initial_init_config_index = init_config_in;

    file_prefix = "N_" + num_spins + "_hamiltonian_type_" + std::to_string(hamiltonian_type_in)
                  + "_Delta_" + Delta_in + "_h_" + h_in + "_alpha_" + alpha_in + "_gamma_" +
                  gamma_in + "_ksi_" + ksi_in + "_J_" + J_in + "_dist_dep_offset_" + distance_dep_offset +
                  "_boundary_" + boundary_conditions_in;

    out_dir_subfolder = "/Results/SSE/" + file_prefix + "_" + "T" + "_" + temperature_in + "_" + loop_type +
                        "_input_config_" + std::to_string(init_config_in) +
                        "_initial_input_config_" + std::to_string(initial_init_config_index);

    root_folder = root_folder_in;

}

SimulationParameters::SimulationParameters(const std::string &num_spins, const std::string &alpha_in,
                                           const std::string &J_in, const std::string &temperature_in,
                                           const int &simulation_steps_in, const int &eq_steps_in,
                                           const int &init_config_in, const std::string &root_folder_in,
                                           const double &new_M_multiplier_in, const int &hamiltonian_type_in,
                                           const std::string &boundary_conditions_in) {

    N = std::stoi(num_spins);
    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    h = 0.0;
    h_B = 0.0;
    gamma = 0.0;
    ksi = 0.0;

    setJAndAlpha(alpha_in, J_in);

    hamiltonian_type = hamiltonian_type_in;

    boundary_conditions = std::stoi(boundary_conditions_in);

    std::string Delta_in;
    if (hamiltonian_type == 1) {
        Delta = 1.0;
        Delta_in = "1.0";
    }
    else if (hamiltonian_type == -1) {
        Delta = -1.0;
        Delta_in = "-1.0";
    }
    else if (hamiltonian_type == 0) {
        Delta = 0.0;
        Delta_in = "0.0";
    }
    else {
        throw std::runtime_error("This constructor requires the hamiltonian type to be either 0, -1 or 1\n");
    }

    new_M_multiplier = new_M_multiplier_in;

    loop_type = "deterministic";
    simulation_steps = simulation_steps_in;
    equilibration_steps = eq_steps_in;

    T = std::stod(temperature_in);
    beta = J/T;

    distance_dep_offset = "0";

    file_prefix = "N_" + num_spins + "_hamiltonian_type_" + std::to_string(hamiltonian_type_in)
            + "_Delta_" + Delta_in + "_h_0.0" + "_alpha_" + alpha_in +
            "_gamma_0.0_ksi_0.0" + "_J_" + J_in + "_dist_dep_offset_" + distance_dep_offset +
            "_boundary_" + boundary_conditions_in;

    init_config_index = init_config_in;
    if (init_config_index == -2) {
        throw std::runtime_error("init_config_index should not be -2 for this constructor\n");
    }
    initial_init_config_index = init_config_in;

    out_dir_subfolder = "/Results/SSE/" + file_prefix + "_" + "T" + "_" + temperature_in + "_" + loop_type +
                        "_input_config_" + std::to_string(init_config_in) +
                        "_initial_input_config_" + std::to_string(initial_init_config_index);

    root_folder = root_folder_in;
}

void SimulationParameters::setJAndAlpha(const std::string &alpha_in, const std::string &J_in){

    if (alpha_in.substr(0,3) == "exp"){
        if (J_in != "exp"){
            throw std::runtime_error("Either both J and alpha need to be exp or neither\n");
        }
        alpha = -1.0;
        J = 1.0;
    }
    else if (J_in == "exp"){
        if (alpha_in.substr(0,3) != "exp"){
            throw std::runtime_error("Either both J and alpha need to be exp or neither\n");
        }
        alpha = -1.0;
        J = 1.0;
    }
    else {
        J = std::stod(J_in);
        alpha = std::stod(alpha_in);
    }

}

void SimulationParameters::setInitialConfigurationsFromInput(const std::vector<int> &input_spin_config,
                                                             const int &init_M_in, const int &init_n_in,
                                                             const std::vector<int> &init_op_locs,
                                                             const double &winding_in,
                                                             const int &initial_start_config_index) {

    init_spin_config = input_spin_config;
    init_M = init_M_in;
    init_n = init_n_in;
    init_operator_locations = init_op_locs;
    winding = winding_in;
    initial_init_config_index = initial_start_config_index;

}
