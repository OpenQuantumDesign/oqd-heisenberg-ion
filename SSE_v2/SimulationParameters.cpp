#include "SimulationParameters.h"

SimulationParameters::SimulationParameters(const std::string &num_spins, const std::string &Delta_in, const std::string &h_in,
                                           const std::string &alpha_in, const std::string &gamma_in, const bool &dist_dep_offset,
                                           const std::string &ksi_in,const std::string &J_in, const double &temperature_in,
                                           const std::string &loop_type_in,const int &simulation_steps_in, const int &eq_steps_in,
                                           const double &new_M_multiplier_in, const int &init_config_in,
                                           const std::string &root_folder_in) {

    N = std::stoi(num_spins);
    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    Delta = std::stod(Delta_in);
    h = std::stod(h_in);
    h_B = h/(J * ((double)N - 1.0));
    alpha = std::stod(alpha_in);
    gamma = std::stod(gamma_in);
    ksi = std::stod(ksi_in);
    J = std::stod(J_in);

    loop_type = loop_type_in;
    simulation_steps = simulation_steps_in;
    equilibration_steps = eq_steps_in;
    //init_M = init_M_in;

    new_M_multiplier = new_M_multiplier_in;

    T = temperature_in;
    beta = J/T;

    std::string distance_dep_offset;
    if (dist_dep_offset){
        distance_dep_offset = "1";
    }
    else {
        distance_dep_offset = "0";
    }

    file_prefix = "N_" + num_spins +"_Delta_" + Delta_in + "_h_" + h_in + "_alpha_" + alpha_in + "_gamma_" +
            gamma_in + "_ksi_" + ksi_in + "_J_" + J_in + "_dist_dep_offset_" + distance_dep_offset;

    out_dir_subfolder = "/Results/SSE/" + file_prefix + "_" + loop_type + "_input_config_" + std::to_string(init_config_in);
    root_folder = root_folder_in;

    init_config_index = init_config_in;

}
