#include "SimulationParameters.h"

SimulationParameters::SimulationParameters(const int &num_spins, const double &Delta_in, const double &h_in,
                                           const double &alpha_in, const double &gamma_in, const double &ksi_in,
                                           const double &J_in, const int &num_chars_in, const int &init_M_in,
                                           const double &temperature_in, const std::string &loop_type_in,
                                           const int &simulation_steps_in, const int &eq_steps_in,
                                           const double &new_M_multiplier_in, const int &init_config_in,
                                           const std::string &root_folder_in) {

    N = num_spins;
    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    Delta = Delta_in;
    h = h_in;
    h_B = h/(J * ((double)N - 1.0));
    alpha = alpha_in;
    gamma = gamma_in;
    ksi = ksi_in;
    J = J_in;

    num_chars_input_files = num_chars_in;

    loop_type = loop_type_in;
    simulation_steps = simulation_steps_in;
    equilibration_steps = eq_steps_in;
    init_M = init_M_in;

    new_M_multiplier = new_M_multiplier_in;

    T = temperature_in;
    beta = J/T;

    file_prefix = "N_" + std::to_string(N) +"_Delta_" +
                  std::to_string(Delta).substr(0, num_chars_input_files) +
                  "_h_" + std::to_string(h).substr(0, num_chars_input_files) +
                  "_alpha_" + std::to_string(alpha).substr(0, num_chars_input_files) +
                  "_gamma_" + std::to_string(gamma).substr(0, num_chars_input_files) +
                  "_ksi_" + std::to_string(ksi).substr(0, num_chars_input_files) +
                  "_J_" + std::to_string(J).substr(0, num_chars_input_files);

    out_dir_subfolder = "/Results/SSE/" + file_prefix + "_" + loop_type;
    root_folder = root_folder_in;

    init_config_index = init_config_in;

}
