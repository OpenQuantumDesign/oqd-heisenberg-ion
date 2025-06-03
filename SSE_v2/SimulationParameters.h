#ifndef SSE_V2_SIMULATIONPARAMETERS_H
#define SSE_V2_SIMULATIONPARAMETERS_H
#include <string>
#include <iostream>
#include <fstream>

class SimulationParameters {

public:
    int N;
    double T;
    double J;
    double beta;
    double gamma;
    double alpha;
    double ksi;
    double Delta;
    double h;
    double h_B;
    int num_bonds;

    int num_chars_input_files;

    std::string file_prefix;
    std::string out_dir_subfolder;
    std::string root_folder;

    int simulation_steps;
    int equilibration_steps;
    double new_M_multiplier;

    std::string loop_type;

    int init_M;
    int init_config_index;

    SimulationParameters(const int &num_spins, const double &Delta_in, const double &h_in,
                         const double &alpha_in, const double &gamma_in, const double &ksi_in,
                         const double &J_in, const int &num_chars_in, const int &init_M_in,
                         const double &temperature_in, const std::string &loop_type_in,
                         const int &simulation_steps_in, const int &eq_steps_in,
                         const double &new_M_multiplier_in, const int &init_config_in,
                         const std::string &root_folder_in);

};


#endif //SSE_V2_SIMULATIONPARAMETERS_H
