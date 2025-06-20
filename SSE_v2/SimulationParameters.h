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

    //int init_M;
    int init_config_index;

    SimulationParameters(const std::string &num_spins, const std::string &Delta_in, const std::string &h_in,
                         const std::string &alpha_in, const std::string &gamma_in, const bool &dist_dep_offset,
                         const std::string &ksi_in,const std::string &J_in, const double &temperature_in,
                         const std::string &loop_type_in,const int &simulation_steps_in, const int &eq_steps_in,
                         const double &new_M_multiplier_in, const int &init_config_in,
                         const std::string &root_folder_in);

};


#endif //SSE_V2_SIMULATIONPARAMETERS_H
