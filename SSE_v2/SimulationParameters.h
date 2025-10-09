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

    int hamiltonian_type;

    std::string file_prefix;
    std::string out_dir_subfolder;
    std::string root_folder;

    int simulation_steps;
    int equilibration_steps;
    double new_M_multiplier;

    std::string loop_type;
    std::string distance_dep_offset;

    //int init_M;
    int init_config_index;
    int initial_init_config_index;

    int init_M;
    int init_n;

    int boundary_conditions;

    SimulationParameters(const std::string &num_spins, const std::string &Delta_in, const std::string &h_in,
                         const std::string &alpha_in, const std::string &gamma_in, const bool &dist_dep_offset,
                         const std::string &ksi_in,const std::string &J_in, const std::string &temperature_in,
                         const std::string &loop_type_in,const int &simulation_steps_in, const int &eq_steps_in,
                         const double &new_M_multiplier_in, const int &init_config_in,
                         const std::string &root_folder_in, const std::vector<int> &input_spin_config,
                         const int &init_M, const int &init_n, const std::vector<int> &init_op_locs,
                         const double &winding_in,const int &initial_start_config_index, const int &hamiltonian_type_in,
                         const std::string &boundary_conditions_in);

    SimulationParameters(const std::string &num_spins, const std::string &Delta_in,
                         const std::string &h_in, const std::string &alpha_in,
                         const std::string &gamma_in, const bool &dist_dep_offset,
                         const std::string &ksi_in,const std::string &J_in,
                         const std::string &temperature_in, const std::string &loop_type_in,
                         const int &simulation_steps_in, const int &eq_steps_in,
                         const double &new_M_multiplier_in, const int &init_config_in,
                         const std::string &root_folder_in, const int &hamiltonian_type_in,
                         const std::string &boundary_conditions_in);

    SimulationParameters(const std::string &num_spins, const std::string &alpha_in,
                       const std::string &J_in, const std::string &temperature_in,
                       const int &simulation_steps_in, const int &eq_steps_in,
                       const int &init_config_in, const std::string &root_folder_in,
                       const double &new_M_multiplier_in, const int &hamiltonian_type_in,
                       const std::string &boundary_conditions_in);

    void setInitialConfigurationsFromInput(const std::vector<int> &input_spin_config, const int &init_M_in,
                                      const int &init_n_in, const std::vector<int> &init_op_locs,
                                      const double &winding_in, const int &initial_start_config_index);

    void setJAndAlpha(const std::string &alpha_in, const std::string &J_in);

    std::vector<int> init_spin_config;
    std::vector<int> init_operator_locations;
    double winding;

};


#endif //SSE_V2_SIMULATIONPARAMETERS_H
