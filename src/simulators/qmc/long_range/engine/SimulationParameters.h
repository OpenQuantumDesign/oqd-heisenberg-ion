#ifndef SSE_V2_SIMULATIONPARAMETERS_H
#define SSE_V2_SIMULATIONPARAMETERS_H
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <map>

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
    std::string interaction_type;

    std::string simulation_subfolder;
    std::string root_folder;

    int simulation_steps;
    int equilibration_steps;
    double new_M_multiplier;

    std::string loop_type;
    bool distance_dep_offset;

    //int init_M;
    int init_config_index;
    int initial_init_config_index;

    int init_M;
    int init_n;

    std::string boundary_conditions;

    bool track_spin_configs;
    bool write_final_SSE_configs;

    std::string uuid;

    std::string init_config_file_path;

    std::vector<int> init_spin_config;
    std::vector<int> init_operator_locations;
    double winding;

    explicit SimulationParameters(std::map<std::string, std::string> &input_key_vals);

    static void extractIntegerEntry(const std::string &key_str, const std::string &val_str, int &member_var,
                             const bool &enforce_minimum, const int &min_val=0);

    static void extractDoubleEntry(const std::string &key_str, const std::string &val_str, double &member_var,
                            const bool &enforce_minimum, const double &min_val=0.0);

    static void extractBoolEntry(const std::string &key_str, const std::string &val_str, bool &member_var);

    static void extractStringEntry(const std::string &key_str, const std::string &val_str, std::string &member_var);

    static void extractListInts(const std::string &key_str, const std::string &val_str, std::vector<int> &member_var,
                         const int &list_size, const bool &enforce_minimum, const int &min_val=0);

    void extractHamiltonianType(const std::string &key_str, const std::string &val_str);

    void extractBoundaryConditions(const std::string &key_str, const std::string &val_str);

    void extractLoopType(const std::string &key_str, const std::string &val_str);

    void extractInteractionType(const std::string &key_str, const std::string &val_str);

    void extractInitialConditionsFromFile(std::string &file_path);

    static void writeNumericEntry(const std::string &key_str, const int &val, std::ofstream &file_stream);

    static void writeNumericEntry(const std::string &key_str, const double &val, std::ofstream &file_stream);

    static void writeStringEntry(const std::string &key_str, const std::string &val, std::ofstream &file_stream);

    static void writeBoolEntry(const std::string &key_str, const bool &val, std::ofstream &file_stream);

};


#endif //SSE_V2_SIMULATIONPARAMETERS_H
