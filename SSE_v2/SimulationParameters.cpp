#include "SimulationParameters.h"

SimulationParameters::SimulationParameters(std::map<std::string, std::string> &input_key_vals) {

    extractStringEntry("Root Folder", input_key_vals["Root Folder"], root_folder);
    extractStringEntry("UUID", input_key_vals["UUID"], uuid);
    simulation_subfolder = root_folder + "/" + uuid + "/";

    std::string file_path = simulation_subfolder + "/" + "Simulation Specs.txt";
    std::ofstream ofs(file_path);

    writeStringEntry("Output Folder", simulation_subfolder, ofs);

    extractIntegerEntry("N", input_key_vals["N"], N, true, 1);
    writeNumericEntry("N", N, ofs);

    num_bonds = (int)(((double)N * ((double)N-1.0))/2.0);

    extractLoopType("Loop Type", input_key_vals["Loop Type"]);
    writeStringEntry("Loop Type", loop_type, ofs);

    extractHamiltonianType("Hamiltonian Type", input_key_vals["Hamiltonian Type"]);
    writeNumericEntry("Hamiltonian Type", hamiltonian_type, ofs);

    extractBoundaryConditions("Boundary Conditions", input_key_vals["Boundary Conditions"]);
    writeStringEntry("Boundary Conditions", boundary_conditions, ofs);

    if (loop_type == "directed_loops") {
        extractDoubleEntry("gamma", input_key_vals["gamma"], gamma,
                           true,0.0);
        writeNumericEntry("gamma", gamma, ofs);
        extractDoubleEntry("ksi", input_key_vals["ksi"], ksi,
                           true, 0.0);
        writeNumericEntry("ksi", ksi, ofs);
        extractBoolEntry("Distance Dependent Offset", input_key_vals["Distance Dependent Offset"],
                         distance_dep_offset);
        writeBoolEntry("Distance Dependent Offset", distance_dep_offset, ofs);
    }

    if (hamiltonian_type == 2 || hamiltonian_type == 3) {
        extractDoubleEntry("Delta", input_key_vals["Delta"], Delta, false);
    }
    else {
        Delta = (double)hamiltonian_type;
    }
    writeNumericEntry("Delta", Delta, ofs);

    if (hamiltonian_type == 3) {
        extractDoubleEntry("h", input_key_vals["h"], h,
                           true, 0.0);
        h_B = h/(J * ((double)N - 1.0));
    }
    else {
        h = 0.0;
        h_B = 0.0;
    }
    writeNumericEntry("h", h, ofs);

    extractInteractionType("Interaction Type", input_key_vals["Interaction Type"]);
    writeStringEntry("Interaction Type", interaction_type, ofs);
    if (interaction_type == "Power-Law") {
        extractDoubleEntry("alpha", input_key_vals["alpha"], alpha,
                           true, 0.0);
        writeNumericEntry("alpha", alpha, ofs);
    }
    extractDoubleEntry("J", input_key_vals["J"], J,
                       true, 0.0);
    writeNumericEntry("J", J, ofs);

    extractIntegerEntry("Simulation Steps", input_key_vals["Simulation Steps"],
                        simulation_steps, true, 1);
    writeNumericEntry("Simulation Steps", simulation_steps, ofs);

    extractIntegerEntry("Equilibration Steps", input_key_vals["Equilibration Steps"],
                        simulation_steps, true, 1);
    writeNumericEntry("Equilibration Steps", equilibration_steps, ofs);

    extractDoubleEntry("Operator List Update Multiplier",
                       input_key_vals["Operator List Update Multiplier"],
                       new_M_multiplier, true, 1.0);
    writeNumericEntry("Operator List Update Multiplier", new_M_multiplier, ofs);

    extractDoubleEntry("T", input_key_vals["T"], T, true, 0.0);
    writeNumericEntry("T", T, ofs);
    beta = J/T;

    extractBoolEntry("Track Spin Configurations", input_key_vals["Track Spin Configurations"],
                     track_spin_configs);
    writeBoolEntry("Track Spin Configurations", track_spin_configs, ofs);

    extractBoolEntry("Write Final Spin Configurations",
                     input_key_vals["Write Final Spin Configurations"], write_final_SSE_configs);
    writeBoolEntry("Write Final Spin Configurations", write_final_SSE_configs, ofs);

    extractIntegerEntry("Initial Configuration Index", input_key_vals["Initial Configuration Index"],
                        init_config_index, true, -2);
    writeNumericEntry("Initial Configuration Index", init_config_index, ofs);

    extractIntegerEntry("Initial Operator List Size", input_key_vals["Initial Configuration Index"],
                        init_M, true, -1);
    if (init_M == -1) { init_M = 50; }
    writeNumericEntry("Initial Operator List Size", init_M, ofs);

    if (init_config_index == -2) {
        extractStringEntry("Initial Configurations File Path",
                           input_key_vals["Initial Configuration Index"], init_config_file_path);
        extractInitialConditionsFromFile(init_config_file_path);
        writeStringEntry("Initial Configurations File Path", init_config_file_path, ofs);
    }

    ofs.close();
}

void SimulationParameters::extractIntegerEntry(const std::string &key_str, const std::string &val_str,
                                                      int &member_var, const bool &enforce_minimum,
                                                      const int &min_val) {

    if (val_str.empty()) {
        throw std::runtime_error("No value found for key: " + key_str + ".\n");
    }

    size_t pos;
    member_var = std::stoi(val_str, &pos);
    if (pos != val_str.length()){
        throw std::runtime_error("Value could not be converted to an integer for key: " + key_str + ".\n");
    }

    if (enforce_minimum && member_var < min_val){
        throw std::runtime_error("Value for key: " + key_str + " is below the expected minimum: " +
        std::to_string(min_val) + ".\n");
    }

}

void SimulationParameters::extractDoubleEntry(const std::string &key_str, const std::string &val_str,
                                              double &member_var, const bool &enforce_minimum,
                                              const double &min_val) {

    if (val_str.empty()) {
        throw std::runtime_error("No value found for key: " + key_str + ".\n");
    }

    size_t pos;
    member_var = std::stod(val_str, &pos);
    if (pos != val_str.length()){
        throw std::runtime_error("Value could not be converted to a double for key: " + key_str + ".\n");
    }

    if (enforce_minimum && member_var < min_val){
        throw std::runtime_error("Value for key: " + key_str + " is below the expected minimum: " +
                                 std::to_string(min_val) + ".\n");
    }
}

void SimulationParameters::extractBoolEntry(const std::string &key_str, const std::string &val_str,
                                            bool &member_var) {

    if (val_str.empty()) {
        throw std::runtime_error("No value found for key: " + key_str + ".\n");
    }

    if (val_str == "True") {
        member_var = true;
    }
    else if (val_str == "False") {
        member_var = false;
    }
    else {
        throw std::runtime_error("Value could not be converted to a bool for key: " + key_str + ".\n");
    }

}

void SimulationParameters::extractStringEntry(const std::string &key_str, const std::string &val_str,
                                              std::string &member_var) {

    if (val_str.empty()) {
        throw std::runtime_error("No value found for key: " + key_str + ".\n");
    }
    else {
        member_var = val_str;
    }

}

void SimulationParameters::extractListInts(const std::string &key_str, const std::string &val_str,
                                           std::vector<int> &member_var, const int &list_size,
                                           const bool &enforce_minimum, const int &min_val) {

    if (val_str.empty()) {
        throw std::runtime_error("No value found for key: " + key_str + ".\n");
    }

    std::istringstream list_stream(val_str);
    std::string list_val_str;
    int list_val;
    int entry_count = 0;
    while (std::getline(list_stream, list_val_str, ',')) {
        extractIntegerEntry(key_str + " Index " + std::to_string(entry_count), list_val_str,
                            list_val, enforce_minimum, min_val);
        entry_count++;
        member_var.push_back(list_val);
    }

    if (entry_count+1 != list_size) {
        throw std::runtime_error("Expected size for " + key_str + " does not match length of list\n");
    }

}

void SimulationParameters::extractLoopType(const std::string &key_str, const std::string &val_str) {

    extractStringEntry(key_str, val_str, loop_type);

    if (loop_type != "deterministic" && loop_type != "heat_bath" && loop_type != "directed_loops") {
        throw std::runtime_error("Loop Type must be 'deterministic', heat_bath or directed_loops\n");
    }
}

void SimulationParameters::extractHamiltonianType(const std::string &key_str, const std::string &val_str) {

    extractIntegerEntry(key_str, val_str, hamiltonian_type, false);

    if (hamiltonian_type > 3 || hamiltonian_type < -1) {
        throw std::runtime_error("Hamiltonian type can only be -1, 0, 1, 2 or 3\n");
    }

    if (loop_type == "deterministic") {
        if (hamiltonian_type != -1 && hamiltonian_type != 0 && hamiltonian_type != 1) {
            throw std::runtime_error("deterministic loops are only compatible with Hamiltonian types -1, 0 and 1\n");
        }
    }

    if (loop_type == "directed_loops" || loop_type == "heat_bath") {
        if (hamiltonian_type != 2 && hamiltonian_type != 3) {
            throw std::runtime_error("directed loops and heat bath probabilities are only compatible with Hamiltonian "
                                     "types 2 and 3\n");
        }
    }
}

void SimulationParameters::extractInteractionType(const std::string &key_str, const std::string &val_str) {

    extractStringEntry(key_str, val_str, interaction_type);

    if (interaction_type != "Power-Law" && interaction_type != "Matrix-Input") {
        throw std::runtime_error("Supported interaction types are: 'Power-Law' and 'Matrix-Input'\n");
    }
}

void SimulationParameters::extractBoundaryConditions(const std::string &key_str, const std::string &val_str) {

    extractStringEntry(key_str, val_str, boundary_conditions);

    if (boundary_conditions != "Open" && boundary_conditions != "Periodic") {
        throw std::runtime_error("Boundary conditions can either be 'Open' or 'Periodic'\n");
    }

}

void SimulationParameters::extractInitialConditionsFromFile(std::string &file_path) {

    std::ifstream file_stream(file_path);
    std::map<std::string, std::string> input_kev_val_pairs;

    if (!file_stream.good())
        throw std::runtime_error("Could not open file with path: " + file_path);

    std::string line;
    int line_count = 0;

    while (std::getline(file_stream, line)) {

        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        line_count += 1;

        std::stringstream row_stream(line);
        std::string line_key;
        std::string line_val;

        std::getline(row_stream, line_key, '\t');
        std::getline(row_stream, line_val, '\t');

        input_kev_val_pairs[line_key] = line_val;
    }

    file_stream.close();

    if (input_kev_val_pairs["N"] != std::to_string(N)) {
        throw std::runtime_error("Number of sites in SSE configuration inputs do not equal expected number of sites\n");
    }

    extractIntegerEntry("M", input_kev_val_pairs["M"], init_M, true, 1);
    extractIntegerEntry("n", input_kev_val_pairs["n"], init_n, true, 0);
    extractDoubleEntry("W", input_kev_val_pairs["W"], winding, false);

    extractListInts("Operator Locations List", input_kev_val_pairs["Operator Locations List"],
                    init_operator_locations,init_M, false);

    extractListInts("Spin Configurations List", input_kev_val_pairs["Spin Configurations List"],
                    init_spin_config,N, false);
}

void SimulationParameters::writeNumericEntry(const std::string &key_str, const int &val, std::ofstream &file_stream) {

    file_stream << key_str << "\t" << std::to_string(val) << "\n";

}

void SimulationParameters::writeNumericEntry(const std::string &key_str, const double &val, std::ofstream &file_stream) {

    file_stream << key_str << "\t" << std::to_string(val) << "\n";

}

void SimulationParameters::writeStringEntry(const std::string &key_str, const std::string &val,
                                            std::ofstream &file_stream) {

    file_stream << key_str << "\t" << val << "\n";

}

void SimulationParameters::writeBoolEntry(const std::string &key_str, const bool &val,
                                            std::ofstream &file_stream) {

    if (val) {
        file_stream << key_str << "\t" << "True" << "\n";
    }
    else {
        file_stream << key_str << "\t" << "False" << "\n";
    }
}