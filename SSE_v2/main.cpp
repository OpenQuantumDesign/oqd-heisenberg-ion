#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include "VertexTypes.h"
#include "SimulationParameters.h"
#include "ConfigurationGenerator.h"
#include "ProbabilityTables.h"
#include "Estimators.h"

void checkRootFolder(std::map<std::string, std::string> &input_key_vals) {

    if (!input_key_vals.contains("Root Folder") || input_key_vals["Root Folder"].empty()) {
        throw std::runtime_error("No root folder provided\n");
    }
    else if(!std::filesystem::exists(input_key_vals["Root Folder"])) {
        throw std::runtime_error("Root folder path does not exist\n");
    }

}

void checkUUIDs(std::map<std::string, std::string> &input_key_vals, const int num_parameter_sets) {

    if (!input_key_vals.contains("UUID") || input_key_vals["UUID"].empty()) {
        throw std::runtime_error("No UUID provided\n");
    }
    else {
        std::unordered_set<std::string> uuid_collection;
        std::string uuid_str = input_key_vals["UUID"];
        std::stringstream uuid_ss(uuid_str);
        std::string single_uuid;
        while(std::getline(uuid_ss, single_uuid, ',')) {
            if (uuid_collection.contains(single_uuid)) {
                throw std::runtime_error("UUIDs are not unique\n");
            }
            else {
                uuid_collection.insert(single_uuid);
            }
        }
        size_t num_uuids = uuid_collection.size();
        if (num_uuids != num_parameter_sets) {
            throw std::runtime_error("Number of UUIDs should be equal to the number of parameter sets provided\n");
        }
    }
}

void checkNumberOfThreads(std::map<std::string, std::string> &input_key_vals) {

    if (input_key_vals.contains("Number of Threads")) {
        std::string num_threads_str = input_key_vals["Number of Threads"];
        if (num_threads_str.empty()) {
            input_key_vals["Number of Threads"] = "1";
        }
        else {
            size_t pos;
            std::stoi(num_threads_str, &pos);
            if (pos != num_threads_str.length()){
                throw std::runtime_error("Value could not be converted to an integer for key: 'Number of Threads'\n");
            }
        }
    }
    else {
        input_key_vals["Number of Threads"] = "1";
    }
}

void readInputFile(const std::string &input_file_path, std::map<std::string, std::string> input_key_vals,
                   std::map<std::string, bool> is_parameter_iterable){

    std::ifstream file_stream(input_file_path);

    if (!file_stream.good())
        throw std::runtime_error("Could not open file with path: " + input_file_path);

    std::string line;
    int num_parameter_sets = 1;
    int line_num_param_sets_identified = 1;
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

        input_key_vals[line_key] = line_val;

        std::stringstream val_stream(line_val);
        std::string entry;
        int count_entries = 0;
        while (std::getline(val_stream, entry, ',')) {
            count_entries++;
        }

        if (count_entries == 1){
            is_parameter_iterable[line_key] = false;
            continue;
        }
        else {
            if (num_parameter_sets == 1){
                num_parameter_sets = count_entries;
                line_num_param_sets_identified = line_count;
            }
            else if (num_parameter_sets == count_entries){
                continue;
            }
            else {
                throw std::runtime_error("Inconsistent numbers of entries in line numbers: " +
                                         std::to_string(line_num_param_sets_identified) +
                                         " and " + std::to_string(line_count) + "\n");
            }
            is_parameter_iterable[line_key] = true;
        }
    }

    input_key_vals["Number of Parameter Sets"] = std::to_string(num_parameter_sets);

    file_stream.close();

    checkRootFolder(input_key_vals);
    checkUUIDs(input_key_vals, num_parameter_sets);
    checkNumberOfThreads(input_key_vals);
}

void gatherThreadParameterSets(const int &num_sets_per_thread, std::map<std::string, bool> &specify_iterables,
                               std::map<std::string, std::string> &input_settings,
                               std::vector<std::map<std::string, std::string>> &thread_parameter_sets) {

    for (int j = 0; j < num_sets_per_thread; j++) {
        std::map<std::string, std::string> sim_parameters;
        for (const std::pair<const std::string, std::string>& key_value_pair : input_settings) {
            std::string key = key_value_pair.first;
            if (key == "Number of Threads") {
                continue;
            }
            std::string val = key_value_pair.second;
            if (specify_iterables[key]) {
                std::stringstream val_stream(val);
                std::string entry;
                std::string remaining_string;
                std::getline(val_stream, entry, ',');
                std::getline(val_stream, remaining_string);
                sim_parameters[key] = entry;
                input_settings[key] = remaining_string;
            }
            else {
                sim_parameters[key] = val;
            }
        }
        thread_parameter_sets.push_back(sim_parameters);
    }
}

void gatherParameterSets(std::map<std::string, bool> &specify_iterables,
                         std::map<std::string, std::string> &input_settings,
                         std::vector<std::vector<std::map<std::string, std::string>>> &parameter_sets) {

    int num_total_parameter_sets = std::stoi(input_settings["Number of Parameter Sets"]);
    int num_threads = std::stoi(input_settings["Number of Threads"]);
    if (num_threads > num_total_parameter_sets) {
        num_threads = num_total_parameter_sets;
    }
    int num_sets_per_thread = num_total_parameter_sets / num_threads;
    int extra_sets = num_total_parameter_sets - num_sets_per_thread * num_threads;

    for (int i = 0; i < extra_sets; i ++) {
        std::vector<std::map<std::string, std::string>> thread_parameter_sets;
        gatherThreadParameterSets(num_sets_per_thread + 1, specify_iterables, input_settings,
                                  thread_parameter_sets);
        parameter_sets.push_back(thread_parameter_sets);
    }

    for (int i = 0; i < num_threads - extra_sets; i ++) {
        std::vector<std::map<std::string, std::string>> thread_parameter_sets;
        gatherThreadParameterSets(num_sets_per_thread, specify_iterables, input_settings,
                                  thread_parameter_sets);
        parameter_sets.push_back(thread_parameter_sets);
    }
}

void simulateParameterSets(std::vector<std::map<std::string, std::string>> parameter_sets) {

    for (std::map<std::string, std::string>& simulation_specs : parameter_sets){

        SimulationParameters sim_params = SimulationParameters(simulation_specs);

        VertexTypes vertex_types = VertexTypes(sim_params.hamiltonian_type);

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

        Estimators estimators = Estimators(sim_params, false);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);
        mc_simulator.generateConfigurations(prob_tables, estimators, sim_params, vertex_types);

        estimators.outputStepData(sim_params);

        if (sim_params.track_spin_configs) {
            estimators.outputDiagnostics(sim_params);
        }

        if (sim_params.write_final_SSE_configs) {
            mc_simulator.writeFinalConfigurations(sim_params);
        }
    }
}

int main(int arg_count, char *args[]) {

    if (arg_count < 2) {
        std::cerr << "Input file path not provided. Correct usage: " << args[0] << " <input_file_path>" << std::endl;
        return 1;
    }
    else {
        std::string input_file_path = args[1];
        std::map<std::string, std::string> input_settings;
        std::map<std::string, bool> specify_iterables;
        readInputFile(input_file_path, input_settings, specify_iterables);

        std::vector<std::vector<std::map<std::string, std::string>>> parameter_sets;
        std::vector<std::thread> threads_list;

        gatherParameterSets(specify_iterables, input_settings, parameter_sets);

        threads_list.reserve(parameter_sets.size());
        for (const std::vector<std::map<std::string, std::string>>& thread_parameter_sets : parameter_sets) {
            threads_list.emplace_back(simulateParameterSets, thread_parameter_sets);
        }

        for(std::thread& thread : threads_list){
            thread.join();
        }

        // Fix the structure factor calculation and test it

        // Add a logger

        // Clean up configuration generator

        // Expose seeds for all generators

    }

    return 0;

}