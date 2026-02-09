#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include "VertexTypes.h"
#include "SimulationParameters.h"
#include "ConfigurationGenerator.h"
#include "ProbabilityTables.h"
#include "Estimators.h"
#include <spdlog/spdlog.h>
#include "spdlog/sinks/basic_file_sink.h"

void checkRootFolder(std::map<std::string, std::string> &input_key_vals, const auto &logger) {

    logger->info("Validating root directories");

    if (!input_key_vals.contains("simulation_folder") || input_key_vals["simulation_folder"].empty()) {
        logger->error("No root folder provided");
        logger->flush();
        throw std::runtime_error("No root folder provided\n");
    }
    else if(!std::filesystem::exists(input_key_vals["simulation_folder"])) {
        logger->error("Simulation folder path does not exist");
        logger->flush();
        throw std::runtime_error("Simulation folder path does not exist\n");
    }

    logger->info("Root directories validated");

}

void checkUUIDs(std::map<std::string, std::string> &input_key_vals, const int num_parameter_sets, const auto &logger) {

    logger->info("Checking UUID collection");

    if (!input_key_vals.contains("uuid") || input_key_vals["uuid"].empty()) {
        logger->error("No UUIDs provided");
        logger->flush();
        throw std::runtime_error("No UUID provided\n");
    }
    else {
        std::unordered_set<std::string> uuid_collection;
        std::string uuid_str = input_key_vals["uuid"];
        std::stringstream uuid_ss(uuid_str);
        std::string single_uuid;
        while(std::getline(uuid_ss, single_uuid, ',')) {
            if (uuid_collection.contains(single_uuid)) {
                logger->error("UUIDs are not unique");
                logger->flush();
                throw std::runtime_error("UUIDs are not unique\n");
            }
            else {
                uuid_collection.insert(single_uuid);
            }
        }
        size_t num_uuids = uuid_collection.size();
        if (num_uuids != num_parameter_sets) {
            logger->error("Number of UUIDs should be equal to the number of parameter sets provided");
            logger->flush();
            throw std::runtime_error("Number of UUIDs should be equal to the number of parameter sets provided\n");
        }
    }

    logger->info("UUID collection validated");
}

void checkNumberOfThreads(std::map<std::string, std::string> &input_key_vals, const auto &logger) {

    logger->info("Checking number of threads input");

    if (input_key_vals.contains("number_of_threads")) {
        std::string num_threads_str = input_key_vals["number_of_threads"];
        if (num_threads_str.empty() || num_threads_str == "-1") {
            input_key_vals["number_of_threads"] = "1";
        }
        else {
            size_t pos;
            std::stoi(num_threads_str, &pos);
            if (pos != num_threads_str.length()){
                logger->error("Value could not be converted to an integer for key: 'number_of_threads");
                logger->flush();
                throw std::runtime_error("Value could not be converted to an integer for key: 'number_of_threads'\n");
            }
        }
    }
    else {
        input_key_vals["number_of_threads"] = "1";
    }

    logger->info("Number of threads input checked");
}

void readInputFile(const std::string &input_file_path, std::map<std::string, std::string> &input_key_vals,
                   std::map<std::string, bool> &is_parameter_iterable, const auto &logger){

    logger->info("Reading input file: " + input_file_path);

    std::ifstream file_stream(input_file_path);

    if (!file_stream.good()) {
        logger->error("Could not open file with path: " + input_file_path);
        logger->flush();
        throw std::runtime_error("Could not open file with path: " + input_file_path);
    }

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
        }
        else {
            if (num_parameter_sets == 1){
                num_parameter_sets = count_entries;
                line_num_param_sets_identified = line_count;
            }
            else if (num_parameter_sets != count_entries) {
                logger->error("Inconsistent numbers of entries in line numbers: " +
                                         std::to_string(line_num_param_sets_identified) +
                                         " and " + std::to_string(line_count));
                logger->flush();
                throw std::runtime_error("Inconsistent numbers of entries in line numbers: " +
                                         std::to_string(line_num_param_sets_identified) +
                                         " and " + std::to_string(line_count) + "\n");
            }
            is_parameter_iterable[line_key] = true;
        }
    }

    input_key_vals["number_of_parameter_sets"] = std::to_string(num_parameter_sets);
    logger->info("Number of parameter sets: " + std::to_string(num_parameter_sets));

    file_stream.close();

    logger->info("Input file read");

    checkRootFolder(input_key_vals, logger);
    checkUUIDs(input_key_vals, num_parameter_sets, logger);
    checkNumberOfThreads(input_key_vals, logger);
}

void gatherThreadParameterSets(const int &num_sets_per_thread, std::map<std::string, bool> &specify_iterables,
                               std::map<std::string, std::string> &input_settings,
                               std::vector<std::map<std::string, std::string>> &thread_parameter_sets) {

    for (int j = 0; j < num_sets_per_thread; j++) {
        std::map<std::string, std::string> sim_parameters;
        for (const std::pair<const std::string, std::string>& key_value_pair : input_settings) {
            std::string key = key_value_pair.first;
            if (key == "number_of_threads") {
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
                         std::vector<std::vector<std::map<std::string, std::string>>> &parameter_sets,
                         const auto &logger) {

    logger->info("Gathering parameter sets");

    int num_total_parameter_sets = std::stoi(input_settings["number_of_parameter_sets"]);
    int num_threads = std::stoi(input_settings["number_of_threads"]);
    if (num_threads > num_total_parameter_sets) {
        num_threads = num_total_parameter_sets;
    }
    logger->info("Number of threads: " + std::to_string(num_threads));

    int num_sets_per_thread = num_total_parameter_sets / num_threads;
    logger->info("Number of threads per set: " + std::to_string(num_sets_per_thread));

    int extra_sets = num_total_parameter_sets - num_sets_per_thread * num_threads;
    logger->info("Number of extra sets: " + std::to_string(extra_sets));

    for (int i = 0; i < extra_sets; i ++) {
        std::vector<std::map<std::string, std::string>> thread_parameter_sets;
        gatherThreadParameterSets(num_sets_per_thread + 1, specify_iterables,
                                  input_settings,thread_parameter_sets);
        parameter_sets.push_back(thread_parameter_sets);
    }

    for (int i = 0; i < num_threads - extra_sets; i ++) {
        std::vector<std::map<std::string, std::string>> thread_parameter_sets;
        gatherThreadParameterSets(num_sets_per_thread, specify_iterables, input_settings,
                                  thread_parameter_sets);
        parameter_sets.push_back(thread_parameter_sets);
    }

    logger->info("Parameter sets gathered");
}

void simulateParameterSets(std::vector<std::map<std::string, std::string>> parameter_sets) {

    for (std::map<std::string, std::string>& simulation_specs : parameter_sets){

        std::string parameter_set_log_dir = simulation_specs["run_folder"] + "/parameter_set.log";
        std::string uuid_str = simulation_specs["uuid"];

        auto parameter_set_logger = spdlog::basic_logger_st(uuid_str,
            parameter_set_log_dir, true);

        parameter_set_logger->info("Parameter set logger initialized");

        SimulationParameters sim_params = SimulationParameters(simulation_specs, parameter_set_logger);

        parameter_set_logger->info("Simulation parameters initialized");

        parameter_set_logger->info("Initializing vertex types for SSE calculation");
        VertexTypes vertex_types = VertexTypes(sim_params.hamiltonian_type);
        parameter_set_logger->info("Vertex types initialized");

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types, parameter_set_logger);

        Estimators estimators = Estimators(sim_params, false, parameter_set_logger);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables, parameter_set_logger);
        mc_simulator.generateConfigurations(prob_tables, estimators, sim_params, vertex_types);

        parameter_set_logger->info("SSE calculation finished. Writing outputs.");

        estimators.outputStepData(sim_params);

        if (sim_params.track_spin_configs) {
            estimators.outputDiagnostics(sim_params);
        }

        if (sim_params.write_final_SSE_configs) {
            mc_simulator.writeFinalConfigurations(sim_params);
        }

        parameter_set_logger->info("SSE output files generated");
    }
}

int main(int arg_count, char *args[]) {

    if (arg_count < 2) {
        std::cerr << "Input file path not provided. Correct usage: " << args[0] << " <input_file_path>" << std::endl;
        return 1;
    }
    else {
        std::string input_file_path = args[1];

        std::filesystem::path file_path(args[1]);
        std::filesystem::path dir = file_path.parent_path();

        auto global_logger = spdlog::basic_logger_st("global_logger",
            dir/ "global.log", true);
        global_logger->flush_on(spdlog::level::info);

        global_logger->info("Global logger initialized");

        std::map<std::string, std::string> input_settings;
        std::map<std::string, bool> specify_iterables;

        readInputFile(input_file_path, input_settings, specify_iterables, global_logger);

        std::vector<std::vector<std::map<std::string, std::string>>> parameter_sets;
        std::vector<std::thread> threads_list;

        gatherParameterSets(specify_iterables, input_settings, parameter_sets, global_logger);

        global_logger->info("Executing threads");

        threads_list.reserve(parameter_sets.size());
        for (const std::vector<std::map<std::string, std::string>>& thread_parameter_sets : parameter_sets) {
            threads_list.emplace_back(simulateParameterSets, thread_parameter_sets);
        }

        for(std::thread& thread : threads_list){
            thread.join();
        }

        global_logger->info("All threads finished");

        // Fix the structure factor calculation and test it

        // Add a logger

        // Expose seeds for all generators

    }

    return 0;

}