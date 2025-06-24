#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include "VertexTypes.h"
#include "SimulationParameters.h"
#include "ConfigurationGenerator.h"
#include "ProbabilityTables.h"
#include "Estimators.h"

void run_simulation_XXZ(std::string input_h, std::string input_Delta, std::string input_gamma, int start_index) {

    std::string input_N = "3";
    std::string input_alpha = "1.0";

    std::string input_ksi = "0.0";
    std::string input_J = "1.0";
    int input_init_config_index = start_index;
    int input_eq_steps = 500000;
    int input_mc_steps = 500000;
    double input_new_M_multiplier = 1.1;
    std::string input_loop_type = "directed_loops";
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.1;
    bool dist_dep_offset = true;

    int hamiltonian_type = 2;
    VertexTypes vertex_types = VertexTypes(hamiltonian_type);

    std::cout << "Thread started" << "\n";

    SimulationParameters sim_params = SimulationParameters(input_N,input_Delta,
                                       input_h,input_alpha,input_gamma,
                                       dist_dep_offset,input_ksi,input_J,
                                       input_temperature,input_loop_type,
                                       input_mc_steps,input_eq_steps,
                                       input_new_M_multiplier,
                                       input_init_config_index,input_root_folder,
                                       hamiltonian_type);

    ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

    Estimators estimators = Estimators(sim_params, false);

    ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);
    mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

    estimators.outputStepData(sim_params);
}

void run_simulation_XXZ_output_start(std::string input_h, std::vector<std::string> input_Delta_list, std::string input_gamma,
                               int start_config) {

    std::string input_N = "3";
    std::string input_alpha = "1.0";

    std::string input_ksi = "0.0";
    std::string input_J = "1.0";
    //int input_init_M = 10;
    int input_init_config_index = start_config;
    //input_init_config_index = 1;
    //int input_eq_steps = 500000;
    //int input_mc_steps = 500000;
    int input_eq_steps = 500000;
    int input_mc_steps = 500000;
    double input_new_M_multiplier = 1.1;
    std::string input_loop_type = "directed_loops";
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.1;
    bool dist_dep_offset = true;
    std::string input_Delta = "-1.0";
    int hamiltonian_type_in = 2;

    std::cout << "Thread started" << "\n";

    VertexTypes vertex_types = VertexTypes(hamiltonian_type_in);

    SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                           input_alpha,input_gamma, dist_dep_offset,
                                                           input_ksi,input_J,input_temperature,
                                                           input_loop_type,input_mc_steps,
                                                           input_eq_steps,input_new_M_multiplier,
                                                           input_init_config_index,input_root_folder,
                                                           hamiltonian_type_in);

    ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

    Estimators estimators = Estimators(sim_params, false);

    ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);

    mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

    estimators.outputStepData(sim_params);

    std::vector<int> input_spin_config = mc_simulator.spin_configuration;
    int init_M_in = mc_simulator.M;
    int init_n_in = mc_simulator.n;
    int winding_in = mc_simulator.num_winding;
    std::vector<int> init_op_locs = mc_simulator.operator_locations;
    input_init_config_index = -2;

    for (int i=1; i<input_Delta_list.size(); i++) {

        std::string input_Delta = input_Delta_list.at(i);
        SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                               input_alpha,input_gamma,
                                                               dist_dep_offset,input_ksi,input_J,
                                                               input_temperature,
                                                               input_loop_type,
                                                               input_mc_steps,
                                                               input_eq_steps,
                                                               input_new_M_multiplier,
                                                               input_init_config_index,
                                                               input_root_folder, input_spin_config,
                                                               init_M_in, init_n_in,
                                                               init_op_locs, winding_in,
                                                               start_config,
                                                               hamiltonian_type_in);

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

        Estimators estimators = Estimators(sim_params, false);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);
        mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

        estimators.outputStepData(sim_params);

        input_spin_config = mc_simulator.spin_configuration;
        init_M_in = mc_simulator.M;
        init_n_in = mc_simulator.n;
        winding_in = mc_simulator.num_winding;
        init_op_locs = mc_simulator.operator_locations;
    }
}

void run_simulation_edge_cases(std::string input_alpha, int hamiltonian_type) {

    std::string input_N = "3";
    std::string input_J = "1.0";
    int input_eq_steps = 500000;
    int input_mc_steps = 500000;
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.1;

    int input_init_config_index;
    if (hamiltonian_type == -1) {
        input_init_config_index = 3;
    }
    else {
        input_init_config_index = 1;
    }

    double new_M_multiplier = 1.25;

    VertexTypes vertex_types = VertexTypes(hamiltonian_type);

    SimulationParameters sim_params = SimulationParameters(input_N, input_alpha,input_J,
                                                           input_temperature,
                                                           input_mc_steps,
                                                           input_eq_steps,
                                                           input_init_config_index,
                                                           input_root_folder,
                                                           new_M_multiplier,
                                                           hamiltonian_type);

    ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

    Estimators estimators = Estimators(sim_params, false);

    ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);

    std::cout << "thread started" << "\n";
    if (hamiltonian_type == 0) {
        mc_simulator.simulateDeterministicXY(prob_tables, estimators, sim_params, vertex_types);
    }
    else {
        mc_simulator.simulateDeterministicIsotropic(prob_tables, estimators, sim_params, vertex_types);
    }

    estimators.outputStepData(sim_params);
}

int main_XXZ_previous_config_start() {

    /*
    std::vector<std::string> Delta_list_all = {"-1.0", "-2.0", "-3.0", "-4.0", "-5.0", "-6.0", "-7.0", "-8.0",
                                               "-9.0", "-10.0", "-11.0", "-12.0", "-13.0", "-14.0", "-15.0", "-16.0",
                                               "-17.0", "-18.0", "-19.0", "-20.0"};
    */
    std::vector<std::string> Delta_list_all = {"-1.0", "1.0"};

    std::vector<std::thread> threads_list;

    //int num_start_configs = 5;
    int num_start_configs = 1;

    for (int i=0; i<num_start_configs; i++) {
        threads_list.emplace_back(run_simulation_XXZ_output_start, "0.0", Delta_list_all, "0.0", i);
    }


    for(int i = 0; i < threads_list.size(); i++){
        threads_list.at(i).join();
    }

    return 0;
}

int main_XXZ() {


    std::vector<std::thread> threads_list;

    //int num_start_configs = 5;
    int num_start_configs = 1;

    threads_list.emplace_back(run_simulation_XXZ, "0.0", "-1.0", "0.0", 2);
    //threads_list.emplace_back(run_simulation_XXZ, "0.0", "1.0", "0.0", 1);
    //threads_list.emplace_back(run_simulation_XXZ, "0.0", "2.0", "0.0", 1);
    //threads_list.emplace_back(run_simulation_XXZ, "0.0", "3.0", "0.0", 1);
    //threads_list.emplace_back(run_simulation_XXZ, "0.0", "-2.0", "0.0", 2);


    for(int i = 0; i < threads_list.size(); i++){
        threads_list.at(i).join();
    }

    return 0;

}

int main() {

    //std::vector<std::string> alpha_list = {"0.0","1.0","2.0","3.0","4.0","5.0"};
    std::vector<std::string> alpha_list = {"9.0"};

    std::vector<std::thread> threads_list;

    for (int j=0; j<alpha_list.size(); j++) {
        std::string alpha = alpha_list.at(j);
        threads_list.emplace_back(run_simulation_edge_cases, alpha, 0);
        threads_list.emplace_back(run_simulation_edge_cases, alpha, 1);
        threads_list.emplace_back(run_simulation_edge_cases, alpha, -1);
    }

    for(int i = 0; i < threads_list.size(); i++){
        threads_list.at(i).join();
    }

    return 0;
}



