#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include "VertexTypes.h"
#include "SimulationParameters.h"
#include "ConfigurationGenerator.h"
#include "ProbabilityTables.h"
#include "Estimators.h"

/*
void run_simulation(std::string input_h, std::string input_Delta, std::string input_gamma) {

    std::string input_N = "9";
    std::string input_alpha = "1.0";

    std::string input_ksi = "0.0";
    std::string input_J = "1.0";
    //int input_init_M = 10;
    int input_init_config_index;
    //input_init_config_index = 1;
    int input_eq_steps = 100000;
    int input_mc_steps = 100000;
    double input_new_M_multiplier = 1.1;
    std::string input_loop_type = "directed_loops";
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.1;
    bool dist_dep_offset = true;

    VertexTypes vertex_types = VertexTypes();

    std::cout << "Thread started" << "\n";
    if (std::stod(input_Delta) < 0.0) {
        for (int i=0; i<5; i++) {
            input_init_config_index = i;
            SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                                   input_alpha,input_gamma, dist_dep_offset,
                                                                   input_ksi,input_J,input_temperature,
                                                                   input_loop_type,input_mc_steps,
                                                                   input_eq_steps,input_new_M_multiplier,
                                                                   input_init_config_index,input_root_folder);

            ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

            Estimators estimators = Estimators(sim_params, false);

            ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);
            mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

            estimators.outputStepData(sim_params);
            estimators.outputClusterHistogram(sim_params, mc_simulator.cluster_spin_probs);
        }

    }
    else{
        input_init_config_index = 1;
        SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                               input_alpha,input_gamma, dist_dep_offset,
                                                               input_ksi,input_J,input_temperature,
                                                               input_loop_type,input_mc_steps,
                                                               input_eq_steps,input_new_M_multiplier,
                                                               input_init_config_index,input_root_folder);

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

        Estimators estimators = Estimators(sim_params, false);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);

        mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

        estimators.outputStepData(sim_params);
        estimators.outputClusterHistogram(sim_params, mc_simulator.cluster_spin_probs);
    }
}
*/

void run_simulation_delta_list(std::string input_h, std::vector<std::string> input_Delta_list, std::string input_gamma, int start_config) {

    std::string input_N = "9";
    std::string input_alpha = "1.0";

    std::string input_ksi = "0.0";
    std::string input_J = "1.0";
    //int input_init_M = 10;
    int input_init_config_index = start_config;
    //input_init_config_index = 1;
    int input_eq_steps = 500000;
    int input_mc_steps = 500000;
    double input_new_M_multiplier = 1.1;
    std::string input_loop_type = "directed_loops";
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.1;
    bool dist_dep_offset = true;

    std::vector<int> input_spin_config;
    for (int site=0; site<std::stoi(input_N); site++){
        input_spin_config.push_back(1);
    }

    int init_M_in = 0;
    int init_n_in = 0;
    std::vector<int> init_op_locs;
    int winding_in = 0;

    VertexTypes vertex_types = VertexTypes();

    std::cout << "Thread started" << "\n";
    for (int i=0; i<input_Delta_list.size(); i++) {

        std::string input_Delta = input_Delta_list.at(i);
        SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                               input_alpha,input_gamma, dist_dep_offset,
                                                               input_ksi,input_J,input_temperature,
                                                               input_loop_type,input_mc_steps,
                                                               input_eq_steps,input_new_M_multiplier,
                                                               input_init_config_index,input_root_folder,
                                                               input_spin_config, init_M_in, init_n_in,
                                                               init_op_locs, winding_in, start_config);

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

        Estimators estimators = Estimators(sim_params, false);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params, prob_tables);
        mc_simulator.simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);

        estimators.outputStepData(sim_params);
        estimators.outputClusterHistogram(sim_params, mc_simulator.cluster_spin_probs);

        input_spin_config = mc_simulator.spin_configuration;
        init_M_in = mc_simulator.M;
        init_n_in = mc_simulator.n;
        winding_in = mc_simulator.num_winding;
        init_op_locs = mc_simulator.operator_locations;
        input_init_config_index = -2;
    }
}

int main() {

    std::vector<std::string> h_list = {"0.0"};

    //std::vector<std::string> Delta_list = {"-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0"};
    //std::vector<std::string> Delta_list = {"-2.0"};
    //std::vector<std::string> Delta_list = {"1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0"};
    //std::vector<std::string> Delta_list = {"-2.0", "-1.0", "0.0"};
    //std::vector<std::string> Delta_list = {"-4.0", "-5.0", "-6.0", "-7.0", "-8.0", "-9.0", "-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0"};
    std::vector<std::string> Delta_list = {"-20.0","-19.0","-18.0","-17.0","-16.0"};
    std::vector<std::string> Delta_list1 = {"-15.0", "-14.0", "-13.0", "-12.0", "-11.0"};
    std::vector<std::string> Delta_list2 = {"-10.0", "-9.0", "-6.0", "-7.0", "-8.0"};
    std::vector<std::string> Delta_list3 = {"-5.0", "-4.0", "-3.0", "-2.0", "-1.0"};
    std::vector<std::string> Delta_list4 = {"0.0", "1.0", "2.0", "3.0", "4.0"};
    std::vector<std::string> Delta_list5 = {"5.0", "6.0", "7.0", "8.0", "9.0"};
    std::vector<std::string> Delta_list6 = {"10.0", "11.0", "12.0", "13.0", "14.0"};

    std::vector<std::string> Delta_list_all = {"-1.0", "-2.0", "-3.0", "-4.0", "-5.0", "-6.0", "-7.0", "-8.0",
                                               "-9.0", "-10.0", "-11.0", "-12.0", "-13.0", "-14.0", "-15.0", "-16.0",
                                               "-17.0", "-18.0", "-19.0", "-20.0"};
    //Delta_list.insert(Delta_list.end(), Delta_list2.begin(), Delta_list2.end());
    //Delta_list.insert(Delta_list.end(), Delta_list3.begin(), Delta_list3.end());
    //Delta_list.insert(Delta_list.end(), Delta_list4.begin(), Delta_list4.end());

    std::vector<std::thread> threads_list;

    for (int j=0; j<Delta_list.size(); j++){

        for (int i = 0; i < h_list.size(); i++) {
            std::string input_h = h_list.at(i);

            //run_simulation(input_h, input_Delta);
            std::string input_Delta = Delta_list.at(j);
            std::string input_Delta1 = Delta_list1.at(j);
            std::string input_Delta2 = Delta_list2.at(j);
            std::string input_Delta3 = Delta_list3.at(j);
            std::string input_Delta4 = Delta_list4.at(j);
            std::string input_Delta5 = Delta_list5.at(j);
            std::string input_Delta6 = Delta_list6.at(j);
            //threads_list.emplace_back(run_simulation_delta_list, input_h, Delta_list_all, "0.0", 0);
            //threads_list.emplace_back(run_simulation_delta_list, input_h, Delta_list_all, "0.0", 1);
            //threads_list.emplace_back(run_simulation, input_h, input_Delta2, "0.0");
            //threads_list.emplace_back(run_simulation, input_h, input_Delta3, "0.0");
            //threads_list.emplace_back(run_simulation, input_h, input_Delta4, "0.1");
            //threads_list.emplace_back(run_simulation, input_h, input_Delta5, "0.1");
            //threads_list.emplace_back(run_simulation, input_h, input_Delta6, "0.1");

            /*
            if (input_Delta1 == "-13.0") {
                run_simulation(input_h, input_Delta1, "0.1");
            }
             */
        }
    }

    for (int i=0; i<5; i++) {
        threads_list.emplace_back(run_simulation_delta_list, "0.0", Delta_list_all, "0.0", i);
    }


    for(int i = 0; i < threads_list.size(); i++){
        threads_list.at(i).join();
    }
};



