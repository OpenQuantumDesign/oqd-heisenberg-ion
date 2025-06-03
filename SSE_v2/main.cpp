#include <iostream>
#include <chrono>
#include <filesystem>
#include <thread>
#include "VertexTypes.h"
#include "SimulationParameters.h"
#include "ConfigurationGenerator.h"
#include "ProbabilityTables.h"
#include "Estimators.h"

int main() {

    std::vector<double> h_list = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

    int input_N = 8;
    double input_Delta = 1.1;
    double input_alpha = 1.0;
    double input_gamma = 0.0;
    double input_ksi = 0.0;
    double input_J = 1.0;
    int input_num_chars = 3;
    int input_init_M = 10;
    int input_eq_steps = 10000;
    int input_mc_steps = 100000;
    double input_new_M_multiplier = 1.25;
    std::string input_loop_type = "directed_loops";
    std::string input_root_folder = "/Users/shaeermoeed/Github/Heisenberg_Ion/";
    double input_temperature = 0.5;
    int input_init_config_index = 0;

    for (int i = 0; i < h_list.size(); i++) {
        double input_h = h_list.at(i);

        SimulationParameters sim_params = SimulationParameters(input_N, input_Delta, input_h,
                                                               input_alpha,input_gamma, input_ksi,
                                                               input_J, input_num_chars,input_init_M,
                                                               input_temperature, input_loop_type,
                                                               input_mc_steps, input_eq_steps,
                                                               input_new_M_multiplier,
                                                               input_init_config_index,
                                                               input_root_folder);

        VertexTypes vertex_types = VertexTypes();

        ProbabilityTables prob_tables = ProbabilityTables(sim_params, vertex_types);

        Estimators estimators = Estimators(sim_params, false);

        ConfigurationGenerator mc_simulator = ConfigurationGenerator(sim_params);
        mc_simulator.simulateProbabilisticLoopsXXZh(prob_tables, estimators, sim_params, vertex_types);

    }
};



