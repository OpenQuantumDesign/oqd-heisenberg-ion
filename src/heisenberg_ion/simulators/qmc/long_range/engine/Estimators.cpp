#include "Estimators.h"

void Estimators::updateAllPropertiesProbabilistic(const int &int_step_n, const std::vector<int> &spin_configs,
                                     const SimulationParameters &sim_params, const double &spectrum_offset,
                                     const double &winding, const bool &skip_loop_update_step,
                                     const std::vector<std::vector<int>> &lattice_sites) {

    if (!skip_loop_update_step){
        step_energy.push_back(-int_step_n/sim_params.beta + spectrum_offset);
        double M_z = 0.0;
        for (int i = 0; i < sim_params.N; i++) {
        M_z += spin_configs.at(i);
        }
        M_z *= 0.5/(double)sim_params.N;
        step_magnetization.push_back(M_z);
        step_spin_stiffness.push_back(pow((double)winding,2)/(sim_params.beta * sim_params.N));
    }
}

void Estimators::updateAllPropertiesDeterministic(const int &int_step_n, const std::vector<int> &spin_configs,
                                                  const SimulationParameters &sim_params, const double &spectrum_offset,
                                                  const double &winding,
                                                  const std::vector<std::vector<int>> &lattice_sites) {

    step_energy.push_back(-int_step_n/sim_params.beta + spectrum_offset);
    double M_z = 0.0;
    for (int i = 0; i < sim_params.N; i++) {
        M_z += spin_configs.at(i);
    }
    M_z *= 0.5/(double)sim_params.N;
    step_magnetization.push_back(M_z);
    step_spin_stiffness.push_back(pow((double)winding,2)/(sim_params.beta * sim_params.N));
}

Estimators::Estimators(const SimulationParameters &sim_params, const bool &track_spin_configs_in) {

    simulation_steps = sim_params.simulation_steps;
    track_spin_configs = track_spin_configs_in;
    out_folder_path = sim_params.simulation_subfolder + "/qmc_output/";

    if (not std::filesystem::is_directory(out_folder_path)){
        std::filesystem::create_directory(out_folder_path);
    }
}

void Estimators::outputStepData(const SimulationParameters &sim_params) const {

    std::string file_path = out_folder_path + "/" + "estimators.csv";
    std::string header = "# MC Step Outputs\n";
    header += "MC Step, Energy, Magnetization, Spin Stiffness\n";
    std::ofstream ofs(file_path);
    ofs << header;
    for (int i = 0; i < step_energy.size(); i++){
        ofs << std::to_string(i+1)
            << "," << std::to_string(step_energy.at(i))
            << "," << std::to_string(step_magnetization.at(i))
            << "," << std::to_string(step_spin_stiffness.at(i)) << "\n";
    }
    ofs.close();

}

/*
void Estimators::outputPairCorrelations(const SimulationParameters &sim_params) const {

    std::string filePath = out_folder_path + "/" + "Pair Correlation Outputs.csv";
    std::string header = sim_params.file_prefix + "_" + sim_params.loop_type + "_" + "pair correlations" + "\n";
    header += "MC Step, bond";
    for (int i = 0; i < sim_params.num_bonds; i++) {
        header += std::to_string(i);
    }
    std::ofstream ofs(filePath);
    ofs << header << "\n";
    for (int i = 0; i < step_energy.size(); i++){
        ofs << std::to_string(i+1);
        for (int j = 0; j< sim_params.num_bonds; j++) {
            ofs << "," << step_pair_correlations.at(i).at(j);
        }
        ofs << "\n";
    }
    ofs.close();
}
 */

void Estimators::outputDiagnostics(const SimulationParameters &sim_params) const {

    std::string file_path = out_folder_path + "/" + "spin_configurations.csv";
    std::string header = "# Spin Configurations\n";
    std::ofstream ofs(file_path);
    ofs << header << "\n";
    for (int i = 0; i < step_spin_configs.size(); i++){
        ofs << std::to_string(step_spin_configs.at(i).at(0));
        for (int j = 1; j < sim_params.N; j++){
            ofs << "," << std::to_string(step_spin_configs.at(i).at(j));
        }
        ofs << "\n";
    }
    ofs.close();
}

