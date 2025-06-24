#include "Estimators.h"

void Estimators::updateAllPropertiesProbabilistic(const int &int_step_n, const std::vector<int> &spin_configs,
                                     const SimulationParameters &sim_params, const double &spectrum_offset,
                                     const int &winding, const bool &skip_loop_update_step) {

    if (!skip_loop_update_step){
        step_energy.push_back(-int_step_n/sim_params.beta + spectrum_offset);
        double M_z = 0.0;
        for (int i = 0; i < sim_params.N; i++) {
        M_z += spin_configs.at(i);
        }
        M_z *= 0.5/(double)sim_params.N;
        step_magnetization.push_back(M_z);
        step_spin_stiffness.push_back(pow((double)winding,2) * 1.5/(sim_params.beta * sim_params.N));
    }
}

void Estimators::updateAllPropertiesDeterministic(const int &int_step_n, const std::vector<int> &spin_configs,
                                                  const SimulationParameters &sim_params, const double &spectrum_offset,
                                                  const int &winding) {

    step_energy.push_back(-int_step_n/sim_params.beta + spectrum_offset);
    double M_z = 0.0;
    for (int i = 0; i < sim_params.N; i++) {
        M_z += spin_configs.at(i);
    }
    M_z *= 0.5/(double)sim_params.N;
    step_magnetization.push_back(M_z);
    step_spin_stiffness.push_back(pow((double)winding,2) * 1.5/(sim_params.beta * sim_params.N));
}

Estimators::Estimators(const SimulationParameters &sim_params, const bool &track_spin_configs_in) {

    simulation_steps = sim_params.simulation_steps;
    track_spin_configs = track_spin_configs_in;
    out_folder_path = sim_params.root_folder + sim_params.out_dir_subfolder + "/";

    if (not std::filesystem::is_directory(out_folder_path)){
        std::filesystem::create_directory(out_folder_path);
    }
}

void Estimators::outputStepData(const SimulationParameters &sim_params) const {

    std::string filePath = out_folder_path + "/" + "MC Step Outputs.csv";
    std::string header = sim_params.file_prefix + "_" + sim_params.loop_type + "\n";
    header += "MC Step, Energy, Magnetization, Spin Stiffness\n";
    std::ofstream ofs(filePath);
    ofs << header;
    for (int i = 0; i < step_energy.size(); i++){
        ofs << std::to_string(i+1)
            << "," << std::to_string(step_energy.at(i))
            << "," << std::to_string(step_magnetization.at(i))
            << "," << std::to_string(step_spin_stiffness.at(i))
            << "\n";
    }
    ofs.close();

}

/*
void Estimators::outputClusterHistogram(const SimulationParameters &sim_params, const std::vector<int> &cluster_probs) const {

    std::string filePath = out_folder_path + "/" + "Cluster Histogram.csv";
    std::string header = sim_params.file_prefix + "_" + sim_params.loop_type + "\n";
    header += "Site, Cluster Probability\n";
    std::ofstream ofs(filePath);
    ofs << header;
    for (int i = 0; i < sim_params.N; i++){
        ofs << std::to_string(i+1)
            << "," << std::to_string(cluster_probs.at(i))
            << "\n";
    }
    ofs.close();
}
 */

void Estimators::outputDiagnostics(const SimulationParameters &sim_params) const {

    std::string filePath = out_folder_path + "/" + "MC Spin Configurations.csv";
    std::string header = sim_params.out_dir_subfolder + "\n";
    header += "Spin Configurations\n";
    std::ofstream ofs(filePath);
    ofs << header;
    for (int i = 0; i < step_energy.size(); i++){
        ofs << std::to_string(step_spin_configs.at(i).at(0));
        for (int j = 1; j < sim_params.N; j++){
            ofs << "," << std::to_string(step_spin_configs.at(i).at(j));
        }
        ofs << "\n";
    }
    ofs.close();
}

