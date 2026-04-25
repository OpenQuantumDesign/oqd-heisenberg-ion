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

Estimators::Estimators(const SimulationParameters &sim_params, const bool &track_spin_configs_in,
    std::shared_ptr<spdlog::logger> &logger_ptr) {

    logger = logger_ptr;
    logger->info("Initializing estimators");

    simulation_steps = sim_params.simulation_steps;
    track_spin_configs = track_spin_configs_in;
    out_folder_path = sim_params.simulation_subfolder + "/qmc_output/";

    if (not std::filesystem::is_directory(out_folder_path)){
        std::filesystem::create_directory(out_folder_path);
    }

    if (sim_params.track_spin_configs) {
        initializeShotStorage(sim_params);
        step_spin_configs.resize(chunk_size * sim_params.N);
    }
    else {
        chunk_size = 2 * simulation_steps;
    }
}

void Estimators::outputStepData() const {

    std::string file_path = out_folder_path + "/" + "estimators.csv";
    logger->info("Writing estimator data to: {}", file_path);
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

void Estimators::trackSpinConfigs(const std::vector<int> &spin_configs, const SimulationParameters &sim_params) {
    std::transform(spin_configs.begin(), spin_configs.end(),
        step_spin_configs.begin() + step_spin_config_offset,
        [](int i) { return static_cast<int8_t>(i);});
        step_spin_config_offset += sim_params.N;

    //step_spin_configs_diagnostics.push_back(spin_configs);
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

/*
void Estimators::outputDiagnostics(const SimulationParameters &sim_params) const {

    std::string file_path = out_folder_path + "/" + "spin_configurations_diagnostics.csv";
    logger->info("Writing shot data to: {}", file_path);
    std::string header = "# Spin Configurations\n";
    std::ofstream ofs(file_path);
    ofs << header << "\n";
    for (int i = 0; i < step_spin_configs_diagnostics.size(); i++){
        ofs << std::to_string(step_spin_configs_diagnostics.at(i).at(0));
        for (int j = 1; j < sim_params.N; j++){
            ofs << "," << std::to_string(step_spin_configs_diagnostics.at(i).at(j));
        }
        ofs << "\n";
    }
    ofs.close();
}
*/

void Estimators::initializeShotStorage(const SimulationParameters &sim_params) {

    std::string file_path = out_folder_path + "/" + "spin_configurations.h5";
    logger->info("Creating shot data file at: {}", file_path);
    shot_data_file_id = H5Fcreate(file_path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (shot_data_file_id == H5I_INVALID_HID) {
        logger->error("Shot data file creation failed. File path: {}", file_path);
        throw std::runtime_error("Shot data file creation failed.");
    }

    chunk_size = 250000/sim_params.N;
    chunk_size = std::max(chunk_size, uint64_t(1));
    chunk_size = std::min(chunk_size, static_cast<uint64_t>(sim_params.simulation_steps));
    dataset_create_property_list_id = H5Pcreate(H5P_DATASET_CREATE);
    if (dataset_create_property_list_id == H5I_INVALID_HID) {
        logger->error("Dataset property list creation failed");
        throw std::runtime_error("Dataset property list creation failed");
    }

    hsize_t chunk_dims[2] = {chunk_size, static_cast<uint64_t>(sim_params.N)};

    herr_t status = H5Pset_chunk(dataset_create_property_list_id, 2, chunk_dims);
    if (status < 0) {
        logger->error("H5Pset_chunk failed while initializing shot data storage");
        throw std::runtime_error("H5Pset_chunk failed");
    }

    hsize_t dims[2] = {static_cast<uint64_t>(simulation_steps), static_cast<uint64_t>(sim_params.N)};
    shot_data_dataspace_id = H5Screate_simple(2, dims, nullptr);
    if (shot_data_dataspace_id == H5I_INVALID_HID) {
        logger->error("Shot data space creation failed");
        throw std::runtime_error("Shot data space creation failed");
    }

    shot_data_dataset_id = H5Dcreate(shot_data_file_id, "/shot_data", H5T_STD_I8LE,
        shot_data_dataspace_id, H5P_DEFAULT,
        dataset_create_property_list_id, H5P_DEFAULT);
    if (shot_data_dataset_id == H5I_INVALID_HID) {
        logger->error("Shot dataset creation failed");
        throw std::runtime_error("Shot dataset creation failed");
    }
}

void Estimators::populateShotDataFile(const SimulationParameters &sim_params) {

    hsize_t subset_dims[2] = {chunk_size, static_cast<uint64_t>(sim_params.N)};

    hid_t memspace_id = H5Screate_simple(2, subset_dims, nullptr);

    hsize_t count[2] = {chunk_size, static_cast<uint64_t>(sim_params.N)};

    hsize_t offset[2] = {static_cast<uint64_t>(num_updates * chunk_size), 0};

    hsize_t stride[2] = {1,1};
    hsize_t block[2] = {1,1};

    herr_t status = H5Sselect_hyperslab(shot_data_dataspace_id, H5S_SELECT_SET,
        offset, stride, count, block);
    if (status < 0) {
        logger->error("H5Sselect_hyperslab failed while populating shot data file");
        throw std::runtime_error("H5Sselect_hyperslab failed");
    }

    status = H5Dwrite(shot_data_dataset_id, H5T_NATIVE_INT8, memspace_id,
        shot_data_dataspace_id, H5P_DEFAULT, step_spin_configs.data());
    if (status < 0) {
        logger->error("H5Dwrite failed while populating shot data file");
        throw std::runtime_error("H5Dwrite failed while populating shot data file");
    }

    num_updates++;
    step_spin_configs.clear();
    step_spin_config_offset = 0;

    status = H5Sclose(memspace_id);
    if (status < 0) {
        logger->error("H5Sclose failed while closing memspace_id after populating shot data file");
        throw std::runtime_error("H5Sclose failed while closing memspace_id after populating shot data file");
    }
}

void Estimators::closeShotDataContainers() const {



    herr_t status = H5Dclose(shot_data_dataset_id);
    if (status < 0) {
        logger->error("H5Dclose failed while closing shot_data_dataset_id after simulation");
        throw std::runtime_error("H5Dclose failed while closing shot_data_dataset_id after simulation");
    }

    status = H5Sclose(shot_data_dataspace_id);
    if (status < 0) {
        logger->error("H5Sclose failed while closing shot_data_dataspace_id after simulation");
        throw std::runtime_error("H5Sclose failed while closing shot_dataspace_id after simulation");
    }

    status = H5Pclose(dataset_create_property_list_id);
    if (status < 0) {
        logger->error("H5Pclose failed while closing dataset_create_property_list_id after simulation");
        throw std::runtime_error("H5Pclose failed while closing dataset_create_property_list_id after simulation");
    }

    status  = H5Fclose(shot_data_file_id);
    if (status < 0) {
        logger->error("H5Fclose failed while closing shot_data_file_id after simulation");
        throw std::runtime_error("H5Fclose failed while closing shot_data_file_id after simulation");
    }

}

