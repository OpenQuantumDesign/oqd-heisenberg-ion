#ifndef cpp_qmc_ESTIMATORS_H
#define cpp_qmc_ESTIMATORS_H
#include <hdf5.h>
#include<vector>
#include<filesystem>
#include "SimulationParameters.h"

class Estimators {
private:
    //std::vector<double> step_correlation;
    std::vector<double> step_magnetization;
    std::vector<double> step_energy;

    std::shared_ptr<spdlog::logger> logger;

    std::string out_folder_path;

    int simulation_steps;

    std::vector<int8_t> step_spin_configs;
    //std::vector<std::vector<int>> step_spin_configs_diagnostics;
    int step_spin_config_offset = 0;
    bool track_spin_configs;

    std::vector<double> step_spin_stiffness;

    std::vector<double> step_S_k_pi;

    hid_t shot_data_file_id, shot_data_dataset_id, shot_data_dataspace_id, dataset_create_property_list_id;
    int num_updates = 0;

    void initializeShotStorage(const SimulationParameters &sim_params);

public:

    uint64_t chunk_size;

    Estimators(const SimulationParameters &sim_params, const bool &track_spin_configs,
        std::shared_ptr<spdlog::logger> &logger_ptr);

    void updateAllPropertiesProbabilistic(const int &i_step_n, const std::vector<int> &spin_configs,
                             const SimulationParameters &sim_params, const double &spectrum_offset,
                             const double &winding, const bool &skip_loop_update_step,
                             const std::vector<std::vector<int>> &lattice_sites);

    void updateAllPropertiesDeterministic(const int &i_step_n, const std::vector<int> &spin_configs,
                                          const SimulationParameters &sim_params, const double &spectrum_offset,
                                          const double &winding, const std::vector<std::vector<int>> &lattice_sites);

    void trackSpinConfigs(const std::vector<int> &spin_configs, const SimulationParameters &sim_params);

    void outputStepData() const;

    //void outputDiagnostics(const SimulationParameters &sim_params) const;

    void populateShotDataFile(const SimulationParameters &sim_params);

    //void outputPairCorrelations(const SimulationParameters &sim_params) const;

    void closeShotDataContainers() const;

    //void outputClusterHistogram(const SimulationParameters &sim_params, const std::vector<int> &cluster_probs) const;

};


#endif //cpp_qmc_ESTIMATORS_H
