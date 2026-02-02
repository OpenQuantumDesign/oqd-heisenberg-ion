#ifndef cpp_qmc_ESTIMATORS_H
#define cpp_qmc_ESTIMATORS_H
#include<vector>
#include "SimulationParameters.h"

class Estimators {
private:
    //std::vector<double> step_correlation;
    std::vector<double> step_magnetization;
    std::vector<double> step_energy;

    std::string out_folder_path;

    int simulation_steps;

    std::vector<std::vector<int>> step_spin_configs;
    bool track_spin_configs;

    std::vector<double> step_spin_stiffness;

    std::vector<double> step_S_k_pi;

public:
    Estimators(const SimulationParameters &sim_params, const bool &track_spin_configs);

    void updateAllPropertiesProbabilistic(const int &i_step_n, const std::vector<int> &spin_configs,
                             const SimulationParameters &sim_params, const double &spectrum_offset,
                             const double &winding, const bool &skip_loop_update_step,
                             const std::vector<std::vector<int>> &lattice_sites);

    void updateAllPropertiesDeterministic(const int &i_step_n, const std::vector<int> &spin_configs,
                                          const SimulationParameters &sim_params, const double &spectrum_offset,
                                          const double &winding, const std::vector<std::vector<int>> &lattice_sites);

    void trackSpinConfigs(const std::vector<int> &spin_configs) {
        step_spin_configs.push_back(spin_configs);
    }

    void outputStepData(const SimulationParameters &sim_params) const;

    void outputDiagnostics(const SimulationParameters &sim_params) const;

    void outputPairCorrelations(const SimulationParameters &sim_params) const;

    //void outputClusterHistogram(const SimulationParameters &sim_params, const std::vector<int> &cluster_probs) const;

};


#endif //cpp_qmc_ESTIMATORS_H
