#ifndef cpp_qmc_CONFIGURATIONGENERATOR_H
#define cpp_qmc_CONFIGURATIONGENERATOR_H
#include<vector>
#include<unordered_set>
#include<random>
#include "ProbabilityTables.h"
#include "Estimators.h"
#include "VertexTypes.h"
#include <spdlog/spdlog.h>

class ConfigurationGenerator {
private:
    int num_legs_per_vertex = 4;
    int num_composite_leg_indices = 16;
    int num_total_legs;

    std::shared_ptr<spdlog::logger> logger;

    int hamiltonian_type;

    std::vector<int> p_list;
    std::vector<int> b_list;

    std::unordered_set<int> connected_legs;

    std::vector<int> first_vertex_leg;
    std::vector<int> last_vertex_leg;

    std::vector<int> vertex_configuration;
    std::vector<int> linked_list;

    bool skip_loop_update;

    double beta;

    // 1 => all spins set to 1, -1 => all spins set to -1, 0 => alternating 1s and 0s, 2 => all spins random
    int init_config_index;
    int init_M;

    int N_l;
    int max_loop_size;
    int count_non_skipped_loop_updates;
    int N;
    int N_over_2;
    int num_free_spins;
    int num_bonds;

    std::vector<int> n_list;
    std::vector<int> cumulative_loop_size_list;
    std::vector<int> cumulative_cluster_size_list;
    int average_cumulative_loop_size;

    std::mt19937_64 initial_config_generator;
    std::mt19937_64 diagonal_update_generator;
    std::mt19937_64 exit_leg_generator;
    std::mt19937_64 disconnected_spin_flip_generator;
    std::mt19937_64 off_diagonal_spin_flip_generator;
    std::mt19937_64 loop_start_position_generator;
    std::mt19937_64 metropolis_insert_generator;
    std::mt19937_64 metropolis_bond_generator;
    std::mt19937_64 metropolis_remove_generator;
    //std::mt19937_64 twist_generator;

    int initial_config_seed;
    int diagonal_update_seed;
    int exit_leg_seed;
    int disconnected_spin_flip_seed;
    int off_diagonal_update_seed;
    int metropolis_insert_seed;
    int metropolis_bond_generator_seed;
    int metropolis_remove_seed;
    //int twist_seed = 9;

    std::uniform_real_distribution<double> metropolis_acceptance_distribution =
            std::uniform_real_distribution<double>(0.0,1.0);

    std::uniform_int_distribution<> binary_dist = std::uniform_int_distribution<>(0, 1);
    std::uniform_int_distribution<> loop_start_pos_dist;

    int equilibration_steps;
    int mc_steps;
    double a_parameter;

    int cumulative_loop_size;

    std::vector<int> spin_labels;

    std::vector<double> out_leg_probs;

    void setInitialConfiguration(const SimulationParameters &sim_params);

    void diagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void diagonalUpdatesXY(const ProbabilityTables &prob_tables);

    void diagonalUpdatesIsotropic(const ProbabilityTables &prob_tables);

    void offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void offDiagonalUpdatesXY(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void offDiagonalUpdatesIsotropic(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void populateLinkedList(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void getNewOperatorLocations(const int &num_fill_zeros);

    void mapVerticesToOperatorLocations(const VertexTypes &vertex_types,
                                        const ProbabilityTables &prob_tables);

    void initializeOffDiagonalUpdates();

    void populateOperatorLocations(const int &num_fill_zeros);

    void randomSpinFlipsXXZ();

    void flipAllSpins();

    void flipFreeSpin(const int &i);

    static double computeAverage(std::vector<int> &vector_entries, const double &num_samples_in);

    void simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables, Estimators &estimators,
                                        const SimulationParameters &sim_params, const VertexTypes &vertex_types);

    void simulateDeterministicIsotropic(const ProbabilityTables &prob_tables, Estimators &estimators,
                                        const SimulationParameters &sim_params, const VertexTypes &vertex_types);

    void simulateDeterministicXY(const ProbabilityTables &prob_tables, Estimators &estimators,
                                 const SimulationParameters &sim_params, const VertexTypes &vertex_types);

    void simulateProbabilisticLoopsXXZ(const ProbabilityTables &prob_tables, Estimators &estimators,
                                       const SimulationParameters &sim_params, const VertexTypes &vertex_types);

public:
    std::vector<int> spin_configuration;
    std::vector<int> operator_locations;
    double num_winding;
    int M;
    int n;

    ConfigurationGenerator(const SimulationParameters &sim_params, const ProbabilityTables &prob_tables,
        std::shared_ptr<spdlog::logger> &logger_ptr);

    void generateConfigurations(const ProbabilityTables &prob_tables, Estimators &estimators,
                                const SimulationParameters &sim_params, const VertexTypes &vertex_types);

    void writeFinalConfigurations(const SimulationParameters &sim_params);
};


#endif //cpp_qmc_CONFIGURATIONGENERATOR_H
