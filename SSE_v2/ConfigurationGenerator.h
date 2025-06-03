#ifndef SSE_V2_CONFIGURATIONGENERATOR_H
#define SSE_V2_CONFIGURATIONGENERATOR_H
#include<vector>
#include<random>
#include "ProbabilityTables.h"
#include "Estimators.h"
#include "VertexTypes.h"

class ConfigurationGenerator {
private:
    std::vector<int> spin_configuration;
    std::vector<int> operator_locations;

    int num_legs_per_vertex = 4;
    int num_composite_leg_indices = 16;
    int num_total_legs;

    std::vector<int> p_list;
    std::vector<int> b_list;

    std::vector<int> first_vertex_leg;
    std::vector<int> last_vertex_leg;

    std::vector<int> vertex_configuration;
    std::vector<int> linked_list;

    double beta;

    // 1 => all spins set to 1, -1 => all spins set to -1, 0 => alternating 1s and 0s, 2 => all spins random
    int init_config_index;
    int init_M;

    int M;
    int n;
    int N_l;
    int max_loop_size;
    int N;

    std::vector<int> n_list;
    std::vector<int> cumulative_loop_size_list;
    int average_cumulative_loop_size;

    std::mt19937_64 initial_config_generator;
    std::mt19937_64 diagonal_update_generator;
    std::mt19937_64 off_diagonal_update_generator;
    std::mt19937_64 disconnected_spin_flip_generator;
    std::mt19937_64 loop_start_pos_generator;
    std::mt19937_64 metropolis_generator;

    int init_config_seed=1;
    int diagonal_update_seed=2;
    int off_diagonal_update_seed=3;
    int disconnected_spin_flip_seed=4;
    int loop_start_pos_seed=5;
    int metropolis_seed=6;

    int num_free_flips;
    int equilibration_steps;
    int mc_steps;
    double a_parameter;

    int cumulative_loop_size;

    std::vector<int> spin_labels;

    void setInitialConfiguration();

    void diagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void initializeOffDiagonalUpdates();

    void populateOperatorLocations(const int &num_fill_zeros);

public:
    explicit ConfigurationGenerator(const SimulationParameters &sim_params);

    void simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables, Estimators &estimators,
                                        const SimulationParameters &sim_params, const VertexTypes &vertex_types);
};


#endif //SSE_V2_CONFIGURATIONGENERATOR_H
