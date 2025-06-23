#ifndef SSE_V2_CONFIGURATIONGENERATOR_H
#define SSE_V2_CONFIGURATIONGENERATOR_H
#include<vector>
#include<unordered_set>
#include<random>
#include "ProbabilityTables.h"
#include "Estimators.h"
#include "VertexTypes.h"

class ConfigurationGenerator {
private:
    int num_legs_per_vertex = 4;
    int num_composite_leg_indices = 16;
    int num_total_legs;

    std::vector<int> p_list;
    std::vector<int> b_list;
    std::vector<int> id_list;
    std::vector<int> diag_op_p_list;
    std::vector<int> diag_order_op_p_list;
    std::vector<int> diag_disorder_op_p_list;
    std::vector<int> off_diag_op_p_list;
    std::vector<double> cluster_bond_prob_list;
    std::unordered_set<int> sliced_by_cluster;

    int count_p = 0;

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
    int num_clusters;
    int max_loop_size;
    int count_non_skipped_loop_updates;
    int N;
    int num_free_spins;
    int num_bonds;

    std::vector<int> n_list;
    std::vector<int> cumulative_loop_size_list;
    std::vector<int> cumulative_cluster_size_list;
    int average_cumulative_loop_size;

    int cluster_affected_spin_count;
    int count_non_skipped_cluster_updates;

    std::mt19937_64 initial_config_generator;
    std::mt19937_64 diagonal_update_generator;
    std::mt19937_64 off_diagonal_update_generator;
    std::mt19937_64 disconnected_spin_flip_generator;
    std::mt19937_64 loop_start_pos_generator;
    std::mt19937_64 metropolis_generator_1;
    std::mt19937_64 metropolis_generator_2;
    std::mt19937_64 metropolis_generator_3;
    std::mt19937_64 cluster_start_pos_generator;
    std::mt19937_64 metropolis_generator_4;

    int init_config_seed=1;
    int diagonal_update_seed=2;
    int off_diagonal_update_seed=3;
    int disconnected_spin_flip_seed=4;
    int loop_start_pos_seed=5;
    int metropolis_seed_1 =6;
    int metropolis_seed_2 = 7;
    int metropolis_seed_3 = 8;
    int cluster_start_pos_generator_seed = 9;
    int metropolis_seed_4 = 10;

    std::uniform_real_distribution<double> metropolis_acceptance_distribution = std::uniform_real_distribution<double>(0.0,1.0);
    std::uniform_int_distribution<> spin_flip_dist = std::uniform_int_distribution<>(0, 1);
    std::uniform_int_distribution<> loop_start_pos_dist;

    int num_free_flips;
    int equilibration_steps;
    int mc_steps;
    double a_parameter;

    int cumulative_loop_size;

    std::vector<int> spin_labels;

    std::vector<double> out_leg_probs;

    void setInitialConfiguration(const SimulationParameters &sim_params);

    void diagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void multiBranchClusterUpdate(const ProbabilityTables &prob_tables, const VertexTypes &vertexTypes);

    void populateLinkedList(const ProbabilityTables &prob_tables, const VertexTypes &vertex_types);

    void initializeOffDiagonalUpdates();

    void populateOperatorLocations(const int &num_fill_zeros);

    void randomSpinFlipsXXZ();

    void flipAllSpins();

    static double computeAverage(std::vector<int> &vector_entries, const double &num_samples_in);

public:
    std::vector<int> spin_configuration;
    std::vector<int> operator_locations;
    int num_winding;
    int M;
    int n;

    std::vector<int> cluster_spin_probs;

    ConfigurationGenerator(const SimulationParameters &sim_params, const ProbabilityTables &prob_tables);

    void simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables, Estimators &estimators,
                                        const SimulationParameters &sim_params, const VertexTypes &vertex_types);

    void simulateProbabilisticLoopsXXZ(const ProbabilityTables &prob_tables, Estimators &estimators,
                                        const SimulationParameters &sim_params, const VertexTypes &vertex_types);
};


#endif //SSE_V2_CONFIGURATIONGENERATOR_H
