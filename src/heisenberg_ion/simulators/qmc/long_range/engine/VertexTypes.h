#ifndef cpp_qmc_VERTEXTYPES_H
#define cpp_qmc_VERTEXTYPES_H
#include<vector>
#include<map>
#include <spdlog/spdlog.h>

class VertexTypes {
private:
    std::map<int, int> config_to_index_mapping;
    std::vector<std::vector<int>> index_to_config_mapping;

    void populateAllowedExitLegs();
    void populateAllowedExitLegsXY();

    std::vector<int> flipInputOutputLegs(const int &l_e, const int &l_x, const std::vector<int> &old_config) const;

public:

    int num_vertices;
    int num_legs_per_vertex;
    int num_composite_leg_indices;

    VertexTypes(const int &hamiltonian_type);

    void setVertexMappings();
    void setVertexMappingsIsotropicAFM();
    void setVertexMappingsIsotropicFM();

    std::vector<int> allowed_exit_legs;
    std::vector<int> flip_left_half_vertex_map;
    std::vector<int> flip_right_half_vertex_map;
    std::vector<int> flip_full_vertex_map;

    std::vector<int> is_off_diag;
    std::vector<int> twist_mapping;

    int getVertexTypeIndex(const std::vector<int> &config) const;

    std::vector<int> getVertexConfig(const int &vertex_index) const {
        return index_to_config_mapping.at(vertex_index);
    };

    int getFlippedSpinsVertexIndex(const int &l_e, const int &l_x, const int &old_vertex_index) const;
};

#endif //cpp_qmc_VERTEXTYPES_H