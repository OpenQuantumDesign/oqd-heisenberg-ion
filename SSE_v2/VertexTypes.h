#ifndef SSE_V2_VERTEXTYPES_H
#define SSE_V2_VERTEXTYPES_H
#include<vector>
#include<map>

class VertexTypes {
private:
    std::map<std::vector<int>, int> config_to_index_mapping;
    std::vector<std::vector<int>> index_to_config_mapping;

    void populateAllowedExitLegs();

    std::vector<int> flipInputOutputLegs(const int &l_e, const int &l_x, const std::vector<int> &old_config) const;

public:

    int num_vertices;
    int num_legs_per_vertex;
    int num_composte_leg_indices;

    VertexTypes();

    std::vector<int> allowed_exit_legs;

    std::vector<int> is_off_diag;

    int getVertexTypeIndex(const std::vector<int> &config) const {
        int vertex_type = config_to_index_mapping.at(config);
        return vertex_type;
    };

    std::vector<int> getVertexConfig(const int &vertex_index) const {
        return index_to_config_mapping.at(vertex_index);
    };

    int getFlippedSpinsVertexIndex(const int &l_e, const int &l_x, const int &old_vertex_index) const;
};

#endif //SSE_V2_VERTEXTYPES_H