#include "VertexTypes.h"

// Can implement more constructors here for special cases as necessary
VertexTypes::VertexTypes(){

    num_vertices = 6;
    num_legs_per_vertex = 4;
    num_composte_leg_indices = 16;

    std::vector<int> config_1 = {1, 1, 1, 1};
    std::vector<int> config_2 = {1, -1, 1, -1};
    std::vector<int> config_3 = {-1, 1, -1, 1};
    std::vector<int> config_4 = {-1, -1, -1, -1};
    std::vector<int> config_5 = {1, -1, -1, 1};
    std::vector<int> config_6 = {-1, 1, 1, -1};

    config_to_index_mapping[config_1] = 0;
    config_to_index_mapping[config_2] = 1;
    config_to_index_mapping[config_3] = 2;
    config_to_index_mapping[config_4] = 3;
    config_to_index_mapping[config_5] = 4;
    config_to_index_mapping[config_6] = 5;

    is_off_diag = {0,0,0,0,1,1};

    index_to_config_mapping = {config_1, config_2, config_3, config_4, config_5, config_6};

    populateAllowedExitLegs();
}

// Can also construct a map that does this via indices like in the python version
int VertexTypes::getFlippedSpinsVertexIndex(const int &l_e, const int &l_x, const int &old_vertex_index) const
{

    std::vector<int> old_config = getVertexConfig(old_vertex_index);

    std::vector<int> new_config = flipInputOutputLegs(l_e, l_x, old_config);

    int new_index = getVertexTypeIndex(new_config);

    return new_index;
};

std::vector<int> VertexTypes::flipInputOutputLegs(const int &l_e, const int &l_x, const std::vector<int> &old_config) const {


    std::vector<int> new_config = {old_config.at(0), old_config.at(1),
                                   old_config.at(2), old_config.at(3)};

    new_config.at(l_e) = -new_config.at(l_e);
    new_config.at(l_x) = -new_config.at(l_x);

    return new_config;

}

void VertexTypes::populateAllowedExitLegs() {

    for (int i = 0; i < num_vertices; i++) {
        std::vector<int> old_config = getVertexConfig(i);
        for (int j = 0; j < num_legs_per_vertex; j++) {
            for (int k = 0; k < num_legs_per_vertex; k++) {
                std::vector<int> new_config = flipInputOutputLegs(j, k, old_config);
                for (int l = 0; l < num_vertices; l++) {
                    if (new_config == index_to_config_mapping.at(l)) {
                        allowed_exit_legs.push_back(k);
                        break;
                    }
                }
            }
        }
    }
}

// Need to add the following code to this
// map from the vertex type index and input leg to allowed exit legs
// std::vector<std::vector<std::vector<int>>> allowed_exit_legs;
// for each vertex type
// for each input leg
// for each exit leg
// flip input and exit legs
// if exit allowed (can be checked using config_to_index_mapping)
// add composite index = vertex type index * 16 + input leg * 4 + exit leg to allowed_exit_legs
// then get the output leg indices for prob_table indexing in config_generator from this vector