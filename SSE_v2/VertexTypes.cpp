#include "VertexTypes.h"

// Can implement more constructors here for special cases as necessary
VertexTypes::VertexTypes(const int &hamiltonian_type){

    num_legs_per_vertex = 4;
    num_composite_leg_indices = 16;

    if (hamiltonian_type == 0) {
        setVertexMappings();
        populateAllowedExitLegsXY();
    }
    else if (hamiltonian_type == 1) {
        setVertexMappingsIsotropicAFM();
        allowed_exit_legs = {3, 2, 1, 0};
    }
    else if (hamiltonian_type == -1) {
        setVertexMappingsIsotropicFM();
        allowed_exit_legs = {1, 0, 3, 2};
    }
    else if (hamiltonian_type == 2 || hamiltonian_type == 3) {
        setVertexMappings();
        populateAllowedExitLegs();
    }
}

int VertexTypes::getVertexTypeIndex(const std::vector<int> &config) const {
    int config_key = 8*config[0] + 4*config[1] + 2*config[2] + config[3];
    int vertex_type = config_to_index_mapping.at(config_key);
    return vertex_type;
}

void VertexTypes::setVertexMappings() {

    num_vertices = 6;

    std::vector<int> config_1 = {1, 1, 1, 1};
    std::vector<int> config_2 = {1, -1, 1, -1};
    std::vector<int> config_3 = {-1, 1, -1, 1};
    std::vector<int> config_4 = {-1, -1, -1, -1};
    std::vector<int> config_5 = {1, -1, -1, 1};
    std::vector<int> config_6 = {-1, 1, 1, -1};

    config_to_index_mapping = { {15,0}, {5,1}, {-5,2},
                                {-15,3}, {3,4}, {-3,5}};

    is_off_diag = {0,0,0,0,1,1};
    twist_mapping = {0,0,0,0,-1,1};

    index_to_config_mapping = {config_1, config_2, config_3, config_4, config_5, config_6};

}

void VertexTypes::setVertexMappingsIsotropicAFM() {

    num_vertices = 4;

    std::vector<int> config_1 = {1, 1, 1, 1};
    std::vector<int> config_4 = {-1, -1, -1, -1};
    std::vector<int> config_5 = {1, -1, -1, 1};
    std::vector<int> config_6 = {-1, 1, 1, -1};

    config_to_index_mapping = { {15,0}, {-15,1}, {3,2}, {-3,3}};

    is_off_diag = {0,0,1,1};
    twist_mapping = {0,0,-1,1};

    index_to_config_mapping = {config_1, config_4, config_5, config_6};

}

void VertexTypes::setVertexMappingsIsotropicFM() {

    num_vertices = 4;

    std::vector<int> config_2 = {1, -1, 1, -1};
    std::vector<int> config_3 = {-1, 1, -1, 1};
    std::vector<int> config_5 = {1, -1, -1, 1};
    std::vector<int> config_6 = {-1, 1, 1, -1};

    config_to_index_mapping = {{5,0}, {-5,1},
                               {3,2}, {-3,3}};

    is_off_diag = {0,0,1,1};
    twist_mapping = {0,0,-1,1};

    index_to_config_mapping = {config_2, config_3, config_5, config_6};

}

// Can also construct a map that does this via indices like in the python version
int VertexTypes::getFlippedSpinsVertexIndex(const int &l_e, const int &l_x, const int &old_vertex_index) const
{

    std::vector<int> old_config = getVertexConfig(old_vertex_index);

    old_config.at(l_e) = -old_config.at(l_e);
    old_config.at(l_x) = -old_config.at(l_x);

    int new_index = getVertexTypeIndex(old_config);

    return new_index;
}

std::vector<int> VertexTypes::flipInputOutputLegs(const int &l_e, const int &l_x,
                                                  const std::vector<int> &old_config) const {


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

void VertexTypes::populateAllowedExitLegsXY() {

    for (int i = 0; i < num_vertices; i++) {
        std::vector<int> old_config = getVertexConfig(i);
        for (int j = 0; j < num_legs_per_vertex; j++) {
            for (int k = 0; k < num_legs_per_vertex; k++) {
                if (k != j) {
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
}