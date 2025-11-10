import numpy as np

num_vertices_isotropic = 4
num_vertices_general = 6

num_legs_per_vertex = 4

num_diagonal_vertices_isotropic = 2
num_diagonal_vertices_general = 4

num_legs_indices = num_legs_per_vertex**2

anti_ferro_vertex_map = {(1,-1,1,-1):0,
              (-1,1,-1,1):1,
              (1,-1,-1,1):2,
              (-1,1,1,-1):3}

ferro_vertex_map = {(1,1,1,1):0,
              (-1,-1,-1,-1):1,
              (1,-1,-1,1):2,
              (-1,1,1,-1):3}

general_vertex_map = {(1,1,1,1):0,
              (1,-1,1,-1):1,
              (-1,1,-1,1):2,
              (-1,-1,-1,-1):3,
              (1,-1,-1,1):4,
              (-1,1,1,-1):5}

general_operator_type = (0,0,0,0,1,1)
isotropic_operator_type = (0,0,1,1)

anti_ferro_leg_spin = [(1,-1,1,-1), (-1,1,-1,1), (1,-1,-1,1), (-1,1,1,-1)]
ferro_leg_spin = [(1,1,1,1), (-1,-1,-1,-1), (1,-1,-1,1), (-1,1,1,-1)]
general_leg_spin = [(1,1,1,1), (1,-1,1,-1), (-1,1,-1,1), (-1,-1,-1,-1), (1,-1,-1,1), (-1,1,1,-1)]
# <1|S_z|1> = 1/2, <-1|S_z|-1> = -1/2

def get_new_vertex(l_e, l_x, vertex_type, vertex_mapping, leg_spin):

    if l_e == l_x:
        return vertex_type
    else:
        spin_configs = list(leg_spin[vertex_type])
        spin_configs[l_e] = -spin_configs[l_e]
        spin_configs[l_x] = -spin_configs[l_x]

        spin_configs_tuple = tuple(spin_configs)

        if spin_configs_tuple in leg_spin:
            return vertex_mapping[tuple(spin_configs)]
        else:
            return -1

def generate_new_vertex_type_array(vertex_mapping, num_vertices, leg_spin):

    num_cols = num_legs_per_vertex**2
    new_vertex_types = np.zeros((num_vertices,num_cols), dtype=int)
    for i in range(num_vertices):
        for j in range(num_cols):
            l_e = j // num_legs_per_vertex
            l_x = j % num_legs_per_vertex
            new_vertex_types[i,j] = get_new_vertex(l_e, l_x, i, vertex_mapping, leg_spin)

    return new_vertex_types

def generate_XY_exit_legs(out_vertex_type_map):

    xy_exit_leg_combinations = np.zeros((num_vertices_general*num_legs_per_vertex, 2), dtype=int)
    for i in range(num_vertices_general):
        for j in range(num_legs_per_vertex):
            composite_index = num_legs_per_vertex * i + j
            count = 0
            for k in range(num_legs_per_vertex):
                col_index = num_legs_per_vertex * j + k
                if out_vertex_type_map[i, col_index] != -1 and j != k:
                    xy_exit_leg_combinations[composite_index, count] = k
                    count += 1
    
    return xy_exit_leg_combinations

anti_ferro_isotropic_exit_leg_map = (1, 0, 3, 2)
ferro_isotropic_exit_leg_map = (3, 2, 1, 0)

isotropic_exit_leg_maps = [ferro_isotropic_exit_leg_map, anti_ferro_isotropic_exit_leg_map]
isotropic_vertex_map = [ferro_vertex_map, anti_ferro_vertex_map]
isotropic_leg_spin = [ferro_leg_spin, anti_ferro_leg_spin]

isotropic_interaction_type_map = {"FMHeisenbergFMZ": 0, "AFMHeisenbergFMZ": 1}
