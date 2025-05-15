import numpy as np

def geometry(N, num_bonds):

    sites = np.zeros((num_bonds,2), dtype=int)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            sites[b,0] = i
            sites[b,1] = j
            b += 1

    return sites

num_vertices = 6
num_legs_per_vertex = 4
num_diagonal_vertices = 4

vertex_map = {(1,1,1,1):0,
              (1,-1,1,-1):1,
              (-1,1,-1,1):2,
              (-1,-1,-1,-1):3,
              (1,-1,-1,1):4,
              (-1,1,1,-1):5}

operator_type = (0,0,0,0,1,1)

leg_spin = [(1,1,1,1), (1,-1,1,-1), (-1,1,-1,1), (-1,-1,-1,-1), (1,-1,-1,1), (-1,1,1,-1)]
# <1|S_z|1> = 1/2, <-1|S_z|-1> = -1/2

def get_new_vertex(l_e, l_x, vertex_type):

    if l_e == l_x:
        return vertex_type
    else:
        spin_configs = list(leg_spin[vertex_type])
        spin_configs[l_e] = -spin_configs[l_e]
        spin_configs[l_x] = -spin_configs[l_x]

        spin_configs_tuple = tuple(spin_configs)

        if spin_configs_tuple in leg_spin:
            return vertex_map[tuple(spin_configs)]
        else: 
            return -1

def generate_new_vertex_type_array():

    num_cols = num_legs_per_vertex**2
    new_vertex_types = np.zeros((num_vertices,num_cols), dtype=int)
    for i in range(num_vertices):
        for j in range(num_cols):
            l_e = j // num_legs_per_vertex
            l_x = j % num_legs_per_vertex
            new_vertex_types[i,j] = get_new_vertex(l_e, l_x, i)

    return new_vertex_types

new_vertex_map = generate_new_vertex_type_array()

# TODO: Can produce diag update table and vertex weights table together.
def generate_diag_update_table(N, num_bonds, epsilon, Delta, alpha, h_B):

    # Use the two-step method with max norm as discussed here: https://scipost.org/SciPostPhysCore.7.2.016/pdf

    prob_table = np.zeros((num_diagonal_vertices, num_bonds))
    max_over_states = np.zeros(num_bonds)

    b = 0

    norm = 0.0
    for i in range(N):
        for j in range(i+1,N):

            abs_index_diff = j-i
            diff_pow_alpha = abs_index_diff**alpha

            #C_b = (np.abs(Delta)/4.0) * (1.0/diff_pow_alpha) + h_B + epsilon

            prob_table[0, b] = (Delta/4.0) * (1.0/diff_pow_alpha) + h_B
            prob_table[1, b] = - (Delta/4.0) * (1.0/diff_pow_alpha)
            prob_table[2, b] = prob_table[1, b]
            prob_table[3, b] = (Delta/4.0) * (1.0/diff_pow_alpha) - h_B

            C_b = np.max(prob_table[:,b]) + epsilon
            prob_table[:,b] += C_b

            max_over_states[b] = np.max(prob_table[:,b])

            prob_table[:,b] /= max_over_states[b]

            norm += max_over_states[b]
            b += 1

    max_over_states[:] /= norm
    
    return prob_table, max_over_states, norm

def generate_vertex_weights(epsilon, h_B, Delta, alpha, N, num_bonds):

    weights = np.zeros((num_vertices, num_bonds))
    b = 0
    spectrum_offset = 0.0
    for i in range(N):
        for j in range(i+1,N):

            abs_index_diff = j-i
            diff_pow_alpha = abs_index_diff**alpha

            #C_b = ((np.abs(Delta)/4.0) * (1.0/diff_pow_alpha)) + h_B + epsilon

            weights[0, b] = ((Delta/4.0) * (1.0/diff_pow_alpha)) + h_B
            weights[1, b] = -((Delta/4.0) * (1.0/diff_pow_alpha))
            weights[2, b] = weights[1, b]
            weights[3, b] = ((Delta/4.0) * (1.0/diff_pow_alpha)) - h_B
            weights[4, b] = 0.5/(diff_pow_alpha)
            weights[5, b] = weights[4, b]

            C_b = np.max(weights[0:4,b]) + epsilon

            weights[0:4,b] += C_b

            spectrum_offset += C_b

            b+=1

    return weights, spectrum_offset

# Generating this table might be slow for large systems because size grows as N^2. Simpler but slightly slower approach would be to
# compute 3 non-zero heat bath probabilities on the fly 
def generate_heat_bath_prob_table(num_bonds, vertex_weights):

    num_legs_indices = num_legs_per_vertex**2
    num_rows = num_vertices * num_legs_indices
    heat_bath_prob_table = np.zeros((num_rows, num_bonds))

    for bond in range(num_bonds):
        for vertex_enum in range(num_vertices):
            for l_e in range(num_legs_per_vertex):

                norm = 0.0
                count_invalid_vertices = 0

                for l_x in range(num_legs_per_vertex):

                    composite_leg_index = num_legs_per_vertex*l_e + l_x
                    row_index = num_legs_indices*vertex_enum + composite_leg_index

                    #new_vertex = get_new_vertex(v_map, l_spin, l_e, l_x, vertex_enum)
                    new_vertex = new_vertex_map[vertex_enum, composite_leg_index]

                    if new_vertex == -1:
                        count_invalid_vertices += 1
                        heat_bath_prob_table[row_index, bond] = 0.0
                    else:
                        heat_bath_prob_table[row_index, bond] = vertex_weights[new_vertex, bond]
                        norm += vertex_weights[new_vertex, bond]

                heat_bath_prob_table[row_index-num_legs_per_vertex+1:row_index+1, bond] /= norm
               
                
    return heat_bath_prob_table

