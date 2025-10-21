import numpy as np
from preprocessing.prob_tables.utils.math_utils import *
from preprocessing.prob_tables.utils.vertex_utils import *

def compute_offset(gamma, Delta, J_ij, h_B):

    if h_B < 0.0: 
        raise Exception("h_B needs to be greater than or equal to 0")
    else:
        Delta_over_four_J_ij = (Delta/4.0) * J_ij
        
        if Delta_over_four_J_ij > h_B:
            offset_b = Delta_over_four_J_ij

        elif Delta < 0.0:
            offset_b = h_B - Delta_over_four_J_ij
        
        else:
            offset_b = h_B

        offset_b += gamma
        return offset_b

# Generating this table might be slow for large systems because size grows as N^2. Simpler but slightly slower approach would be to
# compute 3 non-zero heat bath probabilities on the fly 
def compute_prob_tables_heat_bath(num_bonds, sites, J_ij_vector, gamma, h_B, Delta):

    num_rows = num_vertices * num_legs_indices
    heat_bath_prob_table = np.zeros((num_rows, num_bonds))

    diag_prob_table = np.zeros((num_diagonal_vertices, num_bonds))
    max_over_states = np.zeros(num_bonds)
    vertex_weights = np.zeros((num_vertices, num_bonds))
    spectrum_offset = 0.0
    max_diag_norm = 0.0

    for bond in range(num_bonds):

        i_b = sites[bond, 0]
        j_b = sites[bond, 1]
        #abs_index_diff_pow_alpha = (j_b - i_b)**alpha
        #abs_index_diff_pow_alpha = (distances[bond])**alpha
        J_ij = J_ij_vector[bond]

        vertex_weights[0, bond] = ((Delta/4.0) * J_ij) + h_B
        vertex_weights[1, bond] = -((Delta/4.0) * J_ij)
        vertex_weights[2, bond] = vertex_weights[1, bond]
        vertex_weights[3, bond] = ((Delta/4.0) * J_ij) - h_B
        vertex_weights[4, bond] = 0.5 * J_ij
        vertex_weights[5, bond] = vertex_weights[4, bond]

        #diag_prob_table[:, bond] = vertex_weights[0:4, bond]
        diag_prob_table[0, bond] = vertex_weights[0, bond]
        diag_prob_table[1, bond] = vertex_weights[1, bond]
        diag_prob_table[2, bond] = vertex_weights[2, bond]
        diag_prob_table[3, bond] = vertex_weights[3, bond]

        offset = compute_offset(gamma, Delta, J_ij, h_B)

        vertex_weights[0:4,bond] += offset
        spectrum_offset += offset

        diag_prob_table[:,bond] += offset

        diag_prob_table[0,bond] = set_probability(diag_prob_table[0,bond])
        diag_prob_table[1,bond] = set_probability(diag_prob_table[1,bond])
        diag_prob_table[2,bond] = set_probability(diag_prob_table[2,bond])
        diag_prob_table[3,bond] = set_probability(diag_prob_table[3,bond])

        max_over_states[bond] = np.max(diag_prob_table[:,bond])
        diag_prob_table[:,bond] /= max_over_states[bond]
        max_diag_norm += max_over_states[bond]

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
                        heat_bath_prob_table[row_index, bond] = set_probability(vertex_weights[new_vertex, bond])
                        norm += vertex_weights[new_vertex, bond]

                heat_bath_prob_table[row_index-num_legs_per_vertex+1:row_index+1, bond] /= norm

    max_over_states[:] /= max_diag_norm
               
    return diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, heat_bath_prob_table