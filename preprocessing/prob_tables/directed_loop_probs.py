import numpy as np
from utils import set_probability
from vertex_props import num_vertices
from vertex_props import num_legs_indices
from vertex_props import num_diagonal_vertices
from vertex_props import num_legs_per_vertex
from vertex_props import vertex_map
from vertex_props import leg_spin
from vertex_props import vertical_swap_mapping
from vertex_props import horizontal_swap_mapping
from vertex_props import composed_swaps_mapping

def compute_transition_weights(gamma, Delta, J_ij, h_B, ksi, dist_dep_gamma):

    if h_B < 0.0:
        raise Exception("h_B needs to be greater than or equal to 0")
    else:
        Delta_over_four_J_ij = (Delta/4.0) * J_ij

        Delta_positive = ((Delta + 1.0)/(2.0)) * J_ij
        Delta_negative = ((Delta - 1.0)/(2.0)) * J_ij
        
        if Delta_over_four_J_ij > h_B:
            offset_b = Delta_over_four_J_ij

            if h_B >= Delta_negative:
                b_3_p = 0.0
                epsilon = -Delta_negative/2.0 + h_B/2.0 + gamma
            else:
                b_3_p = Delta_negative - h_B + ksi
                epsilon = gamma

            if h_B <= -Delta_negative:
                b_3 = 0.0
            else:
                b_3 = Delta_negative + h_B + ksi

            a_p = -Delta_negative/2.0 + h_B/2.0 + b_3_p/2.0
            b_p = Delta_positive/2.0 - h_B/2.0 - b_3_p/2.0
            c_p = Delta_negative/2.0 + epsilon - h_B/2.0 - b_3_p/2.0

            a = -Delta_negative/2.0 - h_B/2.0 + b_3/2.0
            b = Delta_positive/2.0 + h_B/2.0 - b_3/2.0
            c = epsilon + Delta_negative/2.0 + h_B/2.0 -b_3/2.0

            b_1_p = 0.0
            b_2_p = 0.0
            b_1 = 0.0
            b_2 = 0.0

        elif Delta < 0.0:

            if h_B == 0.0:
                offset_b = -(Delta/4.0)*J_ij
                if Delta <= -1.0:
                    if dist_dep_gamma:
                        epsilon = gamma - (Delta/10.0)*J_ij
                        c_p = gamma - (Delta/10.0)*J_ij
                        #epsilon = gamma + 0.1/(r_b_pow_alpha)
                        #c_p = gamma + 0.1/(r_b_pow_alpha)
                    else:
                        epsilon = gamma
                        c_p = gamma
                    c = c_p
                    a_p = (1.0/2.0)*J_ij
                    b_p = 0.0
                    a = a_p
                    b = b_p
                    b_2_p = -((1.0 + Delta)/2.0)*J_ij
                    b_2 = b_2_p
                    b_3_p = 0.0
                    b_1_p = 0.0
                    b_1 = b_1_p
                    b_3 = b_3_p
                else:
                    b_2_p = 0.0
                    b_p = ((1.0 + Delta)/(4.0))*J_ij
                    a_p = ((1.0 - Delta)/(4.0))*J_ij
                    c_p = gamma
                    epsilon = ((1.0 + Delta)/(4.0)) * J_ij + gamma
                    c = c_p
                    a = a_p
                    b = b_p
                    b_1_p = 0.0
                    b_3_p = 0.0
                    b_2 = b_2_p
                    b_1 = b_1_p
                    b_3 = b_3_p
            else:
                offset_b = h_B - Delta_over_four_J_ij

                if h_B <= Delta_positive:
                    b_2_p = 0.0
                    epsilon = Delta_positive/2.0 - h_B/2.0 + gamma
                else:
                    b_2_p = h_B - Delta_positive + ksi
                    epsilon = gamma
                
                if h_B <= -Delta_positive:
                    b_2 = -h_B - Delta_positive + ksi
                else:
                    b_2 = 0.0

                if h_B <= -Delta_negative:
                    b_3 = 0.0
                else: 
                    b_3 = h_B + Delta_negative + ksi
                
                a_p = -Delta_negative/2.0 + h_B/2.0 - b_2_p/2.0
                b_p = Delta_positive/2.0 - h_B/2.0 + b_2_p/2.0
                c_p = epsilon - Delta_positive/2.0 + h_B/2.0 - b_2_p/2.0

                a = -Delta_negative/2.0 - h_B/2.0 + b_3/2.0 - b_2/2.0
                b = Delta_positive/2.0 + h_B/2.0 + b_2/2.0 - b_3/2.0
                c = 3.0*h_B/2.0 + epsilon - Delta_positive/2.0 - b_2/2.0 - b_3/2.0

                b_1 = 0.0
                b_1_p = 0.0
                b_3_p = 0.0
        
        else:
            offset_b = h_B
            one_over_four_J_ij = (1.0/4.0) * J_ij

            if h_B <= Delta_negative:
                b_3_p = Delta_negative - h_B + ksi
                epsilon = gamma
            else:
                b_3_p = 0.0
                if h_B <= Delta_positive and h_B <= 2.0*one_over_four_J_ij:
                    epsilon = one_over_four_J_ij - h_B/2.0 + gamma
                else:
                    epsilon = gamma

            if h_B <= Delta_positive:
                b_2_p = 0.0
            else:
                b_2_p = h_B - Delta_positive + ksi

            if h_B < -Delta_negative:
                b_3 = 0
            else:
                b_3 = h_B + Delta_negative + ksi
            
            a_p = -Delta_negative/2.0 + h_B/2.0 + b_3_p/2.0 - b_2_p/2.0
            b_p = Delta_positive/2.0 - h_B/2.0 - b_3_p/2.0 + b_2_p/2.0
            c_p = epsilon - one_over_four_J_ij + h_B/2.0 - b_3_p/2.0 - b_2_p/2.0

            a = -Delta_negative/2.0 - h_B/2.0 + b_3/2.0
            b = Delta_positive/2.0 + h_B/2.0 - b_3/2.0
            c = epsilon - one_over_four_J_ij + 3.0*h_B/2.0 - b_3/2.0

            b_1 = 0.0
            b_2 = 0.0
            b_1_p = 0.0

        offset_b += epsilon
        return (a,b,c,a_p,b_p,c_p,b_1,b_2,b_3,b_1_p,b_2_p,b_3_p,epsilon,offset_b)
    
def update_directed_loop_probs(vertex_enum, l_e, l_x, bond, transition_weight, vertex_weights, prob_table):

    init_composite_leg_index, init_row_index = get_composite_row_prob_index(vertex_enum, l_e, l_x)
    normalization = vertex_weights[vertex_enum, bond]
    if normalization != 0.0:
        prob_table[init_row_index, bond] = set_probability(transition_weight/normalization)

    new_vertex_enum, new_l_e, new_l_x = get_symmetric_indices(vertex_enum, l_e, l_x, vertical_swap_mapping)
    new_composite_leg_index, new_row_index = get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
    prob_table[new_row_index, bond] = prob_table[init_row_index, bond]

    new_vertex_enum, new_l_e, new_l_x = get_symmetric_indices(vertex_enum, l_e, l_x, horizontal_swap_mapping)
    new_composite_leg_index, new_row_index = get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
    prob_table[new_row_index, bond] = prob_table[init_row_index, bond]

    new_vertex_enum, new_l_e, new_l_x = get_symmetric_indices(vertex_enum, l_e, l_x, composed_swaps_mapping)
    new_composite_leg_index, new_row_index = get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
    prob_table[new_row_index, bond] = prob_table[init_row_index, bond]

    return prob_table

def get_composite_row_prob_index(vertex_enum, entrance_leg_enum, exit_leg_enum):

    composite_leg_index = num_legs_per_vertex*entrance_leg_enum + exit_leg_enum
    row_index = num_legs_indices*vertex_enum + composite_leg_index

    return composite_leg_index, row_index

def get_symmetric_indices(vertex_enum, entrance_leg_enum, exit_leg_enum, symmetry_leg_mapping):

    init_spin_tuple = leg_spin[vertex_enum]
    new_spin_tuple = (init_spin_tuple[symmetry_leg_mapping[0]], init_spin_tuple[symmetry_leg_mapping[1]], 
                      init_spin_tuple[symmetry_leg_mapping[2]], init_spin_tuple[symmetry_leg_mapping[3]])
    new_vertex_enum = vertex_map[new_spin_tuple]

    new_entrance_leg_enum = symmetry_leg_mapping[entrance_leg_enum]
    new_exit_leg_enum = symmetry_leg_mapping[exit_leg_enum]

    return new_vertex_enum, new_entrance_leg_enum, new_exit_leg_enum

def compute_prob_tables_directed_loops(num_bonds, sites, J_ij_vector, gamma, h_B, Delta, ksi, dist_dep_gamma):

    num_rows = num_vertices * num_legs_indices
    directed_loop_prob_table = np.zeros((num_rows, num_bonds))
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

        diag_prob_table[:, bond] = vertex_weights[0:4, bond]

        weight_vars = compute_transition_weights(gamma, Delta, J_ij, h_B, ksi, dist_dep_gamma)

        '''
        a = np.round(weight_vars[0], 15)
        b = np.round(weight_vars[1], 15)
        c = np.round(weight_vars[2], 15)
        a_p = np.round(weight_vars[3], 15)
        b_p = np.round(weight_vars[4], 15)
        c_p = np.round(weight_vars[5], 15)
        b_1 = np.round(weight_vars[6], 15)
        b_2 = np.round(weight_vars[7], 15)
        b_3 = np.round(weight_vars[8], 15)
        b_1_p = np.round(weight_vars[9], 15)
        b_2_p = np.round(weight_vars[10], 15)
        b_3_p = np.round(weight_vars[11], 15)
        epsilon = np.round(weight_vars[12], 15)
        offset = np.round(weight_vars[13], 15)
        '''

        '''
        a = set_probability(weight_vars[0])
        b = set_probability(weight_vars[1])
        c = set_probability(weight_vars[2])
        a_p = set_probability(weight_vars[3])
        b_p = set_probability(weight_vars[4])
        c_p = set_probability(weight_vars[5])
        b_1 = set_probability(weight_vars[6])
        b_2 = set_probability(weight_vars[7])
        b_3 = set_probability(weight_vars[8])
        b_1_p = set_probability(weight_vars[9])
        b_2_p = set_probability(weight_vars[10])
        b_3_p = set_probability(weight_vars[11])
        epsilon = set_probability(weight_vars[12])
        offset = set_probability(weight_vars[13])
        '''
        
        
        a = weight_vars[0]
        b = weight_vars[1]
        c = weight_vars[2]
        a_p = weight_vars[3]
        b_p = weight_vars[4]
        c_p = weight_vars[5]
        b_1 = weight_vars[6]
        b_2 = weight_vars[7]
        b_3 = weight_vars[8]
        b_1_p = weight_vars[9]
        b_2_p = weight_vars[10]
        b_3_p = weight_vars[11]
        epsilon = weight_vars[12]
        offset = weight_vars[13]

        vertex_weights[0:4,bond] += offset
        spectrum_offset += offset

        diag_prob_table[:,bond] += offset

        diag_prob_table[0,bond] = set_probability(diag_prob_table[0,bond])
        diag_prob_table[1,bond] = set_probability(diag_prob_table[1,bond])
        diag_prob_table[2,bond] = set_probability(diag_prob_table[2,bond])
        diag_prob_table[3,bond] = set_probability(diag_prob_table[3,bond])

        vertex_weights[0,bond] = set_probability(vertex_weights[0,bond])
        vertex_weights[1,bond] = set_probability(vertex_weights[1,bond])
        vertex_weights[2,bond] = set_probability(vertex_weights[2,bond])
        vertex_weights[3,bond] = set_probability(vertex_weights[3,bond])
        vertex_weights[4,bond] = set_probability(vertex_weights[4,bond])
        vertex_weights[5,bond] = set_probability(vertex_weights[5,bond])

        max_over_states[bond] = np.max(diag_prob_table[:,bond])
        diag_prob_table[:,bond] /= max_over_states[bond]
        max_diag_norm += max_over_states[bond]

        vertex_enum = 0
        l_e = 0
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, b_3, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, 0.0, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, c, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, b, vertex_weights, directed_loop_prob_table)

        vertex_enum = 1
        l_e = 0
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, b_2_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, a_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, c_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, 0.0, vertex_weights, directed_loop_prob_table)

        l_e = 1
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, a, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, b_2, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, 0.0, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, c, vertex_weights, directed_loop_prob_table)

        vertex_enum = 5
        l_e = 0
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, b_1, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, a, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, 0.0, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, b, vertex_weights, directed_loop_prob_table)

        l_e = 1
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, a_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, b_1_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, b_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, 0.0, vertex_weights, directed_loop_prob_table)

        vertex_enum = 3
        l_e = 0
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 0, bond, b_3_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 1, bond, 0.0, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 2, bond, c_p, vertex_weights, directed_loop_prob_table)
        directed_loop_prob_table = update_directed_loop_probs(vertex_enum, l_e, 3, bond, b_p, vertex_weights, directed_loop_prob_table)

    max_over_states[:] /= max_diag_norm

    return diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, directed_loop_prob_table