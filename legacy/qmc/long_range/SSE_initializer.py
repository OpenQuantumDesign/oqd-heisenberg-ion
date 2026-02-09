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
num_legs_indices = num_legs_per_vertex**2

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

def compute_transition_weights(gamma, Delta, r_b_pow_alpha, h_B, ksi):

    if h_B < 0.0:
        raise Exception("h_B needs to be greater than or equal to 0")
    else:
        Delta_over_four_rb_pow_alpha = Delta/(4.0 * r_b_pow_alpha)

        Delta_positive = (Delta + 1.0)/(2.0 * r_b_pow_alpha)
        Delta_negative = (Delta - 1.0)/(2.0 * r_b_pow_alpha)
        
        if Delta_over_four_rb_pow_alpha > h_B:
            offset_b = Delta_over_four_rb_pow_alpha

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
            offset_b = h_B - Delta_over_four_rb_pow_alpha

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
            one_over_four_rb_pow_alpha = (1.0/(4.0 * r_b_pow_alpha))

            if h_B <= Delta_negative:
                b_3_p = Delta_negative - h_B + ksi
                epsilon = gamma
            else:
                b_3_p = 0.0
                if h_B <= Delta_positive and h_B <= 2.0*one_over_four_rb_pow_alpha:
                    epsilon = one_over_four_rb_pow_alpha - h_B/2.0 + gamma
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
            c_p = epsilon - one_over_four_rb_pow_alpha + h_B/2.0 - b_3_p/2.0 - b_2_p/2.0

            a = -Delta_negative/2.0 - h_B/2.0 + b_3/2.0
            b = Delta_positive/2.0 + h_B/2.0 - b_3/2.0
            c = epsilon - one_over_four_rb_pow_alpha + 3.0*h_B/2.0 - b_3/2.0

            b_1 = 0.0
            b_2 = 0.0
            b_1_p = 0.0

        offset_b += epsilon
        return (a,b,c,a_p,b_p,c_p,b_1,b_2,b_3,b_1_p,b_2_p,b_3_p,epsilon,offset_b)
    
def compute_offset(gamma, Delta, r_b_pow_alpah, h_B):

    if h_B < 0.0: 
        raise Exception("h_B needs to be greater than or equal to 0")
    else:
        Delta_over_four_rb_pow_alpha = Delta/(4.0 * r_b_pow_alpah)
        
        if Delta_over_four_rb_pow_alpha > h_B:
            offset_b = Delta_over_four_rb_pow_alpha

        elif Delta < 0.0:
            offset_b = h_B - Delta_over_four_rb_pow_alpha
        
        else:
            offset_b = h_B

        offset_b += gamma
        return offset_b

# Generating this table might be slow for large systems because size grows as N^2. Simpler but slightly slower approach would be to
# compute 3 non-zero heat bath probabilities on the fly 
def compute_prob_tables_heat_bath(num_bonds, sites, alpha, gamma, h_B, Delta):

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
        abs_index_diff_pow_alpha = (j_b - i_b)**alpha

        vertex_weights[0, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) + h_B
        vertex_weights[1, bond] = -((Delta/4.0) * (1.0/abs_index_diff_pow_alpha))
        vertex_weights[2, bond] = vertex_weights[1, bond]
        vertex_weights[3, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) - h_B
        vertex_weights[4, bond] = 0.5/(abs_index_diff_pow_alpha)
        vertex_weights[5, bond] = vertex_weights[4, bond]

        diag_prob_table[:, bond] = vertex_weights[0:4, bond]

        offset = compute_offset(gamma, Delta, abs_index_diff_pow_alpha, h_B)

        vertex_weights[0:4,bond] += offset
        spectrum_offset += offset

        diag_prob_table[:,bond] += offset

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
                        heat_bath_prob_table[row_index, bond] = vertex_weights[new_vertex, bond]
                        norm += vertex_weights[new_vertex, bond]

                heat_bath_prob_table[row_index-num_legs_per_vertex+1:row_index+1, bond] /= norm

    max_over_states[:] /= max_diag_norm
               
    return diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, heat_bath_prob_table

def compute_prob_tables_directed_loops(num_bonds, sites, alpha, gamma, h_B, Delta, ksi):

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
        abs_index_diff_pow_alpha = (j_b - i_b)**alpha

        vertex_weights[0, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) + h_B
        vertex_weights[1, bond] = -((Delta/4.0) * (1.0/abs_index_diff_pow_alpha))
        vertex_weights[2, bond] = vertex_weights[1, bond]
        vertex_weights[3, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) - h_B
        vertex_weights[4, bond] = 0.5/(abs_index_diff_pow_alpha)
        vertex_weights[5, bond] = vertex_weights[4, bond]

        diag_prob_table[:, bond] = vertex_weights[0:4, bond]

        weight_vars = compute_transition_weights(gamma, Delta, abs_index_diff_pow_alpha, h_B, ksi)

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

        vertex_weights[0:4,bond] += offset
        spectrum_offset += offset

        diag_prob_table[:,bond] += offset

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

def get_composite_row_prob_index(vertex_enum, entrance_leg_enum, exit_leg_enum):

    composite_leg_index = num_legs_per_vertex*entrance_leg_enum + exit_leg_enum
    row_index = num_legs_indices*vertex_enum + composite_leg_index

    return composite_leg_index, row_index

vertical_swap_mapping = [2,3,0,1]
horizontal_swap_mapping = [1,0,3,2]
composed_swaps_mapping = [3,2,1,0]

def get_symmetric_indices(vertex_enum, entrance_leg_enum, exit_leg_enum, symmetry_leg_mapping):

    init_spin_tuple = leg_spin[vertex_enum]
    new_spin_tuple = (init_spin_tuple[symmetry_leg_mapping[0]], init_spin_tuple[symmetry_leg_mapping[1]], 
                      init_spin_tuple[symmetry_leg_mapping[2]], init_spin_tuple[symmetry_leg_mapping[3]])
    new_vertex_enum = vertex_map[new_spin_tuple]

    new_entrance_leg_enum = symmetry_leg_mapping[entrance_leg_enum]
    new_exit_leg_enum = symmetry_leg_mapping[exit_leg_enum]

    return new_vertex_enum, new_entrance_leg_enum, new_exit_leg_enum

def update_directed_loop_probs(vertex_enum, l_e, l_x, bond, transition_weight, vertex_weights, prob_table):

    init_composite_leg_index, init_row_index = get_composite_row_prob_index(vertex_enum, l_e, l_x)
    normalization = vertex_weights[vertex_enum, bond]
    if normalization != 0.0:
        prob_table[init_row_index, bond] = transition_weight/normalization

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