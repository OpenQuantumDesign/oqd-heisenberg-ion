import numpy as np
import os

def geometry(N_1, N_2, geometry_type):
    
    if geometry_type == "1d":
        if N_2 == 0:
            sites, distances, geometry_table, N, num_bonds = geometry_1d(N_1)
        else: 
            raise Exception("Geometry type is 1d but N_2  is non-zero.")
    elif geometry_type == "2d_triangular":
        sites, distances, geometry_table, N, num_bonds = geometry_2d_triangular(N_1, N_2)
    else:
        raise Exception("Geometry type: {} is not implemented. Valid possible choices are: {} and {}.".format(geometry_type, "1d", "2d_triangular"))

    return sites, distances, geometry_table, N, num_bonds

def geometry_1d(N):

    num_bonds = int(N*(N-1)/2)
    sites = np.zeros((num_bonds,2), dtype=int)
    geometry_table = np.zeros((num_bonds,3))
    distances = np.zeros(num_bonds)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            sites[b,0] = i
            sites[b,1] = j
            distances[b] = j-i
            geometry_table[b,0] = i
            geometry_table[b,1] = j
            geometry_table[b,2] = distances[b]
            b += 1

    return sites, distances, geometry_table, N, num_bonds

def geometry_2d_triangular(N_1, N_2):

    N = N_1 * N_2
    num_bonds = int(N*(N-1)/2)
    sites = np.zeros((num_bonds,2), dtype=int)
    distances = np.zeros(num_bonds)
    geometry_table = np.zeros(num_bonds,3)
    b = 0
    a_1 = np.array([1.0,0.0])
    a_2 = np.array([-0.5, np.sqrt(3.0)/2.0])
    for i1 in range(N_1):
        for i2 in range(N_2):
            i = i1*N_2 + i2
            for j1 in range(N_1):
                for j2 in range(N_2):
                    j = j1*N_2 + j2
                    if j > i:
                        sites[b,0] = i
                        sites[b,1] = j
                        distances[b] = np.norm((j1-i1)*a_1 + (j2-i2)*a_2)
                        geometry_table[b,0] = i
                        geometry_table[b,1] = j
                        geometry_table[b,2] = distances[b]
                        b += 1

    return sites, distances, geometry_table, N, num_bonds

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

def compute_transition_weights(gamma, Delta, r_b_pow_alpha, h_B, ksi, dist_dep_gamma):

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

            if h_B == 0.0:
                offset_b = -Delta/(4.0*r_b_pow_alpha)
                if Delta < -1.0:
                    if dist_dep_gamma:
                        epsilon = gamma - Delta/(10.0*r_b_pow_alpha)
                        c_p = gamma - Delta/(10.0*r_b_pow_alpha)
                        #epsilon = gamma + 0.1/(r_b_pow_alpha)
                        #c_p = gamma + 0.1/(r_b_pow_alpha)
                    else:
                        epsilon = gamma
                        c_p = gamma
                    c = c_p
                    a_p = 1.0/(2.0*r_b_pow_alpha)
                    b_p = 0.0
                    a = a_p
                    b = b_p
                    b_2_p = -(1.0 + Delta)/(2*r_b_pow_alpha)
                    b_2 = b_2_p
                    b_3_p = 0.0
                    b_1_p = 0.0
                    b_1 = b_1_p
                    b_3 = b_3_p
                else:
                    b_2_p = 0.0
                    b_p = (1.0 + Delta)/(4.0*r_b_pow_alpha)
                    a_p = (1.0 - Delta)/(4.0*r_b_pow_alpha)
                    c_p = gamma
                    epsilon = ((1.0 + Delta)/(4.0 * r_b_pow_alpha)) + gamma
                    c = c_p
                    a = a_p
                    b = b_p
                    b_1_p = 0.0
                    b_3_p = 0.0
                    b_2 = b_2_p
                    b_1 = b_1_p
                    b_3 = b_3_p
            else:
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
    
def compute_offset(gamma, Delta, r_b_pow_alpha, h_B):

    if h_B < 0.0: 
        raise Exception("h_B needs to be greater than or equal to 0")
    else:
        Delta_over_four_rb_pow_alpha = Delta/(4.0 * r_b_pow_alpha)
        
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
def compute_prob_tables_heat_bath(num_bonds, sites, distances, alpha, gamma, h_B, Delta):

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
        abs_index_diff_pow_alpha = (distances[bond])**alpha

        vertex_weights[0, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) + h_B
        vertex_weights[1, bond] = -((Delta/4.0) * (1.0/abs_index_diff_pow_alpha))
        vertex_weights[2, bond] = vertex_weights[1, bond]
        vertex_weights[3, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) - h_B
        vertex_weights[4, bond] = 0.5/(abs_index_diff_pow_alpha)
        vertex_weights[5, bond] = vertex_weights[4, bond]

        #diag_prob_table[:, bond] = vertex_weights[0:4, bond]
        diag_prob_table[0, bond] = vertex_weights[0, bond]
        diag_prob_table[1, bond] = vertex_weights[1, bond]
        diag_prob_table[2, bond] = vertex_weights[2, bond]
        diag_prob_table[3, bond] = vertex_weights[3, bond]

        offset = compute_offset(gamma, Delta, abs_index_diff_pow_alpha, h_B)

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

def set_probability(val):

    if val < 0.0:
        if np.abs(val) < 1e-15:
            print("Correcting small negative")
            return 0.0
        else:
            print("invalid probability encountered: ", val)
    else:
        return val 

def compute_prob_tables_directed_loops(num_bonds, sites, distances, alpha, gamma, h_B, Delta, ksi, dist_dep_gamma):

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
        abs_index_diff_pow_alpha = (distances[bond])**alpha

        vertex_weights[0, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) + h_B
        vertex_weights[1, bond] = -((Delta/4.0) * (1.0/abs_index_diff_pow_alpha))
        vertex_weights[2, bond] = vertex_weights[1, bond]
        vertex_weights[3, bond] = ((Delta/4.0) * (1.0/abs_index_diff_pow_alpha)) - h_B
        vertex_weights[4, bond] = 0.5/(abs_index_diff_pow_alpha)
        vertex_weights[5, bond] = vertex_weights[4, bond]

        diag_prob_table[:, bond] = vertex_weights[0:4, bond]

        weight_vars = compute_transition_weights(gamma, Delta, abs_index_diff_pow_alpha, h_B, ksi, dist_dep_gamma)

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

def write_prob_tables(gamma, Delta, h, N_1, N_2, lattice_type, alpha, ksi, J, loop_update_type, dist_dep_gamma):

    out_dir = "ProbabilityDensities"
    sites, distances, geometry_table, N, num_bonds = geometry(N_1, N_2, lattice_type)
    num_bonds = int(N*(N-1)/2.0)
    h_B = h/(J*(N-1))

    if loop_update_type == "heat_bath":
        prob_tables = compute_prob_tables_heat_bath(num_bonds, sites, distances, alpha, gamma, h_B, Delta)
        diag_prob_table = prob_tables[0]
        max_over_states = prob_tables[1]
        max_diag_norm = prob_tables[2]
        vertex_weights = prob_tables[3]
        spectrum_offset = prob_tables[4]
        loop_update_prob_table = prob_tables[5]

    elif loop_update_type == "directed_loops": 
        prob_tables = compute_prob_tables_directed_loops(num_bonds, sites, distances, alpha, gamma, h_B, Delta, ksi, dist_dep_gamma)
        diag_prob_table = prob_tables[0]
        max_over_states = prob_tables[1]
        max_diag_norm = prob_tables[2]
        vertex_weights = prob_tables[3]
        spectrum_offset = prob_tables[4]
        loop_update_prob_table = prob_tables[5]

    else:
        raise Exception("Invalid key word argument for loop_update_type")
    
    if dist_dep_gamma:
        dist_dep_offset = 1
    else:
        dist_dep_offset = 0
    file_prefix = "N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_{}_J_{}_dist_dep_offset_{}".format(N, Delta, h, alpha, gamma, ksi, J, dist_dep_offset)

    geometry_file_name = os.path.join(out_dir, "N_{}_geometry.csv".format(N))
    diag_file_name = os.path.join(out_dir, file_prefix + "_diag_probs.csv")
    max_over_states_file_name = os.path.join(out_dir, file_prefix + "_max_over_states.csv")
    loop_update_table_file_name = os.path.join(out_dir, file_prefix + "_{}_off_diag_table.csv".format(loop_update_type))
    vertex_weights_file_name = os.path.join(out_dir, file_prefix + "_vertex_weights.csv".format(loop_update_type))

    np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", header="N={}, NumBonds={}".format(N, num_bonds))
    np.savetxt(diag_file_name, diag_prob_table, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={}".format(N, Delta, h, alpha, gamma, ksi, J))
    np.savetxt(vertex_weights_file_name, vertex_weights, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={}".format(N, Delta, h, alpha, gamma, ksi, J))
    np.savetxt(max_over_states_file_name, max_over_states, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={},norm={}".format(N, Delta, h, alpha, gamma, ksi, J, max_diag_norm))
    np.savetxt(loop_update_table_file_name, loop_update_prob_table, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={},spectrum_offset={},loop_update_type={}".format(N, Delta, h, alpha, gamma, ksi, J, spectrum_offset, loop_update_type))

    return 0

if __name__=="__main__":
    
    #Delta_list = [-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
    #Delta_list = [-3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
    Delta_list = [-20.0,-19.0,-18.0,-17.0,-11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-10.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
    #Delta_list = [-3.0]

    h_list = [0.0]
    N_1 = 9
    N_2 = 0
    alpha = 1.0
    ksi = 0.0
    J = 1.0
    loop_update_type = "directed_loops"
    lattice_type = "1d"
    dist_dep_gamma = True

    for Delta in Delta_list:
        if Delta < 0.0:
            gamma = 0.0
        else:
            gamma = 0.1
        for h in h_list:
            write_prob_tables(gamma, Delta, h, N_1, N_2, lattice_type, alpha, ksi, J, loop_update_type, dist_dep_gamma)