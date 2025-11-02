import numpy as np
from preprocessing.prob_tables.utils.math_utils import *
from preprocessing.prob_tables.utils.vertex_utils import *
from probability_tables import ProbabilityTable

class DirectedLoops(ProbabilityTable):

    def __init__(self, system, gamma, ksi, dist_dep_gamma):

        super.__init__(system, gamma=gamma, ksi=ksi, dist_dep_gamma=dist_dep_gamma)

        self.system = system

        self.gamma = gamma
        self.ksi = ksi
        self.dist_dep_gamma = dist_dep_gamma

        self.h_B = self.system.compute_h_B()

        self.build()


    def build(self):

        num_bonds = self.system.geometry.num_bonds

        self.initialize_tables(num_bonds)
        self.set_vertex_enum_transition_weights_map()

        J_ij_vector = self.system.J_ij_vector

        Delta = self.system.hamiltonian_parameters.Delta
        h_B = self.h_B

        gamma = self.gamma
        ksi = self.ksi
        dist_dep_gamma = self.dist_dep_gamma

        self.compute_prob_tables_directed_loops(num_bonds, J_ij_vector, gamma, h_B, Delta, ksi, dist_dep_gamma)

        pass


    def initialize_tables(self, num_bonds):

        num_rows = num_vertices * num_legs_indices

        self.directed_loop_prob_table = np.zeros((num_rows, num_bonds))
        self.diag_prob_table = np.zeros((num_diagonal_vertices, num_bonds))
        self.max_over_states = np.zeros(num_bonds)
        self.vertex_weights = np.zeros((num_vertices, num_bonds))

        self.spectrum_offset = 0.0
        self.max_diag_norm = 0.0

        return 0
    

    def set_vertex_enum_transition_weights_map(self):
        
        self.vertex_weight_label_map = {}

        self.num_vertex_enums = 4
        self.directed_loop_vertex_enums = ["0", "1", "5", "3"]
        self.vertex_enum_weight_list_counts = [1, 2, 2, 1]

        vertex_enum = "0"
        exit_leg_weights_le_0_v_0 = ["b_3", None, "c", "b"] # l_e = 0
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_0]

        vertex_enum = "1"
        exit_leg_weights_le_0_v_1 = ["b_2_p", "a_p", "c_p", 0.0] # l_e = 0
        exit_leg_weights_le_1_v_1 = ["a", "b_2", 0.0, "c"] # l_e = 1
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_1, exit_leg_weights_le_1_v_1]

        vertex_enum = "5"
        exit_leg_weights_le_0_v_5 = ["b_1", "a", 0.0, "b"] # l_e = 0
        exit_leg_weights_le_1_v_5 = ["a_p", "b_1_p", "b_p", 0.0] # l_e = 1
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_5, exit_leg_weights_le_1_v_5]

        vertex_enum = "3"
        exit_leg_weights_le_0_v_3 = ["b_3_p", 0.0, "c_p", "b_p"] # l_e = 0
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_3]

        return 0
    
    
    def update_directed_loop_probs(self, vertex_enum, l_e, l_x, bond, transition_weight):

        init_composite_leg_index, init_row_index = self.get_composite_row_prob_index(vertex_enum, l_e, l_x)
        normalization = self.vertex_weights[vertex_enum, bond]
        if normalization != 0.0:
            self.directed_loop_prob_table[init_row_index, bond] = set_probability(transition_weight/normalization)

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(vertex_enum, l_e, l_x, vertical_swap_mapping)
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(vertex_enum, l_e, l_x, horizontal_swap_mapping)
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(vertex_enum, l_e, l_x, composed_swaps_mapping)
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

        return 0
    

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
    
    
    def update_directed_loop_table(self, bond, transition_weights):

        for i in range(self.num_vertex_enums):

            vertex_enum = self.directed_loop_vertex_enums[i]
            num_unique_entrance_legs = self.vertex_enum_weight_list_counts[i]

            for l_e in range(num_unique_entrance_legs):

                exit_leg_weight_labels = self.vertex_weight_label_map[vertex_enum][l_e]

                for l_x in range(num_legs_per_vertex):

                    l_x_weight_label = exit_leg_weight_labels[l_x]
                    l_x_weight = transition_weights[l_x_weight_label] if l_x_weight_label is not None else 0.0

                    self.update_directed_loop_probs(int(vertex_enum), l_e, l_x, bond, l_x_weight)

        return 0


    def compute_prob_tables_directed_loops(self, num_bonds, J_ij_vector, gamma, h_B, Delta, ksi, dist_dep_gamma):

        self.transition_weights_calculator = LoopTransitionWeights(gamma, Delta, h_B, ksi, dist_dep_gamma)

        for bond in range(num_bonds):
            
            J_ij = J_ij_vector[bond]

            set_vertex_weights(self.vertex_weights, bond, Delta, J_ij, h_B)

            self.diag_prob_table[:, bond] = self.vertex_weights[0:4, bond]

            self.transition_weights_calculator.compute_transition_weights(J_ij)

            transition_weights = self.transition_weights_calculator.transition_weight_container
            offset = self.transition_weights_calculator.offset_b

            self.vertex_weights[0:4,bond] += offset
            self.spectrum_offset += offset

            self.diag_prob_table[:,bond] += offset

            enforce_positive(self.diag_prob_table)
            enforce_positive(self.vertex_weights)

            self.update_directed_loop_table(bond, transition_weights)

            self.max_over_states[bond] = np.max(self.diag_prob_table[:,bond])
            self.diag_prob_table[:,bond] /= self.max_over_states[bond]
            self.max_diag_norm += self.max_over_states[bond]

        self.max_over_states[:] /= self.max_diag_norm

        return 0
    
class LoopTransitionWeights:

    def __init__(self, gamma, Delta, h_B, ksi, dist_dep_gamma):

        keys = ["a", "b", "c", "a_p", "b_p", "c_p", "b_1", "b_2", "b_3", "b_1_p", "b_2_p", "b_3_p"]
        self.transition_weight_container = {key: None for key in keys}

        self.Delta = Delta
        self.h_B = h_B

        self.gamma = gamma
        self.ksi = ksi
        self.dist_dep_gamma = dist_dep_gamma

        if self.h_B < 0.0:
            raise Exception("h_B needs to be greater than or equal to 0")


    def populate_unprimed_transition_weights(self, a, b, c):

        self.transition_weight_container["a"] = a
        self.transition_weight_container["b"] = b
        self.transition_weight_container["c"] = c

        return 0
    
    def populate_primed_transition_weights(self, a_p, b_p, c_p):

        self.transition_weight_container["a_p"] = a_p
        self.transition_weight_container["b_p"] = b_p
        self.transition_weight_container["c_p"] = c_p

        return 0
    
    def populate_unprimed_bounce_weights(self, b_1, b_2, b_3):

        self.transition_weight_container["b_1"] = b_1
        self.transition_weight_container["b_2"] = b_2
        self.transition_weight_container["b_3"] = b_3

        return 0
    
    def populate_primed_bounce_weights(self, b_1_p, b_2_p, b_3_p):

        self.transition_weight_container["b_1_p"] = b_1_p
        self.transition_weight_container["b_2_p"] = b_2_p
        self.transition_weight_container["b_3_p"] = b_3_p

        return 0
    

    def tranisiton_weights_small_field(self, Delta_over_four_J_ij, Delta_positive, Delta_negative):

        self.offset_b = Delta_over_four_J_ij

        if self.h_B >= Delta_negative:
            b_3_p = 0.0
            epsilon = -Delta_negative/2.0 + self.h_B/2.0 + self.gamma
        else:
            b_3_p = Delta_negative - self.h_B + self.ksi
            epsilon = self.gamma

        if self.h_B <= -Delta_negative:
            b_3 = 0.0
        else:
            b_3 = Delta_negative + self.h_B + self.ksi

        a_p = -Delta_negative/2.0 + self.h_B/2.0 + b_3_p/2.0
        b_p = Delta_positive/2.0 - self.h_B/2.0 - b_3_p/2.0
        c_p = Delta_negative/2.0 + epsilon - self.h_B/2.0 - b_3_p/2.0

        a = -Delta_negative/2.0 - self.h_B/2.0 + b_3/2.0
        b = Delta_positive/2.0 + self.h_B/2.0 - b_3/2.0
        c = epsilon + Delta_negative/2.0 + self.h_B/2.0 -b_3/2.0

        b_1_p = 0.0
        b_2_p = 0.0
        b_1 = 0.0
        b_2 = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

        return 0
    
    def transition_weights_negative_Delta(self, Delta_over_four_J_ij, Delta_positive, Delta_negative, J_ij):

        if self.h_B == 0.0:
            self.offset_b = -(self.Delta/4.0)*J_ij
            if self.Delta <= -1.0:
                if self.dist_dep_gamma:
                    epsilon = self.gamma - (self.Delta/10.0)*J_ij
                    c_p = self.gamma - (self.Delta/10.0)*J_ij
                    #epsilon = gamma + 0.1/(r_b_pow_alpha)
                    #c_p = gamma + 0.1/(r_b_pow_alpha)
                else:
                    epsilon = self.gamma
                    c_p = self.gamma
                c = c_p
                a_p = (1.0/2.0)*J_ij
                b_p = 0.0
                a = a_p
                b = b_p
                b_2_p = -((1.0 + self.Delta)/2.0)*J_ij
                b_2 = b_2_p
                b_3_p = 0.0
                b_1_p = 0.0
                b_1 = b_1_p
                b_3 = b_3_p
            else:
                b_2_p = 0.0
                b_p = ((1.0 + self.Delta)/(4.0))*J_ij
                a_p = ((1.0 - self.Delta)/(4.0))*J_ij
                c_p = self.gamma
                epsilon = ((1.0 + self.Delta)/(4.0)) * J_ij + self.gamma
                c = c_p
                a = a_p
                b = b_p
                b_1_p = 0.0
                b_3_p = 0.0
                b_2 = b_2_p
                b_1 = b_1_p
                b_3 = b_3_p
        else:
            self.offset_b = self.h_B - Delta_over_four_J_ij

            if self.h_B <= Delta_positive:
                b_2_p = 0.0
                epsilon = Delta_positive/2.0 - self.h_B/2.0 + self.gamma
            else:
                b_2_p = self.h_B - Delta_positive + self.ksi
                epsilon = self.gamma
            
            if self.h_B <= -Delta_positive:
                b_2 = -self.h_B - Delta_positive + self.ksi
            else:
                b_2 = 0.0

            if self.h_B <= -Delta_negative:
                b_3 = 0.0
            else: 
                b_3 = self.h_B + Delta_negative + self.ksi
            
            a_p = -Delta_negative/2.0 + self.h_B/2.0 - b_2_p/2.0
            b_p = Delta_positive/2.0 - self.h_B/2.0 + b_2_p/2.0
            c_p = epsilon - Delta_positive/2.0 + self.h_B/2.0 - b_2_p/2.0

            a = -Delta_negative/2.0 - self.h_B/2.0 + b_3/2.0 - b_2/2.0
            b = Delta_positive/2.0 + self.h_B/2.0 + b_2/2.0 - b_3/2.0
            c = 3.0*self.h_B/2.0 + epsilon - Delta_positive/2.0 - b_2/2.0 - b_3/2.0

            b_1 = 0.0
            b_1_p = 0.0
            b_3_p = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

        return 0
    
    def transition_weights_large_field(self, Delta_positive, Delta_negative, J_ij):

        self.offset_b = self.h_B
        one_over_four_J_ij = (1.0/4.0) * J_ij

        if self.h_B <= Delta_negative:
            b_3_p = Delta_negative - self.h_B + self.ksi
            epsilon = self.gamma
        else:
            b_3_p = 0.0
            if self.h_B <= Delta_positive and self.h_B <= 2.0*one_over_four_J_ij:
                epsilon = one_over_four_J_ij - self.h_B/2.0 + self.gamma
            else:
                epsilon = self.gamma

        if self.h_B <= Delta_positive:
            b_2_p = 0.0
        else:
            b_2_p = self.h_B - Delta_positive + self.ksi

        if self.h_B < -Delta_negative:
            b_3 = 0
        else:
            b_3 = self.h_B + Delta_negative + self.ksi
        
        a_p = -Delta_negative/2.0 + self.h_B/2.0 + b_3_p/2.0 - b_2_p/2.0
        b_p = Delta_positive/2.0 - self.h_B/2.0 - b_3_p/2.0 + b_2_p/2.0
        c_p = epsilon - one_over_four_J_ij + self.h_B/2.0 - b_3_p/2.0 - b_2_p/2.0

        a = -Delta_negative/2.0 - self.h_B/2.0 + b_3/2.0
        b = Delta_positive/2.0 + self.h_B/2.0 - b_3/2.0
        c = epsilon - one_over_four_J_ij + 3.0*self.h_B/2.0 - b_3/2.0

        b_1 = 0.0
        b_2 = 0.0
        b_1_p = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

        return 0

    def compute_transition_weights(self, J_ij):
        
        Delta_over_four_J_ij = (self.Delta/4.0) * J_ij

        Delta_positive = ((self.Delta + 1.0)/(2.0)) * J_ij
        Delta_negative = ((self.Delta - 1.0)/(2.0)) * J_ij
        
        if Delta_over_four_J_ij > self.h_B:
            self.tranisiton_weights_small_field(Delta_over_four_J_ij, Delta_positive, Delta_negative)
        elif self.Delta < 0.0:
            self.transition_weights_negative_Delta(Delta_over_four_J_ij, Delta_positive, Delta_negative, J_ij)
        else:
            self.transition_weights_large_field(Delta_positive, Delta_negative, J_ij)
            
        return 0
    
ProbabilityTable.register(DirectedLoops, "DirectedLoops")