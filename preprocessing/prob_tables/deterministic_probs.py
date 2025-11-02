import numpy as np
import os
from probability_tables import ProbabilityTable

class Deterministic(ProbabilityTable):

    def __init__(self, system):

        super.__init__(system)

        self.build()
    
    
    def build(self):

        self.compute_max_over_states(self.system.geometry.num_bonds, self.system.interactions.J_ij_vector)
        self.compute_spectrum_offset(self.system.hamiltonian_parameters.hamiltonian_type)

        return 0
    

    def compute_spectrum_offset(self, hamiltonian_type):

        if hamiltonian_type == 0:
            self.spectrum_offset = self.max_diag_norm
        elif hamiltonian_type == 1 or hamiltonian_type == -1:
            self.spectrum_offset = 0.5*self.max_diag_norm
        else:
            raise ValueError("Invalid Hamiltonian type: {} provided for deterministic probability tables. " \
            "Allowed types are -1, 0 and 1".format(hamiltonian_type))
        
        return 0
    

    def compute_max_over_states(self, num_bonds, J_ij_vector):

        max_over_states = np.zeros(num_bonds)
        max_diag_norm = 0.0

        for bond in range(num_bonds):
        
            J_ij = J_ij_vector[bond]
            max_over_states[bond] = 0.5 * J_ij

            max_diag_norm += 0.5 * J_ij

        self.max_over_states = max_over_states/max_diag_norm
        self.max_diag_norm = max_diag_norm

        return 0
    
    def write_to_files(self, out_dir):

        super().write_to_files(out_dir)

        geometry_file_name = os.path.join(out_dir, "geometry.csv")
        max_over_states_file_name = os.path.join(out_dir, "max_over_states.csv")

        geometry_table = self.system.geometry.geometry_table
        N = self.system.geometry.N
        num_bonds = self.system.geometry.num_bonds

        np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", 
                    header="NumBonds={}".format(num_bonds))
        
        header="norm={},spectrum_offset={},loop_update_type={}".format(self.max_diag_norm, 
            self.spectrum_offset, "Deterministic")
        
        np.savetxt(max_over_states_file_name, self.max_over_states, delimiter=",", header=header)

        return 0


ProbabilityTable.register(Deterministic, "Deterministic")