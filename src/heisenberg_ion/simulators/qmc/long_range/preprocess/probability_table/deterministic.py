import numpy as np
import os
from .base import ProbabilityTable

class Deterministic(ProbabilityTable):

    args = {}
    allowed_hamiltonians = {"XY", "fm_heisenberg_afm_Z", "fm_heisenberg_fm_Z"}

    def __init__(self, system):

        super().__init__(system)

        self.validate_system()

        self.build()

    
    def validate_system(self):

        super().validate_system()

        hamiltonian_name = self.system.hamiltonian_parameters.hamiltonian_name

        if hamiltonian_name not in self.allowed_hamiltonians:
            raise Exception("Inconsistent hamiltonian and sampling types. Deterministic probability tables " \
            "only support the following types: {}".format(self.allowed_hamiltonians))
    
    
    def build(self):

        self.compute_max_over_states(self.system.geometry.num_bonds, self.system.interactions.J_ij_vector)
        self.compute_spectrum_offset(self.system.hamiltonian_parameters.hamiltonian_name)

        return 0
    

    def compute_spectrum_offset(self, hamiltonian_name):

        if hamiltonian_name == "XY":
            self.spectrum_offset = self.max_diag_norm
        elif hamiltonian_name == "fm_heisenberg_afm_Z" or hamiltonian_name == "fm_heisenberg_fm_Z":
            self.spectrum_offset = 0.5*self.max_diag_norm
        else:
            raise ValueError("Invalid Hamiltonian type: {} provided for deterministic probability tables. " \
            "Allowed types are {}".format(hamiltonian_name, self.allowed_hamiltonians))
        
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

        geometry_file_name = os.path.join(self.prob_dir, "geometry.csv")
        max_over_states_file_name = os.path.join(self.prob_dir, "max_over_states.csv")

        geometry_table = self.system.geometry.geometry_table
        num_bonds = self.system.geometry.num_bonds

        np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", 
                    header="NumBonds={}".format(num_bonds))
        
        header="norm={},spectrum_offset={},loop_update_type={}".format(self.max_diag_norm, 
            self.spectrum_offset, "deterministic")
        
        np.savetxt(max_over_states_file_name, self.max_over_states, delimiter=",", header=header)

        return 0