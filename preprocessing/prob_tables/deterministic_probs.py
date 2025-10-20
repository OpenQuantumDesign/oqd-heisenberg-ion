import numpy as np

def compute_probability_tables_deterministic(num_bonds, J_ij_vector, hamiltonian_type):

    max_over_states = np.zeros(num_bonds)
    max_norm = 0.0

    for bond in range(num_bonds):
        
        J_ij = J_ij_vector[bond]
        max_over_states[bond] = 0.5 * J_ij

        max_norm += 0.5 * J_ij
    
    if hamiltonian_type == 0:
        spectrum_offset = max_norm
        Delta = 0.0
    elif hamiltonian_type == 1:
        spectrum_offset = 0.5*max_norm
        Delta = 1.0
    elif hamiltonian_type == -1:
        spectrum_offset = 0.5*max_norm
        Delta = -1.0
    else:
        raise Exception("Invalid hamiltonian type provided\n")
    
    max_over_states[:] /= max_norm

    return max_over_states, max_norm, spectrum_offset, Delta