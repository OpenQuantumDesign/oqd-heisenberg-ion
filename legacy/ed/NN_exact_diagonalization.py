import numpy as np
import math
from matplotlib import pyplot as plt

def construct_Hamiltonian_2(Delta, interaction_type_index, J=1.0):

    sigma_z = np.diag([1,-1])
    H_z = 0.25*np.kron(sigma_z, sigma_z)

    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0

    H_int = 0.5*(np.kron(sigma_plus, sigma_minus) + np.kron(sigma_minus, sigma_plus))

    return -((-1)**interaction_type_index) * (Delta*J*H_z + J*H_int)

def construct_Hamiltonian_3(Delta, boundary_condition, interaction_type_index, J=1.0):

    sigma_z = np.diag([1,-1])
    identity_matrix = np.eye(2)
    J_z = J * 0.25 * Delta
    
    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0
    
    H_int = 0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) + np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) + np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix))

    if boundary_condition == "PBC":
        H_z = (np.kron(np.kron(sigma_z, identity_matrix), sigma_z)) + np.kron(np.kron(sigma_z, sigma_z), identity_matrix) + np.kron(identity_matrix, np.kron(sigma_z, sigma_z))
        H_int = 0.5*(np.kron(np.kron(sigma_plus, identity_matrix), sigma_minus) + np.kron(np.kron(sigma_minus, identity_matrix), sigma_plus)) + 0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) + np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) + np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix))
    elif boundary_condition == "OBC":
        H_z = np.kron(np.kron(sigma_z, sigma_z), identity_matrix) + np.kron(identity_matrix, np.kron(sigma_z, sigma_z))
        H_int = 0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) + np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) + np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix))
    else: 
        raise Exception("Boundary conditions must be either OBC or PBC")
    
    return -((-1)**interaction_type_index) * (J_z*H_z + J*H_int)