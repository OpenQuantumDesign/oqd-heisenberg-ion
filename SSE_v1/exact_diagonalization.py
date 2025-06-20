import numpy as np
import math
from matplotlib import pyplot as plt

# Simple construction for small scale testing
####--------------------------------------------------------------------------------------------------------------------------------------------####
def construct_Hamiltonian_2(Delta, J=1.0, h=1.0):

    sigma_z = np.diag([1,-1])
    H_z = 0.25*np.kron(sigma_z, sigma_z)

    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0

    H_int = 0.5*(np.kron(sigma_plus, sigma_minus) + np.kron(sigma_minus, sigma_plus))

    H_field = 0.5*np.kron(sigma_z, np.eye(2)) + 0.5*np.kron(np.eye(2), sigma_z)

    return -Delta*J*H_z - J*H_int - h*H_field

def construct_Hamiltonian_3(Delta, alpha, J=1.0, h=0.0):

    sigma_z = np.diag([1,-1])
    identity_matrix = np.eye(2)
    J_z = Delta * J * 0.25
    H_z = J_z * (np.kron(np.kron(sigma_z, sigma_z), identity_matrix) + np.kron(identity_matrix, np.kron(sigma_z, sigma_z)) + 
                ((0.5**alpha) * np.kron(np.kron(sigma_z, identity_matrix), sigma_z)))
    
    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0
    
    H_int = J*0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) + 
                   np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + 
                   np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) + 
                   np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix) + 
                   ((0.5**alpha) * (np.kron(np.kron(sigma_plus, identity_matrix), sigma_minus) + 
                                    np.kron(np.kron(sigma_minus, identity_matrix), sigma_plus))))
    
    H_field = h*(0.5*np.kron(sigma_z, np.eye(4)) + 0.5*np.kron(np.eye(4), sigma_z) + 
                 0.5*np.kron(np.kron(identity_matrix, sigma_z), identity_matrix))
    
    return -H_z - H_int -H_field
