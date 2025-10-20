import numpy as np

def get_J_ij_power_law(num_bonds, distances, alpha):

    J_ij_vector = np.zeros(num_bonds)
    for b in range(num_bonds):
        J_ij_vector[b] = 1.0/(distances[b])**alpha

    return J_ij_vector

def get_J_ij_exp(N, num_bonds, filepath):

    J_ij_matrix = np.loadtxt(filepath, delimiter=',', skiprows=1)
    J_ij_vector = np.zeros(num_bonds)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            J_ij_vector[b] = J_ij_matrix[i,j]
            b += 1

    return J_ij_vector