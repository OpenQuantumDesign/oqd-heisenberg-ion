import numpy as np

def get_J_ij_power_law(num_bonds, distances, alpha):

    J_ij_vector = np.zeros(num_bonds)
    for b in range(num_bonds):
        J_ij_vector[b] = 1.0/(distances[b])**alpha

    return J_ij_vector

def get_J_ij_from_matrix(N, num_bonds, J_ij_file):

    J_ij_matrix = np.loadtxt(J_ij_file, delimiter=',', skiprows=1)
    J_ij_vector = np.zeros(num_bonds)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            J_ij_vector[b] = J_ij_matrix[i,j]
            b += 1

    return J_ij_vector

def get_J_ij_vector(interaction_type, geometry, **interaction_args):

    if interaction_type == "Power-Law":
        J_ij_vector = get_J_ij_power_law(geometry.num_bonds, geometry.distances, **interaction_args)
    elif interaction_type == "Input-Matrix":
        J_ij_vector = get_J_ij_from_matrix(geometry.N, geometry.num_bonds, **interaction_args)
    else:
        raise ValueError("Interaction type: {} not recognized. Available types are 'Power-Law' and" \
        " 'Input-Matrix'".format(interaction_type))
    
    return J_ij_vector