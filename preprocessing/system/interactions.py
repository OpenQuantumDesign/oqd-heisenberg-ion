import numpy as np

class Interactions:

    def __init__(self, interaction_type, geometry, **interaction_args):

        if interaction_type == "Power-Law":

            self.alpha = interaction_args["alpha"]
            self.get_J_ij_power_law(geometry.num_bonds, geometry.distances, self.alpha)

        elif interaction_type == "Input-Matrix":
            
            self.file_path = interaction_args["J_ij_File"]
            self.get_J_ij_from_matrix(geometry.N, geometry.num_bonds, self.file_path)

        else:
            raise ValueError("Interaction type: {} not recognized. Available types are 'Power-Law' and" \
            " 'Input-Matrix'".format(interaction_type))
        
    def get_J_ij_power_law(self, num_bonds, distances, alpha):

        self.J_ij_vector = np.zeros(num_bonds)
        for b in range(num_bonds):
            self.J_ij_vector[b] = 1.0/(distances[b])**alpha

        return 0

    def get_J_ij_from_matrix(self, N, num_bonds, J_ij_file):

        J_ij_matrix = np.loadtxt(J_ij_file, delimiter=',', skiprows=1)
        self.J_ij_vector = np.zeros(num_bonds)
        b = 0
        for i in range(N):
            for j in range(i+1,N):
                self.J_ij_vector[b] = J_ij_matrix[i,j]
                b += 1

        return 0