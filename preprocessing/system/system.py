# TODO: Implement a System class that combines interaction vector and geometry based on inputs
from geometry import Geometry
from interactions import *

class System:

    def __init__(self, geometry_type, geometry_args, interaction_type, interaction_args):
        
        self.geometry_type = geometry_type
        self.geometry = Geometry.create(geometry_type, **geometry_args)

        self.interaction_type = interaction_type
        if interaction_type == "Power-Law":
            self.J_ij_vector = get_J_ij_power_law(self.geometry.num_bonds, self.geometry.distances, **interaction_args)
        elif interaction_type == "Input-Matrix":
            self.J_ij_vector = get_J_ij_from_matrix(self.geometry.N, self.geometry.num_bonds, **interaction_args)
        else:
            raise ValueError("Interaction type: {} not recognized. Available types are 'Power-Law' and 'Input-Matrix'".format(interaction_type))

        return 0

