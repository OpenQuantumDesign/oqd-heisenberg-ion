from geometry import Geometry
from hamiltonian import HamiltonianParameters
from interactions import Interactions 

class System:

    def __init__(self, hamiltonian_args, geometry_args, interaction_args):

        self.model_name = hamiltonian_args["Type"]
        self.hamiltonian_parameters = HamiltonianParameters.create(self.model_name, **hamiltonian_args)
        
        self.geometry_type = geometry_args["Type"]
        self.geometry = Geometry.create(self.geometry_type, **geometry_args)

        self.interaction_type = interaction_args["Type"]
        self.interactions = Interactions(self.interaction_type, self.geometry, **interaction_args)

    
    def compute_h_B(self):

        h = self.hamiltonian_parameters.h
        N = self.geometry.N
        J = self.hamiltonian_parameters.J

        return h/(J*(N-1))