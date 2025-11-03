from geometry import Geometry
from hamiltonian import HamiltonianParameters
from interactions import Interactions

class System:

    def __init__(self, hamiltonian_args, geometry_args, interaction_args):

        self.model_name = hamiltonian_args["model_name"]
        hamiltonian_args = self.extract_args(hamiltonian_args)
        self.hamiltonian_parameters = HamiltonianParameters.create(self.model_name, **hamiltonian_args)
        
        self.geometry_name = geometry_args["geometry_name"]
        geometry_args = self.extract_args(geometry_args)
        self.geometry = Geometry.create(self.geometry_name, **geometry_args)

        self.interaction_name = interaction_args["interaction_type"]
        interaction_args = self.extract_args(interaction_args)
        self.interactions = Interactions(self.interaction_type, self.geometry, **interaction_args)

    
    def extract_args(self, kwargs):
        
        args = {key : kwargs[key] for key in kwargs if key != "Name"}
        return args

    
    def compute_h_B(self):

        h = self.hamiltonian_parameters.h
        J = self.hamiltonian_parameters.J
        num_neighbors = self.goemetry.num_neighbors_per_site

        return h/(J*num_neighbors)