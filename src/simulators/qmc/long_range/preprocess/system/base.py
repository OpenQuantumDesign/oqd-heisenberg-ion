from geometry import GeometryFactory
from hamiltonian import HamiltonianFactory
from interactions import Interactions

class System:

    def __init__(self, **kwargs):

        self.model_name = kwargs["HamiltonianType"]

        hamiltonian_args = HamiltonianFactory.extract_args(self.model_name, **kwargs)
        self.hamiltonian_parameters = HamiltonianFactory.create(self.model_name, **hamiltonian_args)
        
        self.geometry_name = kwargs["InteractionRange"] + kwargs["BoundaryType"] \
            + kwargs["SpatialDimension"] + geometry_args["GeometryType"]
        
        geometry_args = GeometryFactory.extract_args(self.geometry_name, **kwargs)
        self.geometry = GeometryFactory.create(self.geometry_name, **geometry_args)

        self.interaction_type = kwargs["InteractionType"]
        
        self.interactions = Interactions(self.geometry, **kwargs)
    
    def compute_h_B(self):

        h = self.hamiltonian_parameters.h
        J = self.hamiltonian_parameters.J
        num_neighbors = self.goemetry.num_neighbors_per_site

        return h/(J*num_neighbors)