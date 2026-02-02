from .geometry import GeometryFactory
from .hamiltonian import HamiltonianFactory
from .interactions import InteractionsFactory

class System:

    def __init__(self, **kwargs):

        self.model_name = kwargs["hamiltonian_name"]
        hamiltonian_args = HamiltonianFactory.extract_args(self.model_name, **kwargs)
        self.hamiltonian_parameters = HamiltonianFactory.create(self.model_name, **hamiltonian_args)
        
        self.geometry_name = kwargs["interaction_range"] + "_" + kwargs["boundary"] + "_"\
            + kwargs["spatial_dimension"] + "_" + kwargs["lattice_type"]
        geometry_args = GeometryFactory.extract_args(self.geometry_name, **kwargs)
        self.geometry = GeometryFactory.create(self.geometry_name, **geometry_args)

        self.interaction_range = kwargs["interaction_range"]
        self.interaction_name = self.get_interaction_name(kwargs)
        self.interaction_args = InteractionsFactory.extract_args(self.interaction_name, **kwargs)
        self.interactions = InteractionsFactory.create(self.interaction_name, self.geometry, **kwargs)


    def get_interaction_name(self, kwarg_dict):

        if self.interaction_range == "long_range":
            return kwarg_dict["interaction_range"]
        else:
            return self.interaction_range
        
    
    def compute_h_B(self):

        h = self.hamiltonian_parameters.h
        J = self.hamiltonian_parameters.J
        num_neighbors = self.goemetry.num_neighbors_per_site

        return h/(J*num_neighbors)

