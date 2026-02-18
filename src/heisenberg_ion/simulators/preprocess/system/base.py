from .geometry import GeometryFactory
from .hamiltonian import HamiltonianFactory
from .interactions import InteractionsFactory


class System:
    """
    Defines the model. This includes the Hamiltonian parameters, interaction strengths and the lattice geometry
    """

    def __init__(self, **kwargs):
        """
        Defines the system from key word arguments

        Args:
            **kwargs(dict): key word arguments required to specify the system
        """

        self.model_name = kwargs["hamiltonian_name"]
        hamiltonian_args = HamiltonianFactory.extract_args(self.model_name, **kwargs)
        self.hamiltonian_parameters = HamiltonianFactory.create(self.model_name, **hamiltonian_args)

        self.geometry_name = (
            kwargs["interaction_range"]
            + "_"
            + kwargs["boundary"]
            + "_"
            + kwargs["spatial_dimension"]
            + "_"
            + kwargs["lattice_type"]
        )
        geometry_args = GeometryFactory.extract_args(self.geometry_name, **kwargs)
        self.geometry = GeometryFactory.create(self.geometry_name, **geometry_args)

        self.interaction_range = kwargs["interaction_range"]
        self.interaction_name = self.get_interaction_name(kwargs)
        self.interaction_args = InteractionsFactory.extract_args(self.interaction_name, **kwargs)
        self.interactions = InteractionsFactory.create(self.interaction_name, self.geometry, **self.interaction_args)

    def get_interaction_name(self, kwarg_dict):
        """
        Extracts sets the interaction_name from inputs if long range interactions are specified by inputs

        Args:
            kwarg_dict (dict): contains the interaction_type key value pair

        Returns:
            (str): interaction name
        """

        if self.interaction_range == "long_range":
            return kwarg_dict["interaction_type"]
        else:
            return self.interaction_range

    def compute_h_B(self):
        """
        computes the required h_B parameter for stochastic series expansion based on the field strength, the energy scale and the geometry

        Returns:
            (float): h_B, the field contribution of a single bond in SSE
        """

        h = self.hamiltonian_parameters.h
        J = self.hamiltonian_parameters.J
        num_neighbors = self.geometry.num_neighbors_per_site

        return h / (J * num_neighbors)

    def update_parameters(self, parameter_dict):
        """
        updates a dictionary with system parameters

        Args:
            parameter_dict (dict): parameter set specified as key word arguments

        Returns:
            (dict): updated parameter set
        """

        self.hamiltonian_parameters.update_parameters(parameter_dict)
        self.geometry.update_parameters(parameter_dict)

        return parameter_dict
