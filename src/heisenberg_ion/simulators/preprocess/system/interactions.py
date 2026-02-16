import shutil as sh

import numpy as np


class InteractionsFactory:
    """
    Factory for generating the required instance of the Interactions subclass. Carries a registry of Interactions subclasses

    Raises:
        Exception: if requested subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[Interactions]): Interactions subclass to be registered
        """

        cls.registry[name] = subclass

    def extract_args(cls, name, **kwargs):
        """
        extracts the arguments associated with a given subclass

        Args:
            name (str): subclass name (must exist in registry)
            **kwargs (dict): key word arguments. Must contain inputs for the specific subclass

        Returns:
            (dict): contains the subclass arguments as key value pairs
        """

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[key])

        return arg_vals

    def create(cls, name, geometry, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            name (str): name of requested subclass

        Raises:
            Exception: if the requested Interactions subclass is not found in the registry

        Returns:
            (Interactions): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"Interactions implementation not found for interaction name: {name}")
        else:
            return cls.registry[name](geometry, **kwargs)


InteractionsFactory.register = classmethod(InteractionsFactory.register)
InteractionsFactory.create = classmethod(InteractionsFactory.create)
InteractionsFactory.extract_args = classmethod(InteractionsFactory.extract_args)


class Interactions:
    """
    Interactions base class. Different types of interactions implemented as subclasses
    """

    def __init__(self, geometry):
        """
        constructor initializes the member variables

        Args:
            geometry (Geometry): contains contains the lattice sites and distances
        """

        self.geometry = geometry
        self.J_ij_matrix = None
        self.J_ij_vector = None
        self.J_ij_file = None

    def get_J_ij(self):
        """
        every Interactions subclass must implement a method to populate the interaction matrix
        """

        pass


class MatrixInput(Interactions):
    """
    Used when the interactions are specified via an input file containing an coupling matrix
    """

    args = {"interaction_matrix_file": str}

    def __init__(self, geometry, interaction_matrix_file):
        """
        populates the J_ij vector from file

        Args:
            geometry (Geometry): contains the lattice sites and distances
            interaction_matrix_file (str): file path to the interaction matrix file
        """

        super().__init__(geometry)

        self.J_ij_file = interaction_matrix_file
        geometry.initialize_tables()
        self.get_J_ij(geometry.N, geometry.num_bonds, self.J_ij_file)

    def get_J_ij(self, N, num_bonds, interaction_matrix_file):
        """
        populates the J_ij vector with the coupling strength for each bond

        Args:
            N (int): number of sites in the lattice
            num_bonds (int): number of bonds in the lattice
            interaction_matrix_file (str): file path to the interactions matrix
        """

        self.J_ij_matrix = np.loadtxt(interaction_matrix_file, delimiter=",", skiprows=1)
        self.J_ij_vector = np.zeros(num_bonds)
        b = 0
        for i in range(N):
            for j in range(i + 1, N):
                self.J_ij_vector[b] = self.J_ij_matrix[i, j]
                b += 1

    def write_to_file(self, target_file):
        """
        copies the input J_ij matrix to specified directory

        Args:
            target_file (str): path to target file
        """

        sh.copyfile(self.J_ij_file, target_file)


class PowerLaw(Interactions):
    """
    Used to generate the coupling using power law J_ij = 1/r_{ij}^alpha
    """

    args = {"alpha": float}

    def __init__(self, geometry, alpha):
        """
        constructs the coupling matrix

        Args:
            geometry (Geometry): contains the lattice sites and distances
            alpha (float): power law interaction strength exponent
        """

        super().__init__(geometry)

        self.alpha = alpha
        self.geometry.initialize_tables()
        self.geometry.build()
        self.get_J_ij(self.geometry.num_bonds, self.geometry.distances, self.alpha)

    def get_J_ij(self, num_bonds, distances, alpha):
        """
        populates the J_ij matrix

        Args:
            num_bonds (int): number of bonds
            distances (numpy.ndarray[float]): num_bonds x 1 array containing the distances between all pairs of interacting sites
            alpha (float): interaction strength exponent
        """

        N = self.geometry.N
        self.J_ij_vector = np.zeros(num_bonds)
        self.J_ij_matrix = np.zeros((N, N))
        b = 0
        for i in range(N):
            for j in range(i + 1, N):
                r_pow_alpha = (distances[b]) ** alpha
                self.J_ij_matrix[i, j] = 1.0 / r_pow_alpha
                self.J_ij_matrix[j, i] = 1.0 / r_pow_alpha
                self.J_ij_vector[b] = 1.0 / r_pow_alpha
                b += 1

    def write_to_file(self, target_file):
        """
        writes the J_ij matrix to file

        Args:
            target_file (str): path to target file
        """

        np.savetxt(
            target_file,
            self.J_ij_matrix,
            delimiter=",",
            header=f"J_ij_matrix, interactions=power_law,alpha={self.alpha}",
        )


class NearestNeighbors(Interactions):
    """
    interactions for a nearest neighbor lattice
    """

    args = {}

    def __init__(self, geometry):
        """
        calls the base class constructor, no interaction matrix is needed

        Args:
            geometry (Geometry): contains the lattice sites and distances
        """

        super().__init__(geometry)


InteractionsFactory.register("power_law", PowerLaw)
InteractionsFactory.register("matrix_input", MatrixInput)
InteractionsFactory.register("nearest_neighbor", NearestNeighbors)
