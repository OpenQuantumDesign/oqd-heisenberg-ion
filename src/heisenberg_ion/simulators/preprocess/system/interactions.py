import shutil as sh

import numpy as np


# Make register and create factory classmethods so subclasses can be added to the registry
# and instantiated agnostically
class InteractionsFactory:
    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def extract_args(cls, name, **kwargs):

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[key])

        return arg_vals

    def create(cls, name, geometry, **kwargs):

        if name not in cls.registry:
            raise ValueError(f"Geometry implementation not found for geometry name: {name}")
        else:
            return cls.registry[name](geometry, **kwargs)


InteractionsFactory.register = classmethod(InteractionsFactory.register)
InteractionsFactory.create = classmethod(InteractionsFactory.create)
InteractionsFactory.extract_args = classmethod(InteractionsFactory.extract_args)


class Interactions:
    def __init__(self, geometry):

        self.geometry = geometry
        self.J_ij_matrix = None
        self.J_ij_vector = None
        self.J_ij_file = None

    def get_J_ij(self):

        pass


class MatrixInput(Interactions):
    args = {"interaction_matrix_file": str}

    def __init__(self, geometry, interaction_matrix_file):

        super().__init__(geometry)

        self.J_ij_file = interaction_matrix_file
        geometry.initialize_tables()
        self.get_J_ij(geometry.N, geometry.num_bonds, self.J_ij_file)

    def get_J_ij(self, N, num_bonds, interaction_matrix_file):

        self.J_ij_matrix = np.loadtxt(interaction_matrix_file, delimiter=",", skiprows=1)
        self.J_ij_vector = np.zeros(num_bonds)
        b = 0
        for i in range(N):
            for j in range(i + 1, N):
                self.J_ij_vector[b] = self.J_ij_matrix[i, j]
                b += 1

        return 0

    def write_to_file(self, target_file):

        sh.copyfile(self.J_ij_file, target_file)

        return 0


class PowerLaw(Interactions):
    args = {"alpha": float}

    def __init__(self, geometry, alpha):

        super().__init__(geometry)

        self.alpha = alpha
        self.geometry.initialize_tables()
        self.geometry.build()
        self.get_J_ij(self.geometry.num_bonds, self.geometry.distances, self.alpha)

    def get_J_ij(self, num_bonds, distances, alpha):

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

        return 0

    def write_to_file(self, target_file):

        np.savetxt(
            target_file,
            self.J_ij_matrix,
            delimiter=",",
            header=f"J_ij_matrix, interactions=power_law,alpha={self.alpha}",
        )


class NearestNeighbors(Interactions):
    args = {}

    def __init__(self, geometry):

        super().__init__(geometry)


InteractionsFactory.register("power_law", PowerLaw)
InteractionsFactory.register("matrix_input", MatrixInput)
InteractionsFactory.register("nearest_neighbor", NearestNeighbors)
