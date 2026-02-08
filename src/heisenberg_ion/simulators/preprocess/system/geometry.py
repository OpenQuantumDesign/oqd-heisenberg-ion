import numpy as np


# ---------------------------------- Base Geometry Class  ----------------------------------
class Geometry:
    args = {}

    def __init__(self):

        self.spatial_dimension = None

        self.boundary = None
        self.interaction_range = None
        self.lattice_type = None

        self.N_config = None
        self.N = None
        self.num_bonds = None
        self.num_neighbors_per_site = None

        self.sites = None
        self.distances = None
        self.geometry_table = None

        self.bipartite = False

    def initialize_tables(self):

        self.geometry_table = np.zeros((self.num_bonds, 3))
        self.distances = np.zeros(self.num_bonds)
        self.sites = np.zeros((self.num_bonds, 2), dtype=int)

        return 0

    def build(self):
        pass

    def update_parameters(self, parameter_dict):

        parameter_dict["spatial_dimension"] = self.spatial_dimension
        parameter_dict["boundary"] = self.boundary
        parameter_dict["interaction_range"] = self.interaction_range
        parameter_dict["N"] = self.N
        parameter_dict["num_bonds"] = self.num_bonds

        return parameter_dict


# Make register and create factory classmethods so subclasses can be added to the registry
# and instantiated agnostically
class GeometryFactory:
    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def extract_args(cls, name, **kwargs):

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[key])

        return arg_vals

    def create(cls, name, **kwargs):

        if name not in cls.registry:
            raise ValueError(f"Geometry implementation not found for geometry name: {name}")
        else:
            return cls.registry[name](**kwargs)


GeometryFactory.register = classmethod(GeometryFactory.register)
GeometryFactory.create = classmethod(GeometryFactory.create)
GeometryFactory.extract_args = classmethod(GeometryFactory.extract_args)


# ------------------ Long Range 1d Chain With Open Boundaries ------------------
class LongRangeOpenChain(Geometry):
    args = {"N": int}

    def __init__(self, N):

        super().__init__()

        self.spatial_dimension = 1

        self.boundary = "open"
        self.interaction_range = "long_range"
        self.lattice_type = "chain"

        self.N = N
        self.N_config = N
        self.num_bonds = int(self.N * (self.N - 1) / 2)
        self.num_neighbors_per_site = N - 1

        # self.initialize_tables()
        # self.build()

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        self.initialize_tables()

        b = 0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                self.sites[b, 0] = i
                self.sites[b, 1] = j

                self.distances[b] = j - i

                self.geometry_table[b, 0] = i
                self.geometry_table[b, 1] = j
                self.geometry_table[b, 2] = self.distances[b]

                b += 1

        return 0


# ------------------ Long Range 1d Chain With Periodic Boundaries ------------------
class LongRangePeriodicChain(Geometry):
    args = {"N": int}

    def __init__(self, N):

        super().__init__()

        self.spatial_dimension = 1

        self.boundary = "periodic"
        self.interaction_range = "long_range"
        self.lattice_type = "chain"

        self.N = N
        self.N_config = N
        self.num_bonds = int(self.N * (self.N - 1) / 2)
        self.num_neighbors_per_site = N - 1

        # self.initialize_tables()
        # self.build()

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        self.initialize_tables()

        b = 0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                self.sites[b, 0] = i
                self.sites[b, 1] = j

                self.geometry_table[b, 0] = i
                self.geometry_table[b, 1] = j

                if (j - i) <= self.N - (j - i):
                    self.distances[b] = j - i
                    self.geometry_table[b, 2] = self.distances[b]
                else:
                    self.distances[b] = self.N - (j - i)
                    self.geometry_table[b, 2] = -self.distances[b]

                b += 1

        return 0


# ------------------ Long Range 2d Triangular Lattice With Open Boundaries ------------------
class LongRangeOpenTriangular(Geometry):
    args = {"N1": int, "N2": int}

    def __init__(self, N1, N2):

        super().__init__()

        self.spatial_dimension = 2

        self.boundary = "periodic"
        self.interaction_range = "long_range"
        self.lattice_type = "chain"

        self.N = N1 * N2
        self.N_config = (N1, N2)
        self.num_bonds = int(self.N * (self.N - 1) / 2)
        self.num_neighbors_per_site = self.N - 1

        # self.initialize_tables()
        # self.build()

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        self.initialize_tables()

        a_1 = np.array([1.0, 0.0])
        a_2 = np.array([-0.5, np.sqrt(3.0) / 2.0])

        N_1 = self.N_config[0]
        N_2 = self.N_config[1]

        b = 0
        for i1 in range(N_1):
            for i2 in range(N_2):
                i = i1 * N_2 + i2

                for j1 in range(N_1):
                    for j2 in range(N_2):
                        j = j1 * N_2 + j2

                        if j > i:
                            self.sites[b, 0] = i
                            self.sites[b, 1] = j

                            self.distances[b] = np.norm((j1 - i1) * a_1 + (j2 - i2) * a_2)

                            self.geometry_table[b, 0] = i
                            self.geometry_table[b, 1] = j
                            self.geometry_table[b, 2] = self.distances[b]

                            b += 1

        return 0


class NearestNeighborPeriodicChain(Geometry):
    args = {"N": int}

    def __init__(self, N):

        self.spatial_dimension = 1

        self.boundary = "periodic"
        self.interaction_range = "nearest_neighbor"
        self.lattice_type = "chain"

        self.N_config = N
        self.N = N
        self.num_bonds = N
        self.num_neighbors_per_site = 2

        if self.N % 2 == 0:
            self.bipartite = True
        else:
            self.bipartite = False

        # self.initialize_tables()
        # self.build()

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        self.initialize_tables()

        self.sites = np.zeros((self.num_bonds, 2), dtype=int)

        for i in range(self.N - 1):
            self.sites[i, 0] = i
            self.sites[i, 1] = i + 1

        self.sites[self.N - 1, 0] = 0
        self.sites[self.N - 1, 1] = self.N - 1


class NearestNeighborOpenChain(Geometry):
    args = {"N": int}

    def __init__(self, N):

        self.spatial_dimension = 1

        self.boundary = "periodic"
        self.interaction_range = "nearest_neighbor"
        self.lattice_type = "chain"

        self.N_config = N
        self.N = N
        self.num_bonds = N
        self.num_neighbors_per_site = 2

        self.bipartite = True

        # self.initialize_tables()
        # self.build()

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        self.initialize_tables()

        self.sites = np.zeros((self.num_bonds, 2), dtype=int)

        for i in range(self.N - 1):
            self.sites[i, 0] = i
            self.sites[i, 1] = i + 1


# Register all implemented geometries in the geometry factory
GeometryFactory.register("long_range_open_1d_rectangular", LongRangeOpenChain)
GeometryFactory.register("long_range_periodic_1d_rectangular", LongRangePeriodicChain)
GeometryFactory.register("long_range_open_2d_triangular", LongRangeOpenTriangular)
GeometryFactory.register("nearest_neighbor_periodic_1d_rectangular", NearestNeighborPeriodicChain)
GeometryFactory.register("nearest_neighbor_open_1d_rectangular", NearestNeighborOpenChain)
