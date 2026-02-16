import numpy as np


# ---------------------------------- Base Geometry Class  ----------------------------------
class Geometry:
    """
    Geometry base class to determine lattice properties
    """

    args = {}

    def __init__(self):
        """
        constructor specifies the member variables that determine the lattice properties
        """

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
        """
        initialize the geometry table, distances and sites as numpy arrays
        """

        self.geometry_table = np.zeros((self.num_bonds, 3))
        self.distances = np.zeros(self.num_bonds)
        self.sites = np.zeros((self.num_bonds, 2), dtype=int)

    def build(self):
        """
        Each Geometry subclass should implement a build function to populate the geometry tables
        """
        pass

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with lattice properties

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        parameter_dict["spatial_dimension"] = self.spatial_dimension
        parameter_dict["boundary"] = self.boundary
        parameter_dict["interaction_range"] = self.interaction_range
        parameter_dict["N"] = self.N
        parameter_dict["num_bonds"] = self.num_bonds

        return parameter_dict


# Make register and create factory classmethods so subclasses can be added to the registry
# and instantiated agnostically
class GeometryFactory:
    """
    Factory for generating the required instance of the Geometry subclass. Carries a registry of Geometry subclasses

    Raises:
        Exception: if requested geometry subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[Geometry]): Geometry subclass to be registered
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

    def create(cls, name, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            name (str): name of requested subclass

        Raises:
            Exception: if the requested Geometry is not found in the registry

        Returns:
            (Geometry): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"Geometry implementation not found for geometry name: {name}")
        else:
            return cls.registry[name](**kwargs)


GeometryFactory.register = classmethod(GeometryFactory.register)
GeometryFactory.create = classmethod(GeometryFactory.create)
GeometryFactory.extract_args = classmethod(GeometryFactory.extract_args)


# ------------------ Long Range 1d Chain With Open Boundaries ------------------
class LongRangeOpenChain(Geometry):
    """
    implements a 1d lattice with long range interactions and open boundaries
    """

    args = {"N": int}

    def __init__(self, N):
        """
        sets the member variables without constructing the tables.

        Args:
            N (int): number of sites
        """

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
        """
        initializes the geometry tables
        """
        return super().initialize_tables()

    def build(self):
        """
        populates the geometry tables
        """

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


# ------------------ Long Range 1d Chain With Periodic Boundaries ------------------
class LongRangePeriodicChain(Geometry):
    """
    implements a 1d lattice with long range interactions and periodic boundaries
    """

    args = {"N": int}

    def __init__(self, N):
        """
        sets the member variables without constructing the tables.

        Args:
            N (int): number of sites
        """

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
        """
        initializes the geometry tables
        """
        return super().initialize_tables()

    def build(self):
        """
        populates the geometry tables
        """

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


# ------------------ Long Range 2d Triangular Lattice With Open Boundaries ------------------
class LongRangeOpenTriangular(Geometry):
    """
    implements a 2d triangular lattice with long range interactions and periodic boundaries (needs to be tested further)
    """

    args = {"N1": int, "N2": int}

    def __init__(self, N1, N2):
        """
        sets the member variables without constructing the tables.

        Args:
            N1 (int): number of sites in the first spatial dimension
            N2 (int): number of sites in the second spatial dimension
        """

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
        """
        initializes the geometry tables
        """
        return super().initialize_tables()

    def build(self):
        """
        populates the geometry tables
        """

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


class NearestNeighborPeriodicChain(Geometry):
    """
    implements a 1d lattice with nearest neighbor interactions and periodic boundaries
    """

    args = {"N": int}

    def __init__(self, N):
        """
        sets the member variables without constructing the tables.

        Args:
            N (int): number of sites
        """

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
        """
        initializes the geometry tables
        """
        return super().initialize_tables()

    def build(self):
        """
        populates the geometry tables
        """

        self.initialize_tables()

        self.sites = np.zeros((self.num_bonds, 2), dtype=int)

        for i in range(self.N - 1):
            self.sites[i, 0] = i
            self.sites[i, 1] = i + 1

        self.sites[self.N - 1, 0] = 0
        self.sites[self.N - 1, 1] = self.N - 1


class NearestNeighborOpenChain(Geometry):
    """
    implements a 1d lattice with nearest neighbor interactions and open boundaries
    """

    args = {"N": int}

    def __init__(self, N):
        """
        sets the member variables without constructing the tables.

        Args:
            N (int): number of sites
        """

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
        """
        initializes the geometry tables
        """
        return super().initialize_tables()

    def build(self):
        """
        populates the geometry tables
        """

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
