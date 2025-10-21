import numpy as np

class Geometry:

    def __init__(self):

        self.spatial_dimension = None

        self.boundary_type = None
        
        self.N_config = None
        self.N = None
        self.num_bonds = None

        self.sites = None
        self.distances = None
        self.geometry_table = None

        return 0

    def initialize_tables(self):

        self.geometry_table = np.zeros((self.num_bonds,3))
        self.distances = np.zeros(self.num_bonds)
        self.sites = np.zeros((self.num_bonds,2), dtype=int)

        return 0

    def build(self):
        pass


class LongRangeOpenChain(Geometry):

    def __init__(self, N):

        super.__init__()

        self.spatial_dimension = 1

        self.boundary_type = "Open"

        self.N = N
        self.N_config = (N)
        self.num_bonds = int(self.N*(self.N-1)/2)

        self.initialize_tables()

        return 0

    def initialize_tables(self):
        return super().initialize_tables()

    def build(self):

        b = 0
        for i in range(self.N):
            for j in range(i+1,self.N):

                self.sites[b,0] = i
                self.sites[b,1] = j

                self.distances[b] = j-i

                self.geometry_table[b,0] = i
                self.geometry_table[b,1] = j
                self.geometry_table[b,2] = self.distances[b]

                b += 1
        
        return 0


class LongRangePeriodicChain(Geometry):

    def __init__(self, N):

        super.__init__()

        self.spatial_dimension = 1

        self.boundary_type = "Periodic"

        self.N = N
        self.N_config = (N)
        self.num_bonds = int(self.N*(self.N-1)/2)

        self.initialize_tables()

        return 0

    def initialize_tables(self):
        return super().initialize_tables()
    
    def build(self):

        b = 0
        for i in range(self.N):
            for j in range(i+1,self.N):

                self.sites[b,0] = i
                self.sites[b,1] = j

                self.geometry_table[b,0] = i
                self.geometry_table[b,1] = j

                if (j-i) <= self.N - (j-i):
                    self.distances[b] = j-i
                    self.geometry_table[b,2] = self.distances[b]
                else:
                    self.distances[b] = N - (j-i)
                    self.geometry_table[b,2] = -self.distances[b]

                b += 1

        return 0


class LongRangeOpenTriangular(Geometry):

    def __init__(self, N_1, N_2):

        super.__init__()

        self.spatial_dimension = 2

        self.boundary_type = "Periodic"

        self.N = N_1 * N_2
        self.N_config = (N_1, N_2)
        self.num_bonds = int(self.N*(self.N-1)/2)

        self.initialize_tables()

        return 0

    def initialize_tables(self):
        return super().initialize_tables()
    
    def build(self):

        a_1 = np.array([1.0,0.0])
        a_2 = np.array([-0.5, np.sqrt(3.0)/2.0])
        
        N_1 = self.N_config[0]
        N_2 = self.N_config[1]

        for i1 in range(N_1):
            for i2 in range(N_2):
                i = i1*N_2 + i2

                for j1 in range(N_1):
                    for j2 in range(N_2):
                        j = j1*N_2 + j2

                        if j > i:
                            self.sites[b,0] = i
                            self.sites[b,1] = j

                            self.distances[b] = np.norm((j1-i1)*a_1 + (j2-i2)*a_2)

                            self.geometry_table[b,0] = i
                            self.geometry_table[b,1] = j
                            self.geometry_table[b,2] = self.distances[b]

                            b += 1

        return 0
    