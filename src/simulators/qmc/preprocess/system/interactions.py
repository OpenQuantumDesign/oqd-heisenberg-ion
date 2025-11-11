import numpy as np

# Make register and create factory classmethods so subclasses can be added to the registry and instantiated agnostically
class InteractionsFactory:

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
        self.J_ij = None

    def get_J_ij_vector(self):

        pass


class MatrixInput(Interactions):

    args = {"JijFile": str}

    def __init__(self, geometry, JijFile):

        super().__init__(geometry)

        self.file_path = JijFile
        self.get_J_ij(geometry.N, geometry.num_bonds, self.file_path)
        

    def get_J_ij(self, N, num_bonds, J_ij_file):
    
        J_ij_matrix = np.loadtxt(J_ij_file, delimiter=',', skiprows=1)
        self.J_ij_vector = np.zeros(num_bonds)
        b = 0
        for i in range(N):
            for j in range(i+1,N):
                self.J_ij_vector[b] = J_ij_matrix[i,j]
                b += 1

        return 0
    

class PowerLaw(Interactions):

    args = {"Alpha" : float}

    def __init__(self, geometry, Alpha):

        super().__init__(geometry)

        self.alpha = Alpha
        self.get_J_ij(geometry.num_bonds, geometry.distances, self.alpha)


    def get_J_ij(self, num_bonds, distances, alpha):

        self.J_ij_vector = np.zeros(num_bonds)
        for b in range(num_bonds):
            self.J_ij_vector[b] = 1.0/(distances[b])**alpha

        return 0
    

class NearestNeighbors(Interactions):

    args = {}

    def __init__(self, geometry):

        super().__init__(geometry)
        

InteractionsFactory.register("PowerLaw", PowerLaw)
InteractionsFactory.register("MatrixInput", MatrixInput)
InteractionsFactory.register("NearestNeighbors", NearestNeighbors)
