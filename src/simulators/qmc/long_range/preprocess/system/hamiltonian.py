class HamiltonianParameters:

    args = {}

    def __init__(self):

        self.hamiltonian_type = None
        self.Delta = None
        self.h = None
        self.J = None

class HamiltonianFactory:

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
            raise ValueError(f"HamiltonianParameters implementation not found for name: {name}")
        else:
            return cls.registry[name](**kwargs)
        
# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically
HamiltonianFactory.register = classmethod(HamiltonianFactory.register)
HamiltonianFactory.create = classmethod(HamiltonianFactory.create)
        
class XY(HamiltonianParameters):

    args = {"hamiltonian_type" : int, "J" : float}

    def __init__(self, hamiltonian_type, J):

        super.__init__()

        if hamiltonian_type != 0:
            raise ValueError("hamiltonian_type must be 0 for the XY model")

        self.hamiltonian_type = 0
        self.Delta = 0.0
        self.h = 0.0
        self.J = J

class Heisenberg(HamiltonianParameters):

    args = {"hamiltonian_type" : int, "J" : float}

    def __init__(self, hamiltonian_type, J):

        super.__init__()

        if hamiltonian_type != 1 or hamiltonian_type != -1:
            raise ValueError("hamiltonian_type must be 1 or -1 for the Heisenberg model")

        self.hamiltonian_type = hamiltonian_type
        self.Delta = -1.0
        self.h = 0.0
        self.J = J

class XXZ(HamiltonianParameters):

    args = {"hamiltonian_type" : int, "Delta": float, "J" : float}

    def __init__(self, hamiltonian_type, Delta, J):

        super.__init__()

        if hamiltonian_type != 2:
            raise ValueError("hamiltonian_type must be 2 for the XXZ model")
        
        self.Delta = Delta
        self.h = 0.0
        self.J = J

class XXZh(HamiltonianParameters):

    args = {"hamiltonian_type" : int, "Delta": float, "h": float, "J" : float}

    def __init__(self, hamiltonian_type, Delta, h, J):

        super.__init__()

        if hamiltonian_type != 3:
            raise ValueError("hamiltonian_type must be 3 for the XXZh model")

        self.Delta = Delta
        self.h = h
        self.J = J

# Register all implemented Hamiltonian parameter types in the base HamiltonianParameters class
HamiltonianFactory.register("XY", XY)
HamiltonianFactory.register("XXZ", XXZ)
HamiltonianFactory.register("XXZh", XXZh)
HamiltonianFactory.register("Heisenberg", Heisenberg)