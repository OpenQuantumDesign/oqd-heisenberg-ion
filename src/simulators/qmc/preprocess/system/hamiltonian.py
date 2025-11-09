class HamiltonianParameters:

    args = {}

    def __init__(self):

        self.hamiltonian_name = None
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
HamiltonianFactory.extract_args = classmethod(HamiltonianFactory.extract_args)


class FMHeisenbergAFMZ(HamiltonianParameters):

    args = {"J" : float}

    def __init__(self, J):

        super.__init__()

        self.hamiltonian_name = "FMHeisenbergAFMZ"
        self.hamiltonian_type = -1
        self.Delta = -1.0
        self.h = 0.0
        self.J = J
        
class XY(HamiltonianParameters):

    args = {"J" : float}

    def __init__(self, J):

        super.__init__()

        self.hamiltonian_name = "XY"
        self.hamiltonian_type = 0
        self.Delta = 0.0
        self.h = 0.0
        self.J = J

class FMHeisenbergFMZ(HamiltonianParameters):

    args = {"J" : float}

    def __init__(self, J):

        super.__init__()

        self.hamiltonian_name = "FMHeisenbergFMZ"
        self.hamiltonian_type = 1
        self.Delta = 1.0
        self.h = 0.0
        self.J = J

class XXZ(HamiltonianParameters):

    args = {"Delta": float, "J" : float}

    def __init__(self, Delta, J):

        super.__init__()
        
        self.Delta = Delta
        self.hamiltonian_type = 2
        self.h = 0.0
        self.J = J

class XXZh(HamiltonianParameters):

    args = {"Delta": float, "h": float, "J" : float}

    def __init__(self, Delta, h, J):

        super.__init__()

        self.Delta = Delta
        self.hamiltonian_type = 3
        self.h = h
        self.J = J

class AFMHeisenbergFMZ(HamiltonianParameters):

    args = {"J" : float}

    def __init__(self, J):

        super.__init__()

        self.hamiltonian_name = "AFMHeisenbergFMZ"
        self.hamiltonian_type = 4
        self.Delta = 1.0
        self.h = 0.0
        self.J = J

# Register all implemented Hamiltonian parameter types in the base HamiltonianParameters class
HamiltonianFactory.register("XY", XY)
HamiltonianFactory.register("XXZ", XXZ)
HamiltonianFactory.register("XXZh", XXZh)
HamiltonianFactory.register("FMHeisenbergFMZ", FMHeisenbergFMZ)
HamiltonianFactory.register("FMHeisenbergAFMZ", FMHeisenbergAFMZ)
HamiltonianFactory.register("AFMHeisenbergFMZ", AFMHeisenbergFMZ)