from deterministic import Deterministic
from heatbath import Heatbath
from directed_loops import DirectedLoops

class ProbabilityTableFactory:

    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def extract_args(cls, name, **kwargs):

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[key])

        return arg_vals

    def create(cls, name, system, **kwargs):

        if name not in cls.registry:
            raise ValueError(f"ProbabilityTable implementation not found for name: {name}")
        else:
            return cls.registry[name](system, **kwargs)
        

# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically
ProbabilityTableFactory.register = classmethod(ProbabilityTableFactory.register)
ProbabilityTableFactory.create = classmethod(ProbabilityTableFactory.create)


ProbabilityTableFactory.register("Deterministic", Deterministic)
ProbabilityTableFactory.register("Heatbath", Heatbath)
ProbabilityTableFactory.register("DirectedLoops", DirectedLoops)
