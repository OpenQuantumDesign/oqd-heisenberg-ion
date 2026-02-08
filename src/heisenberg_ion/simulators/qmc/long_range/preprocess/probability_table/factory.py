from .deterministic import Deterministic
from .directed_loops import DirectedLoops
from .heatbath import Heatbath


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
            raise ValueError(f"Probability Table implementation not found for name: {name}")
        else:
            return cls.registry[name](system, **kwargs)


# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically
ProbabilityTableFactory.register = classmethod(ProbabilityTableFactory.register)
ProbabilityTableFactory.create = classmethod(ProbabilityTableFactory.create)
ProbabilityTableFactory.extract_args = classmethod(ProbabilityTableFactory.extract_args)

ProbabilityTableFactory.register("deterministic", Deterministic)
ProbabilityTableFactory.register("heatbath", Heatbath)
ProbabilityTableFactory.register("directed_loop", DirectedLoops)
