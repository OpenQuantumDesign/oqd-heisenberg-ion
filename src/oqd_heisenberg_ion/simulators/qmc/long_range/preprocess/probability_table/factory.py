from .deterministic import Deterministic
from .directed_loops import DirectedLoops
from .heatbath import Heatbath


class ProbabilityTableFactory:
    """
    Factory for generating the required instance of the ProbabilityTable subclass. Carries a registry of ProbabilityTable subclasses
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[ProbabilityTable]): ProbabilityTable subclass to be registered
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

    def create(cls, name, system, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            name (str): name of requested subclass

        Raises:
            Exception: if the requested ProbabilityTable is not found in the registry

        Returns:
            (ProbabilityTable): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"Probability Table implementation not found for name: {name}")
        else:
            return cls.registry[name](system, **kwargs)


# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically
ProbabilityTableFactory.register = classmethod(ProbabilityTableFactory.register)
ProbabilityTableFactory.create = classmethod(ProbabilityTableFactory.create)
ProbabilityTableFactory.extract_args = classmethod(ProbabilityTableFactory.extract_args)

ProbabilityTableFactory.register("deterministic", Deterministic)
ProbabilityTableFactory.register("heatbath", Heatbath)
ProbabilityTableFactory.register("directed_loop", DirectedLoops)
