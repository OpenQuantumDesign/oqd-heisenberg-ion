from heisenberg_ion.simulators.ed.preprocessor import ExactDiagonalization
from heisenberg_ion.simulators.qmc.long_range.preprocessor import LongRangeQMC
from heisenberg_ion.simulators.qmc.nearest_neighbor.preprocessor import NearestNeighborQMC


class PreprocessorFactory:
    """
    Factory for generating the required instance of the Preprocessor subclass. Carries a registry of Preprocessor subclasses

    Raises:
        Exception: if requested subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        registers the specified subclass with the given name in the factory

        Args:
            name (str): name to be used while registering the specified subclass
            subclass (Type[Preprocessor]): the subclass of preprocessor to be registered
        """

        cls.registry[name] = subclass

    def create(cls, name, parameter_set_list):
        """
        creates an instance of the requested subclass

        Args:
            name (str): requested subclass name
            parameter_set_list (list[dict]): list of parameter sets, each set specified as a dict, needed by the preprocessor constructor

        Raises:
            Exception: if requested Preprocessor subclass is not found

        Returns:
            (Preprocessor): instance of Preprocessor subclass requested
        """

        if name not in cls.registry:
            raise Exception(f"Preprocessor implementation not found for name: {name}")
        else:
            return cls.registry[name](parameter_set_list)


PreprocessorFactory.create = classmethod(PreprocessorFactory.create)
PreprocessorFactory.register = classmethod(PreprocessorFactory.register)

PreprocessorFactory.register("long_range_qmc", LongRangeQMC)
PreprocessorFactory.register("nearest_neighbor_qmc", NearestNeighborQMC)
PreprocessorFactory.register("exact_diagonalization", ExactDiagonalization)
