from heisenberg_ion.simulators.ed.preprocessor import ExactDiagonalization
from heisenberg_ion.simulators.qmc.long_range.preprocessor import LongRangeQMC
from heisenberg_ion.simulators.qmc.nearest_neighbor.preprocessor import NearestNeighborQMC


class PreprocessorFactory:
    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def create(cls, name, parameter_set_list):

        if name not in cls.registry:
            raise ValueError(f"Preprocessor implementation not found for name: {name}")
        else:
            return cls.registry[name](parameter_set_list)


PreprocessorFactory.create = classmethod(PreprocessorFactory.create)
PreprocessorFactory.register = classmethod(PreprocessorFactory.register)

PreprocessorFactory.register("long_range_qmc", LongRangeQMC)
PreprocessorFactory.register("nearest_neighbor_qmc", NearestNeighborQMC)
PreprocessorFactory.register("exact_diagonalization", ExactDiagonalization)
