from src.simulators.qmc.long_range.preprocessor import LongRangeQMC
from simulators.qmc.nearest_neighbor.preprocessor import NearestNeighborQMC
from simulators.ed.preprocessor import ExactDiagonalization

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

PreprocessorFactory.register("LongRangeQMC", LongRangeQMC)
PreprocessorFactory.register("NearestNeighborQMC", NearestNeighborQMC)
PreprocessorFactory.register("ExactDiagonalization", ExactDiagonalization)