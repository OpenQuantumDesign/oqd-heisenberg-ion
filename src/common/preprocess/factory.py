from simulators.qmc.long_range.preprocess.preprocessor import LongRangeQMC
from simulators.qmc.nearest_neighbor.preprocessor import NearestNeighborQMC
from simulators.ed.preprocessor import ExactDiagonalization

class PreprocessFactory:

    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def create(cls, name, parameter_set_list):

        if name not in cls.registry:
            raise ValueError(f"Preprocessor implementation not found for name: {name}")
        else:
            return cls.registry[name](parameter_set_list)
        

PreprocessFactory.create = classmethod(PreprocessFactory.create)
PreprocessFactory.register = classmethod(PreprocessFactory.register)

PreprocessFactory.register("LongRangeQMC", LongRangeQMC)
PreprocessFactory.register("NearestNeighborQMC", NearestNeighborQMC)
PreprocessFactory.register("ExactDiagonalization", ExactDiagonalization)