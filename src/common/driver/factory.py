from simulators.qmc.long_range.driver import LongRangeQMC
from simulators.qmc.nearest_neighbor.driver import NearestNeighborQMC
from simulators.ed.driver import ExactDiagonalization

class DriverFactory:

    registry = {}

    def register(cls, name, subclass):

        cls.registry[name] = subclass
        

    def create(cls, name, root_folder, simulator_inputs):

        if name not in cls.registry:
            raise ValueError(f"Driver implementation not found for name: {name}")
        else:
            return cls.registry[name](root_folder, simulator_inputs)
        

DriverFactory.create = classmethod(DriverFactory.create)
DriverFactory.register = classmethod(DriverFactory.register)


DriverFactory.register("LongRangeQMC", LongRangeQMC)
DriverFactory.register("NearestNeighborQMC", NearestNeighborQMC)
DriverFactory.register("ExactDiagonalization", ExactDiagonalization)