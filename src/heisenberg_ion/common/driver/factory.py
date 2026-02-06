from heisenberg_ion.simulators.ed.driver import ExactDiagonalization
from heisenberg_ion.simulators.qmc.long_range.driver import LongRangeQMC
from heisenberg_ion.simulators.qmc.nearest_neighbor.driver import NearestNeighborQMC


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


DriverFactory.register("long_range_qmc", LongRangeQMC)
DriverFactory.register("nearest_neighbor_qmc", NearestNeighborQMC)
DriverFactory.register("exact_diagonalization", ExactDiagonalization)
