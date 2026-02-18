from heisenberg_ion.simulators.ed.driver import ExactDiagonalization
from heisenberg_ion.simulators.qmc.long_range.driver import LongRangeQMC
from heisenberg_ion.simulators.qmc.nearest_neighbor.driver import NearestNeighborQMC


class DriverFactory:
    """
     Generates a Driver object based on input specifications

     Carries a registry of recognized Driver subclasses, the constructors of which are called to build the Driver objects

    Raises:
        ValueError: if the requested driver implementation is not found

    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds a Driver subclass to the registry

        Args:
            name (str): defines the name of the subclass to be registered
            subclass (Type[Driver]): specifies the subclass to be registered
        """

        cls.registry[name] = subclass

    def create(cls, name, root_folder, simulator_inputs):
        """
        generates an instance of the requested subclass of Driver

        Args:
            name (str): specifies the name of the requested subclass
            root_folder (str): root directory for simulation outputs
            simulator_inputs (dict): contains the key word arguments to instantiate the Driver subclass requested

        Raises:
            Exception: if requested driver not found in registry

        Returns:
            (Driver): instance of requested driver subclass
        """

        if name not in cls.registry:
            raise Exception(f"Driver implementation not found for name: {name}")
        else:
            return cls.registry[name](root_folder, simulator_inputs)


DriverFactory.create = classmethod(DriverFactory.create)
DriverFactory.register = classmethod(DriverFactory.register)


DriverFactory.register("long_range_qmc", LongRangeQMC)
DriverFactory.register("nearest_neighbor_qmc", NearestNeighborQMC)
DriverFactory.register("exact_diagonalization", ExactDiagonalization)
