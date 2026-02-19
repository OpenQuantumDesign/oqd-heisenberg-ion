from oqd_heisenberg_ion.common.driver.base import Driver

from .engine.configuration_generator import ConfigurationGenerator


class NearestNeighborQMC(Driver):
    """
    Driver subclass for nearest neighbor qmc
    """

    def __init__(self, simulation_folder, simulator_inputs):
        """
        initializes the nearest neighbor qmc driver and counts the number of parameter sets

        Args:
            simulation_folder (str): path to the simulation output folder
            simulator_inputs (list[dict]): list of nearest neighbor qmc parameter sets to be simulated
        """

        super().__init__(simulation_folder)

        self.sse_inputs = simulator_inputs

        self.num_parameter_sets = len(simulator_inputs)

    def simulate(self):
        """
        calls the Python ConfigurationGenerator class, the nearest neighbor qmc engine

        Returns:
            (int): exit code, 0 if the simulation terminates gracefully
        """

        for i in range(self.num_parameter_sets):
            run_folder = self.sse_inputs[i].run_folder
            system = self.sse_inputs[i].system
            sampling_args = self.sse_inputs[i].sampling_parameters

            config_generator = ConfigurationGenerator(system, sampling_args, run_folder)
            config_generator.simulate()
            config_generator.write_outputs()

        return 0
