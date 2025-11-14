from heisenberg_ion.common.driver.base import Driver
from .engine.configuration_generator import ConfigurationGenerator

class NearestNeighborQMC(Driver):

    def __init__(self, simulation_folder, simulator_inputs):
        
        super().__init__(simulation_folder)

        self.sse_inputs = simulator_inputs

        self.num_parameter_sets = len(simulator_inputs)


    def simulate(self):

        for i in range(self.num_parameter_sets):

            run_folder = self.sse_inputs[i].run_folder
            system = self.sse_inputs[i].system
            sampling_args = self.sse_inputs[i].sampling_parameters

            config_generator = ConfigurationGenerator(system, sampling_args, run_folder)
            config_generator.simulate()
            config_generator.write_outputs()

        return 0
        