from heisenberg_ion.common.preprocessor.base import Preprocessor
from ...preprocess.system.base import System
from heisenberg_ion.common.inputs.input_parser import InputParser
from .preprocess.sampling_params import SamplingParameters
from .preprocess.sampling_params import SSEParameters
import os

class NearestNeighborQMC(Preprocessor):

    allowed_hamiltonians = {"afm_heisenberg_fm_Z", "XY", "fm_heisenberg_fm_Z"}

    def __init__(self, parameter_set_list):

        super().__init__(parameter_set_list)

        self.driver_inputs = []

        self.build()


    def build(self):

        self.check_single_input("root_folder")
        self.root_folder = self.parameter_set_list[0]["root_folder"]

        self.simulation_folder = self.create_output_folder()

        self.check_unique_uuids()

        self.configure_simulation()


    def configure_simulation(self):

        for i in range(self.num_parameter_sets):

            parameter_set = self.parameter_set_list[i]

            self.configure_parameter_set(parameter_set)

            #self.write_sse_input_file()


    def configure_parameter_set(self, parameter_args):

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config['system']

        misc_args = input_config.simulation_config['misc']
        run_id = self.get_run_id(misc_args)
        misc_args['uuid'] = run_id

        misc_args['simulation_folder'] = self.simulation_folder

        run_folder = self.create_run_folder(misc_args)
        misc_args['run_folder'] = run_folder

        system = System(**system_args)

        self.validate_system(system)

        simulation_args = input_config.simulation_config['simulation']
        sampling_args = input_config.simulation_config['sampling']
        combined_args = simulation_args | sampling_args

        sampling_params = SamplingParameters(**combined_args)

        simulation_parameters = SSEParameters(system, sampling_params, run_folder)

        self.driver_inputs.append(simulation_parameters)

        return 0
    

    def validate_system(self, system):

        if system.model_name not in self.allowed_hamiltonians:
            raise Exception("Unrecognized Hamiltonian name for nearest-neighbor simulator. " \
            "Allowed Hamiltonians are: {}".format(self.allowed_hamiltonians))

        if system.model_name == "afm_heisenberg_fm_Z":
            if not system.geometry.bipartite:
                raise Exception("AFMHeisenbergFMZ Hamiltonian requires a bipartite lattice\n")
            
        if system.hamiltonian_parameters.J <= 0:
            raise Exception("J sets the energy scale. It must be a positive number for QMC\n")
        
        if system.interaction_range != "nearest_neighbor":
            print(system.interaction_range)
            raise Exception("NearestNeighborQMC can only be used for nearest-neighbor interactions\n")