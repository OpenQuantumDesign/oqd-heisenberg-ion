from heisenberg_ion.common.preprocessor.base import Preprocessor
from ..preprocess.system.base import System
from .preprocess.probability_table.factory import ProbabilityTableFactory
from heisenberg_ion.common.inputs.input_parser import InputParser
import os

class LongRangeQMC(Preprocessor):

    def __init__(self, parameter_set_list):

        super().__init__(parameter_set_list)

        self.bin_folder = None
        self.cpp_source_folder = None

        self.build()
    

    def build(self):

        self.check_single_input("root_folder")
        self.root_folder = self.parameter_set_list[0]["root_folder"]

        self.simulation_folder = self.create_output_folder()

        self.check_single_input("number_of_threads", True)

        self.extract_cli_requirements()

        self.check_unique_uuids()

        self.configure_simulation()


    def configure_simulation(self):

        for i in range(self.num_parameter_sets):
            self.configure_parameter_set(self.parameter_set_list[i])

        self.write_sse_input_file()


    def configure_parameter_set(self, parameter_args):

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config['system']

        misc_args = input_config.simulation_config['misc']
        run_id = self.get_run_id(misc_args)
        misc_args['uuid'] = run_id
        parameter_args['uuid'] = run_id

        misc_args['simulation_folder'] = self.simulation_folder
        parameter_args['simulation_folder'] = self.simulation_folder

        run_folder = self.create_run_folder(misc_args)
        misc_args['run_folder'] = run_folder
        parameter_args['run_folder'] = run_folder

        system = System(**system_args)
        parameter_args['hamiltonian_type'] = system.hamiltonian_parameters.hamiltonian_type

        sampling_args = input_config.simulation_config['sampling']
        prob_table_type = sampling_args['loop_type']

        prob_table_args = ProbabilityTableFactory.extract_args(prob_table_type, **sampling_args)
        probability_table = ProbabilityTableFactory.create(prob_table_type, system, **prob_table_args)
        probability_table.write_to_files(run_folder)

        return 0
    

    def write_sse_input_file(self):

        simulation_folder = self.simulation_folder
        sse_input_file = os.path.join(simulation_folder, "sse_inputs.txt")

        with open(sse_input_file, "w") as f:
            for key in self.parameter_set_list[0].keys():
                text_line = key + "\t" + str(self.parameter_set_list[0][key])
                for i in range(1, self.num_parameter_sets):
                    if not key in self.keys_single_parameters:
                        text_line += "," + str(self.parameter_set_list[i][key])
                text_line += "\n"
                f.write(text_line)

        return 0
            

    def extract_cli_requirements(self):

        bin_folder = self.extract_optional_input("bin_folder", True)
        cpp_source_folder = self.extract_optional_input("cpp_source_folder", True)
        if bin_folder == None and cpp_source_folder == None:
            raise Exception("No cpp binaries or source directory provided\n")
        else:
            self.driver_inputs = {"bin_folder": bin_folder, "cpp_source_folder": cpp_source_folder}