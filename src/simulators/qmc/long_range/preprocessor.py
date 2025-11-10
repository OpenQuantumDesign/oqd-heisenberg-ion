from src.common.preprocessor.base import Preprocessor
from preprocess.system.base import System
from preprocess.probability_table.factory import ProbabilityTableFactory
from common.inputs.input_parser import InputParser
import os

class LongRangeQMC(Preprocessor):

    def __init__(self, parameter_set_list):

        super.__init__(parameter_set_list)

        self.bin_folder = None
        self.cpp_source_folder = None

        self.build()
    

    def build(self):

        self.check_single_input("RootFolder")
        self.root_folder = self.parameter_set_list[0]["RootFolder"]

        self.check_single_input("NumberOfThreads", True)

        self.extract_cli_requirements()

        self.check_unique_uuids()

        self.configure_simulation()


    def configure_simulation(self):

        for i in range(self.num_parameter_set):

            parameter_set = self.parameter_set_list[i]

            self.configure_parameter_set(parameter_set)

            self.write_see_input_file()


    def configure_parameter_set(self, parameter_args):

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config['System']

        misc_args = input_config.simulation_config['Misc']
        run_id = self.get_run_id(misc_args)
        misc_args['Uuid'] = run_id

        run_folder = self.create_run_folder(misc_args)
        misc_args['RunFolder'] = run_folder

        system = System(system_args)

        sampling_args = input_config.simulation_config['Sampling']
        prob_table_type = sampling_args['LoopType']

        probability_table = ProbabilityTableFactory.create(prob_table_type, system, sampling_args)
        probability_table.write(run_folder)

        return 0
    

    def write_see_input_file(self):

        root_folder_path = self.parameter_set_list[0]["RootFolder"]
        sse_input_file = os.path.join(root_folder_path, "sse_inputs.txt")

        with open(sse_input_file, "w") as f:
            for key in self.parameter_set_list[0].keys():
                text_line = key + "\t" + str(self.parameter_set_list[0][key])
                for i in range(1, self.num_parameter_sets):
                    text_line += "," + str(self.parameter_set_list[i][key])
                text_line += "\n"
                f.write(text_line)

        return 0
            

    def extract_cli_requirements(self):

        bin_folder = self.extract_optional_input("BinFolder", True)
        cpp_source_folder = self.extract_optional_input("CppSourceFolder", True)
        if bin_folder == None and cpp_source_folder == None:
            raise Exception("No cpp binaries or source directory provided\n")
        else:
            self.driver_inputs = {"BinFolder": bin_folder, "SourceFolder": cpp_source_folder}