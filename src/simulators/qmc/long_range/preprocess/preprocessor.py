from system.base import System
from probability_table.factory import ProbabilityTableFactory
from src.common.preprocess.base import Preprocess
import uuid
import os

class LongRangeQMC(Preprocess):

    def __init__(self):

        self.num_parameter_sets = None
        self.parameter_set_list = None

        self.input_format = None
        self.bin_folder = None
        self.root_folder = None
        self.cpp_source_folder = None


    def simulate(self):

        cpp_driver = LongRangeQMC(self.root_folder, bin_dir=self.bin_folder, source_dir=self.cpp_source_folder)
        cpp_driver.call_bin()

        return 0


    def build(self, **kwargs):

        if "InputFile" in kwargs:
            self.input_format = "InputFile"
            self.build_from_input_file(self, kwargs["InputFile"])
        else:
            self.input_format = "keywordArguments"
            self.build_from_parameters(**kwargs)

        self.configure_simulation()
    

    def build_from_input_file(self, input_file_path):

        self.file_inputs = InputFileReader(input_file_path)

        self.num_parameter_sets = self.file_inputs.num_parameter_sets
        self.parameter_set_list = self.file_inputs.parameter_set_list

        self.check_single_input("RootFolder")
        self.root_folder = self.parameter_set_list[0]["RootFolder"]

        self.check_single_input("NumberOfThreads", True)

        self.extract_cli_requirements()

        self.check_unique_uuids()


    def build_from_parameters(self, **kwargs):

        self.num_parameter_sets = 1
        self.parameter_set_list = [kwargs]

        self.root_folder = self.parameter_set_list[0]["RootFolder"]
        self.extract_cli_requirements()


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
                text_line = key + "\t" + self.parameter_set_list[0][key]
                for i in range(1, self.num_parameter_sets):
                    text_line += "," + self.parameter_set_list[i][key]
                text_line += "\n"
                f.write(text_line)

        return 0
    

    def check_unique_uuids(self):

        if "Uuid" in self.file_inputs.is_param_iterable:
            uuid_set = {self.parameter_set_list[i]["Uuid"] for i in range(self.num_parameter_sets)}
            if len(uuid_set) != self.num_parameter_sets:
                raise Exception("Specified uuids are not unique\n")
            

    def extract_cli_requirements(self):

        self.bin_folder = self.extract_optional_input("BinFolder", True)
        self.cpp_source_folder = self.extract_optional_input("CppSourceFolder", True)
        if self.bin_folder == None and self.cpp_source_folder == None:
            raise Exception("No cpp binaries or source directory provided\n")
            

    def extract_optional_input(self, key, unique=False):

        if unique: 
            if self.check_single_input(key, True):
                return self.parameter_set_list[0][key]
            else:
                return None
        else:
            if self.check_input_provided(key, True):
                return self.parameter_set_list[0][key]
            else:
                return None


    def check_single_input(self, key, optional):

        input_provided = self.check_input_provided(key, optional)
        if input_provided and self.file_inputs.is_param_iterable[key]:
                raise Exception("There should only be one input for key: {}".format(key))

        return input_provided
    

    def check_input_provided(self, key, optional=False):

        input_provided = key in self.parameter_set_list[0]
        if not input_provided and not optional:
            raise Exception("Missing required key: {}".format(key))
        else:
            return input_provided


    def get_run_id(self, misc_args):

        if "Uuid" in misc_args:
            run_id = misc_args["Uuid"]
        else:
            run_id = uuid.uuid4()

        return run_id
    

    def create_run_folder(self, misc_args):

        root_folder = misc_args["RootFolder"]
        run_id = misc_args['Uuid']
        run_folder = os.path.join(root_folder, run_id)
        os.mkdir(run_folder)

        return run_folder