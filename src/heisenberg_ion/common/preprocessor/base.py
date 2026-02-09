import os
import uuid
from datetime import datetime


class Preprocessor:
    keys_single_parameters = {"root_folder", "simulation_folder", "number_of_threads", "output_folder_name"}

    def __init__(self, parameter_set_list):

        self.parameter_set_list = parameter_set_list
        self.num_parameter_sets = len(parameter_set_list)

        self.root_folder = None
        self.simulation_folder = None

        self.driver_inputs = None

    def build(self):

        pass

    def create_output_folder(self):

        output_folder_name_provided = self.check_single_input("output_folder_name", True)

        if output_folder_name_provided:
            self.output_folder_name = self.parameter_set_list[0]["output_folder_name"]
        else:
            self.output_folder_name = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

        self.simulation_folder = os.path.join(self.root_folder, self.output_folder_name)

        os.mkdir(self.simulation_folder)

        return self.simulation_folder

    def create_run_folder(self, misc_args):

        simulation_folder = self.simulation_folder
        run_id = str(misc_args["uuid"])
        run_folder = os.path.join(simulation_folder, run_id)

        os.mkdir(run_folder)

        return run_folder

    def get_run_id(self, misc_args):

        if "uuid" in misc_args:
            run_id = misc_args["uuid"]
        else:
            run_id = uuid.uuid4()

        return run_id

    def check_input_provided(self, key, optional=False):

        input_provided = key in self.parameter_set_list[0]
        if not input_provided and not optional:
            raise Exception(f"Missing required key: {key}")
        else:
            return input_provided

    def check_unique_uuids(self):

        input_provided = self.check_input_provided("uuid", True)
        if input_provided:
            uuid_set = {self.parameter_set_list[i]["uuid"] for i in range(self.num_parameter_sets)}
            if len(uuid_set) != self.num_parameter_sets:
                raise Exception("Specified uuids are not unique\n")

    def check_single_input(self, key, optional=False):

        input_provided = self.check_input_provided(key, optional)
        if input_provided:
            val = self.parameter_set_list[0][key]
            for i in range(self.num_parameter_sets):
                if self.parameter_set_list[i][key] != val:
                    raise Exception(f"There should only be one input for key: {key}")

        return input_provided

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

    def write_input_file(self):

        simulation_folder = self.simulation_folder
        sse_input_file = os.path.join(simulation_folder, "inputs.txt")

        with open(sse_input_file, "w") as f:
            for key in self.parameter_set_list[0].keys():
                text_line = key + "\t" + str(self.parameter_set_list[0][key])
                for i in range(1, self.num_parameter_sets):
                    if key not in self.keys_single_parameters:
                        text_line += "," + str(self.parameter_set_list[i][key])
                text_line += "\n"
                f.write(text_line)

        return 0
