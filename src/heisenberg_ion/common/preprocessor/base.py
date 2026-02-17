import os
import uuid
from datetime import datetime


class Preprocessor:
    """
    Preprocessor base class. Handles the preprocessing step before a driver is called
    Carries a set of keys for which a single value should always be specified.
    These are: "root_folder", "simulation_folder", "number_of_threads", "output_folder_name"

    Raises:
        Exception: if a required parameter is not specified
        Exception: if a uuid is repeated
        Exception: if multiple values are provided for a parameter that should have a single value
    """

    keys_single_parameters = {"root_folder", "simulation_folder", "number_of_threads", "output_folder_name"}

    def __init__(self, parameter_set_list):
        """
        constructor for the Preprocessor base class

        Args:
            parameter_set_list (list[dict]): list of unprocessed parameter sets specified as dicts
        """

        self.parameter_set_list = parameter_set_list
        self.num_parameter_sets = len(parameter_set_list)

        self.root_folder = None
        self.simulation_folder = None

        self.driver_inputs = None

    def build(self):
        """
        Each subclass must implement a build method
        """

        pass

    def create_output_folder(self):
        """
        creates the output folder to contain the simulation outputs inside the root folder

        Returns:
            (str): class attribute simulation_folder, the path to the simulation output folder
        """

        output_folder_name_provided = self.check_single_input("output_folder_name", True)

        if output_folder_name_provided:
            self.output_folder_name = self.parameter_set_list[0]["output_folder_name"]
        else:
            self.output_folder_name = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

        self.simulation_folder = os.path.abspath(os.path.join(self.root_folder, self.output_folder_name))

        os.mkdir(self.simulation_folder)

        return self.simulation_folder

    def create_run_folder(self, misc_args):
        """
        creates the run folder for a given parameter set inside the simulation output folder

        Args:
            misc_args (dict): contains the key "uuid" to specify the run_id corresponding to the associated parameter set

        Returns:
            (str): path to the parameter set run folder
        """

        simulation_folder = self.simulation_folder
        run_id = str(misc_args["uuid"])
        run_folder = os.path.join(simulation_folder, run_id)

        os.mkdir(run_folder)

        return run_folder

    def get_run_id(self, misc_args):
        """
        Extracts the uuid if specified, and creates one if not specified in the arguments

        Args:
            misc_args (dict): contains the key "uuid" to specify the run_id corresponding to the associated parameter set

        Returns:
            (str): the uuid associated with the parameter set
        """

        if "uuid" in misc_args:
            run_id = misc_args["uuid"]
        else:
            run_id = uuid.uuid4()

        return run_id

    def check_input_provided(self, key, optional=False):
        """
        checks whether an input has been specified in the parameter set

        Args:
            key (str): _description_
            optional (bool, optional): determines whether value associated with key is optional. Defaults to False.

        Raises:
            Exception: is a required input is not provided

        Returns:
            (bool): determines whether an input was provided for the given key
        """

        input_provided = key in self.parameter_set_list[0]
        if not input_provided and not optional:
            raise Exception(f"Missing required key: {key}")
        else:
            return input_provided

    def check_unique_uuids(self):
        """
        checks if multiple uuids are provided and if they are all unique. Does not throw an error if no uuid is provided

        Raises:
            Exception: if a uuid is repeated
        """

        input_provided = self.check_input_provided("uuid", True)
        if input_provided:
            uuid_set = {self.parameter_set_list[i]["uuid"] for i in range(self.num_parameter_sets)}
            if len(uuid_set) != self.num_parameter_sets:
                raise Exception("Specified uuids are not unique\n")

    def check_single_input(self, key, optional=False):
        """
        checks if multiple inputs are provided for a given key if it should have a single input

        Args:
            key (str): the parameter key
            optional (bool, optional): specifies if the parameter is optional. Defaults to False.

        Raises:
            Exception: if multiple inputs provided for a key that should have a single input

        Returns:
            (bool): determines whether the specified input was provided
        """

        input_provided = self.check_input_provided(key, optional)
        if input_provided:
            val = self.parameter_set_list[0][key]
            for i in range(self.num_parameter_sets):
                if self.parameter_set_list[i][key] != val:
                    raise Exception(f"There should only be one input for key: {key}")

        return input_provided

    def extract_optional_input(self, key, unique=False):
        """
        extracts an optional input if provided and return None otherwise

        Args:
            key (str): parameter key
            unique (bool, optional): determines if the value needs to be unique. Defaults to False.

        Returns:
            (str): the value associated with key if provided, otherwise returns None
        """

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
        """
        writes a tab-delimited input file for the engine using the parameter set values
        """

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
