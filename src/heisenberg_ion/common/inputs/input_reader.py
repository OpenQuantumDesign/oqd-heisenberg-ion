import copy
from collections.abc import Sequence

from .utils import convert_to_snake_case


class InputReader:
    """
    Extracts provided inputs and generates a list of parameter sets, with each parameter set represented as a dict.
    Inputs need to be provided either in a tab-delimited file or as key word arguments
    """

    def __init__(self, simulator=None, input_file_path=None):
        """
        initializes the InputReader

        Args:
            input_file_path (str, optional): path to input file. Defaults to None.
        """

        self.input_file = input_file_path
        if simulator is not None:
            self.simulator = convert_to_snake_case(simulator)

        self.num_parameter_sets = 1

        self.is_param_iterable = {}

        self.parameter_set_list = []

    def read_kwarg_inputs(self, **kwargs):
        """
        Extracts inputs from key word arguments and builds list of parameter sets

        Args:
            **kwargs (dict): input key word arguments provided by user
        """

        input_config = self.extract_kwarg_inputs(**kwargs)

        self.extract_parameter_set_list(input_config)

    def read_inputs_from_file(self, **overrides):
        """
        Extracts inputs from file and builds list of parameter sets

        Args:
            **overrides (dict): input key word arguments provided by user as overrides
        """

        input_config = self.extract_file_key_value_inputs()

        input_config = self.override_file_inputs(input_config, **overrides)

        self.extract_parameter_set_list(input_config)

        self.simulator = convert_to_snake_case(self.parameter_set_list[0]["simulator"])

    def extract_file_key_value_inputs(self):
        """
        Reads input file to generate a dict of key value pairs specifying parameters. Each value is a list of parameters

        Returns:
            (dict): containing the key value pairs from the tab delimited input file
        """

        input_config = {}

        with open(self.input_file) as f:
            line_count = 0

            for line in f.readlines():
                line_count += 1

                if line.startswith("#") or line.strip() == "":
                    continue

                line_data = line.strip().split("\t")

                key = line_data[0]
                data = line_data[1].strip().split(",")

                input_config = self.record_input(input_config, key, data)

        return input_config

    def extract_kwarg_inputs(self, **kwargs):
        """
        Uses key word arguments to generate a dict of key value pairs specifying parameters. Each value is a list of parameters

        Returns:
            (dict): containing the key value pairs from the tab delimited input file
        """

        input_config = {}

        for key, val in kwargs.items():
            if isinstance(val, Sequence) and not isinstance(val, str):
                input_config = self.record_input(input_config, key, val)
            else:
                input_config = self.record_input(input_config, key, [val])

        return input_config

    def record_input(self, input_config, key, data):
        """
        counts the parameter sets for a given key and records the key and the value in the input_config dictionary

        Args:
            input_config (dict): key value pairs container for input parameters
            key (str): parameter key
            data (str): parameter value

        Returns:
            (dict): containing the key value pairs from the tab delimited input file
        """

        count_entries = len(data)
        self.count_parameter_sets(count_entries, key)

        input_config[key] = data
        self.is_param_iterable[key] = count_entries != 1

        return input_config

    def override_file_inputs(self, input_config, **kwargs):
        """
        overrides file parameter specification with given parameter values

        Args:
            input_config (dict): key value pairs container for input parameters
            **kwargs (dict): key word arguments containing the overriding keys and values

        Returns:
            (dict): key value pairs container for input parameters
        """

        for key, val in kwargs.items():
            data = val.strip().split(",")
            self.record_input(input_config, key, data)

        return input_config

    def count_parameter_sets(self, count_entries, key):
        """
        updates the number of parameter sets if more than one parameter set is specified for given key

        Args:
            count_entries (int): number of entries for given key
            key (str): parameter key

        Raises:
            Exception: if number of values provided for key and number of parameter sets > 1 are inconsistent
        """

        if count_entries != 1:
            if self.num_parameter_sets == 1:
                self.num_parameter_sets = count_entries
                self.key_num_parameters_set = key

            # Can't have n>1 parameter sets in one field and m>1 parameter sets in another
            elif self.num_parameter_sets != count_entries:
                raise Exception(f"Inconistent number of entries for fields: {key} and {self.key_num_parameters_set}\n")

    def extract_parameter_set_list(self, input_config):
        """
        populates the list of parameter sets (class member) from the dictionary of input key value pairs

        Args:
            input_config (dict): contains the key value pairs specifying the parameter set
        """
        if not self.parameter_set_list:
            for i in range(self.num_parameter_sets):
                self.parameter_set_list.append({})
        elif len(self.parameter_set_list) == 1:
            for i in range(len(self.parameter_set_list), self.num_parameter_sets):
                self.parameter_set_list.append(copy.deepcopy(self.parameter_set_list[0]))
        elif len(self.parameter_set_list) != self.num_parameter_sets:
            raise Exception("Inconsistent numbers of parameter sets encountered")

        for key, val in input_config.items():
            if len(val) == 1:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[0]
            else:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[i]
