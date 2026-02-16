from .utils import convert_to_snake_case


class InputFileReader:
    """
    Reads the tab-delimited input file and contains the parameters in a list of dicts. Each dict represents a parameter set
    """

    def __init__(self, input_file_path, **kwargs):
        """
        InputFileReader constructor. Extracts parameter sets from file.

        Args:
            input_file_path (str): path to input file specifying parameter sets
        """

        self.input_file = input_file_path

        self.num_parameter_sets = 1

        self.is_param_iterable = {}
        input_config = self.extract_key_value_inputs()

        self.override_file_inputs(input_config, **kwargs)

        self.extract_parameter_set_list(input_config)

        self.simulator = convert_to_snake_case(self.parameter_set_list[0]["simulator"])

    def extract_key_value_inputs(self):
        """
         Method for extracting key value pairs seperated by a tab from input file

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

                self.record_input(input_config, key, data)

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

        self.parameter_set_list = [{} for i in range(self.num_parameter_sets)]

        for key, val in input_config.items():
            if len(val) == 1:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[0]
            else:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[i]
