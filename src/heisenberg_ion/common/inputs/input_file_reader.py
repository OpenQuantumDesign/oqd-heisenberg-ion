from .utils import convert_to_snake_case


class InputFileReader:
    def __init__(self, input_file_path, **kwargs):

        self.input_file = input_file_path

        self.num_parameter_sets = 1

        self.is_param_iterable = {}
        input_config = self.extract_key_value_inputs()

        self.override_file_inputs(input_config, **kwargs)

        self.extract_parameter_set_list(input_config)

        self.simulator = convert_to_snake_case(self.parameter_set_list[0]["simulator"])

    def extract_key_value_inputs(self):

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

        count_entries = len(data)
        self.count_parameter_sets(count_entries, key)

        input_config[key] = data
        self.is_param_iterable[key] = count_entries != 1

        return input_config

    def override_file_inputs(self, input_config, **kwargs):

        for key, val in kwargs.items():
            data = val.strip().split(",")
            self.record_input(input_config, key, data)

        return input_config

    def count_parameter_sets(self, count_entries, key):

        if count_entries != 1:
            if self.num_parameter_sets == 1:
                self.num_parameter_sets = count_entries
                self.key_num_parameters_set = key

            # Can't have n>1 parameter sets in one field and m>1 parameter sets in another
            elif self.num_parameter_sets != count_entries:
                raise ValueError(f"Inconistent number of entries for fields: {key} and {self.key_num_parameters_set}\n")

        return 0

    def extract_parameter_set_list(self, input_config):

        self.parameter_set_list = [{} for i in range(self.num_parameter_sets)]

        for key, val in input_config.items():
            if len(val) == 1:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[0]
            else:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[i]

        return 0
