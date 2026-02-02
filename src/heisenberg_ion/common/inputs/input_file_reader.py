from .utils import convert_to_snake_case

class InputFileReader:

    def __init__(self, input_file_path):

        self.input_file = input_file_path

        self.num_parameter_sets = 1

        self.is_param_iterable = {}
        input_config = self.extract_key_value_inputs()

        self.extract_parameter_set_list(input_config)

        self.simulator = convert_to_snake_case(self.parameter_set_list[0]["simulator"])

    
    def extract_key_value_inputs(self):

        input_config = {}

        with open(self.input_file, 'r') as f:

            line_count = 0

            for line in f.readlines():

                line_count += 1

                if line.startswith("#") or line.strip() == "":
                    continue

                line_data = line.strip().split("\t")

                data = line_data[1].strip().split(",")
                count_entries = len(data)

                self.count_parameter_sets(count_entries, line_count)

                key = line_data[0]

                input_config[key] = data
                self.is_param_iterable[key] = (count_entries != 1)

        return input_config
    
    
    def count_parameter_sets(self, count_entries, line_number):

        if count_entries != 1:
            if self.num_parameter_sets == 1:
                self.num_parameter_sets = count_entries
                line_num_parameters_set = line_number

            elif self.num_parameter_sets != count_entries:
                raise ValueError("Inconistent number of entries in line: {} and "
                "line: {} \n".format(line_number, line_num_parameters_set))
            
        return 0
    

    def extract_parameter_set_list(self, input_config):


        self.parameter_set_list = [{} for i in range(self.num_parameter_sets)]

        for key,val in input_config.items():
            if len(val) == 1:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[0]
            else:
                for i in range(self.num_parameter_sets):
                    self.parameter_set_list[i][key] = val[i]

        return 0