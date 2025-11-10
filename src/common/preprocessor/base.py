import os
import uuid

class Preprocessor:

    def __init__(self, parameter_set_list):

        self.parameter_set_list = parameter_set_list
        self.num_parameter_sets = len(parameter_set_list)

        self.root_folder = None

        self.driver_inputs = None

    
    def build(self):

        pass


    def create_run_folder(self, misc_args):

        root_folder = misc_args['RootFolder']
        run_id = misc_args['Uuid']
        run_folder = os.path.join(root_folder, run_id)
        os.mkdir(run_folder)

        return run_folder


    def get_run_id(self, misc_args):

        if "Uuid" in misc_args:
            run_id = misc_args["Uuid"]
        else:
            run_id = uuid.uuid4()

        return run_id
    

    def check_input_provided(self, key, optional=False):

        input_provided = key in self.parameter_set_list[0]
        if not input_provided and not optional:
            raise Exception("Missing required key: {}".format(key))
        else:
            return input_provided
    

    def check_unique_uuids(self):

        input_provided = self.check_input_provided("Uuid", True)
        if input_provided:
            uuid_set = {self.parameter_set_list[i]["Uuid"] for i in range(self.num_parameter_sets)}
            if len(uuid_set) != self.num_parameter_sets:
                raise Exception("Specified uuids are not unique\n")


    def check_single_input(self, key, optional):

        input_provided = self.check_input_provided(key, optional)
        if input_provided:
            val = self.parameter_set_list[0][key]
            for i in range(self.num_parameter_sets):
                if self.parameter_set_list[i][key] != val:
                    raise Exception("There should only be one input for key: {}".format(key))

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