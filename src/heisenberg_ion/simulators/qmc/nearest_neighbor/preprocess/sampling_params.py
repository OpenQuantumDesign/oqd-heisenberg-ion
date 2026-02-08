class SamplingParameters:
    def __init__(self, **sampling_args):

        self.T = sampling_args["T"]
        self.beta = 1.0 / self.T

        self.equilibration_steps = sampling_args["equilibration_steps"]
        self.mc_steps = sampling_args["simulation_steps"]

        self.initial_operator_list_size = self.set_optional_parameter("initial_operator_list_size", sampling_args, 50)
        self.operator_list_update_multiplier = self.set_optional_parameter(
            "operator_list_update_multiplier", sampling_args, 1.25
        )
        self.initial_configuration_index = self.set_optional_parameter("initial_configuration_index", sampling_args, 0)

    def set_optional_parameter(self, key, sampling_args, default_val):

        if key in sampling_args:
            var = sampling_args[key]
        else:
            var = default_val

        return var

    def update_parameters(self, parameter_dict):

        parameter_dict["initial_operator_list_size"] = self.initial_operator_list_size
        parameter_dict["operator_list_update_multiplier"] = self.operator_list_update_multiplier
        parameter_dict["initial_configuration_index"] = self.initial_configuration_index

        return parameter_dict


class SSEParameters:
    def __init__(self, system, sampling_parameters, run_folder):

        self.system = system
        self.sampling_parameters = sampling_parameters

        self.run_folder = run_folder
