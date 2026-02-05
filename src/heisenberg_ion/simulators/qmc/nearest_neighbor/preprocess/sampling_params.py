class SamplingParameters:

    def __init__(self, **sampling_args):

        self.T = sampling_args["T"]
        self.beta = 1.0/self.T

        self.equilibration_steps = sampling_args["equilibration_steps"]
        self.mc_steps = sampling_args["simulation_steps"]

        self.M = self.set_optional_parameter("initial_operator_list_size", sampling_args, 50)
        self.a = self.set_optional_parameter("operator_list_update_multiplier", sampling_args, 1.25)
        self.init_config_index = self.set_optional_parameter("initial_configuration_index", sampling_args, 0)


    def set_optional_parameter(self, key, sampling_args, default_val):

        if key in sampling_args:
            var = sampling_args[key]
        else:
            var = default_val

        return var