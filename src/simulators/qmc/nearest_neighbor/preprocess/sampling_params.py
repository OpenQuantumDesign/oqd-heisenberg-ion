class SamplingParameters:

    def __init__(self, sampling_args):

        self.T = sampling_args["T"]
        self.beta = 1.0/self.T

        self.equilibration_steps = sampling_args["EquilibrationSteps"]
        self.mc_steps = sampling_args["MCSteps"]

        self.M = self.set_optional_parameter("InitialOperatorListSize", sampling_args, 50)
        self.a = self.set_optional_parameter("OperatorListUpdateMultiplier", sampling_args, 1.25)
        self.init_config_index = self.set_optional_parameter("InitialConfigurationIndex", sampling_args, 0)


    def set_optional_parameter(self, key, sampling_args, default_val):

        if key in sampling_args:
            var = sampling_args[key]
        else:
            var = default_val

        return var