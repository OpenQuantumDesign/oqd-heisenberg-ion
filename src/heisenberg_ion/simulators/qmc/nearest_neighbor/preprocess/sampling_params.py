class SamplingParameters:
    """
    object defining the sampling parameters for nearest-neighbor SSE
    """

    def __init__(self, **sampling_args):
        """
        extracts the sampling parameters needed for nearest neighbor QMC. Sets defaults for optional parameters if needed

        Args:
            **sampling_args (dict): key word arguments that contain the sampling parameters for nearest neighbor QMC
        """

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
        """
        checks if specified optional sampling parameters for nearest neighbor QMC is provided, if an input is not provided, a default value is returned

        Args:
            key (str): parameter key
            sampling_args (dict): user specified key value arguments for parameter set
            default_val (float): default value of SSE sampling parameter

        Returns:
            (float): value to be used for parameter specified by key
        """

        if key in sampling_args:
            var = sampling_args[key]
        else:
            var = default_val

        return var

    def update_parameters(self, parameter_dict):
        """
        updates the parameter dict associated with a parameter set with nearest neighbor SSE parameters

        Returns:
            (dict): parameter set
        """

        parameter_dict["initial_operator_list_size"] = self.initial_operator_list_size
        parameter_dict["operator_list_update_multiplier"] = self.operator_list_update_multiplier
        parameter_dict["initial_configuration_index"] = self.initial_configuration_index

        return parameter_dict


class SSEParameters:
    """
    Object to collect all SSE parameters that need to be passed to the nearest neighbor QMC engine
    """

    def __init__(self, system, sampling_parameters, run_folder):

        self.system = system
        self.sampling_parameters = sampling_parameters

        self.run_folder = run_folder
