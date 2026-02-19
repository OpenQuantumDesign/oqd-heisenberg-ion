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

        self.initial_config_seed = self.set_optional_parameter("initial_config_seed", sampling_args, 17951893)
        self.disconnected_spin_flip_seed = self.set_optional_parameter(
            "disconnected_spin_flip_seed", sampling_args, 674604219
        )
        self.off_diagonal_update_seed = self.set_optional_parameter(
            "off_diagonal_update_seed", sampling_args, 961025794
        )
        self.metropolis_insert_seed = self.set_optional_parameter("metropolis_insert_seed", sampling_args, 148014634)
        self.metropolis_remove_seed = self.set_optional_parameter("metropolis_remove_seed", sampling_args, 148014634)
        self.diagonal_update_seed = self.set_optional_parameter("diagonal_update_seed", sampling_args, 734634223)
        self.exit_leg_seed = self.set_optional_parameter("exit_leg_seed", sampling_args, 2345346346)

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

        parameter_dict["initial_config_seed"] = self.initial_config_seed
        parameter_dict["metropolis_insert_seed"] = self.metropolis_insert_seed
        parameter_dict["metropolis_remove_seed"] = self.metropolis_remove_seed
        parameter_dict["disconnected_spin_flip_seed"] = self.disconnected_spin_flip_seed
        parameter_dict["off_diagonal_update_seed"] = self.off_diagonal_update_seed
        parameter_dict["exit_leg_seed"] = self.exit_leg_seed
        parameter_dict["diagonal_update_seed"] = self.diagonal_update_seed

        return parameter_dict


class SSEParameters:
    """
    Object to collect all SSE parameters that need to be passed to the nearest neighbor QMC engine
    """

    def __init__(self, system, sampling_parameters, run_folder):

        self.system = system
        self.sampling_parameters = sampling_parameters

        self.run_folder = run_folder
