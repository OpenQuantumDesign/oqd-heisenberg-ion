from heisenberg_ion.common.inputs.input_parser import InputParser
from heisenberg_ion.common.preprocessor.base import Preprocessor

from ...preprocess.system.base import System
from .preprocess.sampling_params import SamplingParameters, SSEParameters


class NearestNeighborQMC(Preprocessor):
    """
    Preprocessor subclass for nearest neighbor QMC. Configures the parameter sets and writes the input file containing the preprocessor updated parameter sets.
    """

    allowed_hamiltonians = {"afm_heisenberg_fm_Z", "XY", "fm_heisenberg_fm_Z"}

    def __init__(self, parameter_set_list):
        """
        executes the preprocessing logic for nearest neighbor qmc

        Args:
            parameter_set_list (list[dict]): list of unparsed parameters sets
        """

        super().__init__(parameter_set_list)

        self.driver_inputs = []

    def preprocess(self):
        """
        configures the input parameter sets for the nearest neighbor qmc calculation,
        checks that a single root folder is specified, creates the simulation output folder
        and validates the uuids.
        """

        self.check_single_input("root_folder")
        self.root_folder = self.parameter_set_list[0]["root_folder"]

        self.simulation_folder = self.create_output_folder()

        self.check_unique_uuids()

        self.configure_simulation()

        return self.driver_inputs

    def configure_simulation(self):
        """
        configures the input parameter sets and writes a tab delimited file containing the configured parameters
        """

        for i in range(self.num_parameter_sets):
            parameter_set = self.parameter_set_list[i]

            self.configure_parameter_set(parameter_set)

        self.write_input_file()

    def configure_parameter_set(self, parameter_args):
        """
        configures a single parameter set, specifies and creates the parameter set output folder, defines and validates the system and the sampling parameters for nearest neighbor qmc

        Args:
            parameter_args (dict): unparsed parameter set
        """

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config["system"]

        misc_args = input_config.simulation_config["misc"]
        run_id = self.get_run_id(misc_args)
        misc_args["uuid"] = run_id
        misc_args["output_folder_name"] = self.output_folder_name

        misc_args["simulation_folder"] = self.simulation_folder

        run_folder = self.create_run_folder(misc_args)
        misc_args["run_folder"] = run_folder

        system = System(**system_args)
        system.geometry.build()

        self.validate_system(system)

        system_args = system.update_parameters(system_args)

        simulation_args = input_config.simulation_config["simulation"]
        sampling_args = input_config.simulation_config["sampling"]
        combined_args = simulation_args | sampling_args

        sampling_params = SamplingParameters(**combined_args)
        sampling_args = sampling_params.update_parameters(sampling_args)

        simulation_parameters = SSEParameters(system, sampling_params, run_folder)

        self.driver_inputs.append(simulation_parameters)

        self.processed_configs.append(input_config.simulation_config)

    def validate_system(self, system):
        """
        validates the system for nearest neighbor qmc

        Args:
            system (System): object defining the system to be simulated

        Raises:
            Exception: if the specified Hamiltonian is not allowed for the nearest neighbor qmc simulator
            Exception: if the hamiltonian is anti-ferromagnetic and the geometry is not bipartite
            Exception: if J <= 0 and the hamiltonian is not anti-ferromagnetic
            Exception: if the interaction range is not nearest neighbor
        """

        if system.model_name not in self.allowed_hamiltonians:
            raise Exception(
                "Unrecognized Hamiltonian name for nearest-neighbor simulator. "
                f"Allowed Hamiltonians are: {self.allowed_hamiltonians}"
            )

        if system.model_name == "afm_heisenberg_fm_Z":
            if not system.geometry.bipartite:
                raise Exception("afm_heisenberg_fm_Z Hamiltonian requires a bipartite lattice\n")

        if system.hamiltonian_parameters.J <= 0 and system.model_name != "afm_heisenberg_fm_Z":
            raise Exception(
                "J sets the energy scale. It must be a positive number for QMC unless the model is afm_heisenberg_fm_Z on a bipartite lattice\n"
            )

        if system.interaction_range != "nearest_neighbor":
            print(system.interaction_range)
            raise Exception("NearestNeighborQMC can only be used for nearest-neighbor interactions\n")
