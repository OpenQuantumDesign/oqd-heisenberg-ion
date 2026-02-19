import os

from oqd_heisenberg_ion.common.inputs.input_parser import InputParser
from oqd_heisenberg_ion.common.preprocessor.base import Preprocessor

from ...preprocess.system.base import System
from .preprocess.probability_table.factory import ProbabilityTableFactory


class LongRangeQMC(Preprocessor):
    """
    Preprocessor subclass for long range QMC. Configures the parameter sets and writes the input file for the engine.
    Also calls the probability table builder to write the tables to files for use in the QMC simulation
    """

    def __init__(self, parameter_set_list):
        """
        executes the preprocessing logic for long range QMC

        Args:
            parameter_set_list (list[dict]): list of parameter sets specified as dicts with unparsed values
        """

        super().__init__(parameter_set_list)

        self.bin_folder = None
        self.cpp_source_folder = None

    def preprocess(self):
        """
        validates user inputs that need to be unique and creates the simulation output folder. Extracts the inputs for the long range QMC driver
        and executes the parameter set configuration logic
        """

        self.check_single_input("root_folder")

        self.simulation_folder = self.create_output_folder()

        self.check_single_input("number_of_threads", True)

        self.extract_cli_requirements()

        self.check_unique_uuids()

        self.configure_simulation()

        return self.driver_inputs

    def configure_simulation(self):
        """
        configures the inputs corresponding to each parameter set requested and writes the input file for the QMC engine
        """

        for i in range(self.num_parameter_sets):
            self.configure_parameter_set(self.parameter_set_list[i])

        self.write_input_file()

    def configure_parameter_set(self, parameter_args):
        """
        implements the configuration logic for a single parameter set. Parses inputs, creates the parameter set output directory and
        calls the probability table builder for each parameter set

        Args:
            parameter_args (dict): contains a single unparsed parameter set
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
        system_args = system.update_parameters(system_args)

        sampling_args = input_config.simulation_config["sampling"]
        prob_table_type = sampling_args["loop_type"]

        prob_table_args = ProbabilityTableFactory.extract_args(prob_table_type, **sampling_args)
        probability_table = ProbabilityTableFactory.create(prob_table_type, system, **prob_table_args)
        probability_table.write_to_files(run_folder)

        self.processed_configs.append(input_config.simulation_config)

    def extract_cli_requirements(self):
        """
        prepares the driver inputs by extracting the binary folder and cpp source folder if they are provided. If not provided, attempts to find the source folder in the expected location
        """

        bin_folder = self.extract_optional_input("bin_folder", True)
        cpp_source_folder = self.extract_optional_input("cpp_source_folder", True)

        if bin_folder is None and cpp_source_folder is None:
            cpp_source_folder = os.path.dirname(os.path.abspath(__file__)) + "/engine/"

        self.driver_inputs = {"bin_folder": bin_folder, "cpp_source_folder": cpp_source_folder}
