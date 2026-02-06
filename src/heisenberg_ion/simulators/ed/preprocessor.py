import os
import shutil as sh

from heisenberg_ion.common.inputs.input_parser import InputParser
from heisenberg_ion.common.preprocessor.base import Preprocessor

from ..preprocess.system.base import System


class ExactDiagonalization(Preprocessor):
    def __init__(self, parameter_set_list):

        super().__init__(parameter_set_list)

        # self.driver_inputs = []

        self.build()

    def build(self):

        self.check_single_input("root_folder")
        self.root_folder = self.parameter_set_list[0]["root_folder"]

        self.simulation_folder = self.create_output_folder()

        self.extract_cli_requirements()

        self.check_unique_uuids()

        self.configure_simulation()

    def configure_simulation(self):

        for i in range(self.num_parameter_sets):
            parameter_set = self.parameter_set_list[i]

            self.configure_parameter_set(parameter_set)

        self.write_input_file()

    def configure_parameter_set(self, parameter_args):

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config["system"]

        misc_args = input_config.simulation_config["misc"]
        run_id = self.get_run_id(misc_args)
        misc_args["uuid"] = run_id
        parameter_args["uuid"] = run_id

        misc_args["simulation_folder"] = self.simulation_folder
        parameter_args["simulation_folder"] = self.simulation_folder

        run_folder = self.create_run_folder(misc_args)
        misc_args["run_folder"] = run_folder
        parameter_args["run_folder"] = run_folder

        system = System(**system_args)
        system_args["B"] = system.hamiltonian_parameters.B
        parameter_args["B"] = system.hamiltonian_parameters.B

        if system.interaction_range == "long_range":
            J_ij_file_path = os.path.join(run_folder, "J_ij_file.csv")
            system.interactions.write_to_file(J_ij_file_path)
            misc_args["J_ij_file"] = J_ij_file_path
            parameter_args["J_ij_file"] = J_ij_file_path

        """
        theta = system_args['theta']

        ed_parameters = EDParameters(system, run_folder, theta)

        self.driver_inputs.append(ed_parameters)
        """

        return 0

    def extract_cli_requirements(self):

        julia_path = self.extract_optional_input("julia_path", True)
        if julia_path is None:
            julia_path = sh.which("julia")
            if julia_path is None:
                raise Exception("No cpp binaries or source directory provided\n")

        self.driver_inputs = {"julia_path": julia_path}


"""
class EDParameters:

    def __init__(self, system, run_folder, theta):

        self.system = system
        self.run_folder = run_folder
        self.theta = theta

        if self.system.interaction_range == "long_range":
            self.write_J_ij_file()


    def write_J_ij_file(self, system, run_folder):

        self.J_ij_file = os.path.join(run_folder, "J_ij_matrix.csv")

        if self.interaction_name == "power_law":
            system.interaction.write_to_file(self.J_ij_file)

        elif system.interaction_name == "matrix_input":
            sh.copy2(system.interaction.file_path, self.J_ij_file)
"""
