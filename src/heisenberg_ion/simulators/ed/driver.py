import os
import subprocess

from heisenberg_ion.common.driver.base import Driver


class ExactDiagonalization(Driver):
    def __init__(self, simulation_folder, simulator_inputs):

        super().__init__(simulation_folder)

        self.input_file = os.path.join(simulation_folder, "inputs.txt")
        self.julia_path = simulator_inputs["julia_path"]
        self.ed_engine = os.path.dirname(os.path.abspath(__file__)) + "/engine.jl"

    def simulate(self):

        try:
            subprocess.run(
                [self.julia_path, self.ed_engine, self.input_file], check=True, capture_output=True, text=True
            )
        except subprocess.CalledProcessError as error:
            print(f"Julia execution failed with exit code: {error.returncode}\n")
            print(f"STDOUT: {error.stdout}")
            print(f"STDERR:{error.stderr}")
            raise

        return 0
