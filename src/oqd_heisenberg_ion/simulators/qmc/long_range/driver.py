import os
import subprocess

from oqd_heisenberg_ion.common.driver.base import Driver


class LongRangeQMC(Driver):
    """
    Driver subclass for long range QMC
    """

    def __init__(self, simulation_folder, simulator_inputs):
        """
        builds and compiles the C++ engine if needed after extracting the build directory and input file for the QMC engine

        Args:
            simulation_folder (str): simulation output folder
            simulator_inputs (dict): specifies either the binaries or the cpp source directory
        """

        super().__init__(simulation_folder)

        self.input_file = os.path.join(os.path.abspath(simulation_folder), "inputs.txt")
        self.build_dir = os.path.join(os.path.abspath(simulation_folder), "build")
        print(self.build_dir)
        os.mkdir(self.build_dir)

        self.bin_dir = simulator_inputs["bin_folder"]
        self.source_dir = simulator_inputs["cpp_source_folder"]

        if self.bin_dir is None:
            self.build_from_cmake(self.source_dir)
            self.compile_source()

    def build_from_cmake(self, cpp_source):
        """
        method for calling cmake via Python subprocess

        Args:
            cpp_source (str): path pointing to the C++ source code
        """

        configure_command = ["cmake", cpp_source, "-DCMAKE_BUILD_TYPE=Release"]
        try:
            subprocess.run(configure_command, cwd=self.build_dir, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            print("Build from cmake failed with exit code: {}\n".format(error.returncode))
            print("STDOUT: {}".format(error.stdout))
            print("STDERR:{}".format(error.stderr))
            raise

    def compile_source(self):
        """
        method for compiling the C++ source code
        """

        compile_command = ["cmake", "--build", ".", "-j"]
        try:
            subprocess.run(compile_command, cwd=self.build_dir, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            print("Compile failed with exit code: {}\n".format(error.returncode))
            print("STDOUT: {}".format(error.stdout))
            print("STDERR:{}".format(error.stderr))
            raise

    def simulate(self):
        """
        method for calling the C++ engine for simulating the system, via a Python subprocess call

        Returns:
            (int): exit code, 0 if the simulation terminates gracefully
        """

        if self.bin_dir is not None:
            copy_command = ["cp", "-r", self.bin_dir, self.build_dir]
            try:
                subprocess.run(copy_command, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as error:
                print("Copying binaries run failed with exit code: {}\n".format(error.returncode))
                print("STDOUT: {}".format(error.stdout))
                print("STDERR:{}".format(error.stderr))
                raise

        try:
            subprocess.run(
                ["./cpp_qmc", self.input_file], cwd=self.build_dir, check=True, capture_output=True, text=True
            )
        except subprocess.CalledProcessError as error:
            print("Cpp run failed with exit code: {}\n".format(error.returncode))
            print("STDOUT: {}".format(error.stdout))
            print("STDERR:{}".format(error.stderr))
            raise

        return 0
