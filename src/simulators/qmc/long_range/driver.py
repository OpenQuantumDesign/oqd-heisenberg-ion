import subprocess
import os
from src.common.driver.base import Driver

class LongRangeQMC(Driver):

    def __init__(self, simulation_folder, bin_dir=None, source_dir=None):

        self.input_file = os.path.join(simulation_folder, "sse_inputs.txt")
        self.build_dir = os.path.join(simulation_folder, "build")
        os.mkdir(self.build_dir)

        self.bin_dir = bin_dir
        self.source_dir = source_dir

        if self.bin_dir is None:
            self.build_from_cmake(self.source_dir)
            self.compile_source(self.source_dir)


    def build_from_cmake(self, cpp_source):

        configure_command = ["cmake", cpp_source]
        build_results = subprocess.run(configure_command, cwd=self.build_dir, check=True, capture_output=True)

        if build_results.returncode != 0:
            raise Exception("Build from cmake failed with exit code: {}.\n" \
            "Stack trace: {}".format(build_results.returncode, build_results.stderr))


    def compile_source(self):

        compile_command = ["cmake", "--build", ".", "-j"]
        compile_results = subprocess.run(compile_command, cwd=self.build_dir, check=True, capture_output=True)

        if compile_results.returncode != 0:
            raise Exception("Compile failed with exit code: {}.\n" \
            "Stack trace: {}".format(compile_results.returncode, compile_results.stderr))


    def simulate(self):

        copy_command = ["cp", "-r", self.bin_dir, self.build_dir]
        copy_results = subprocess.run(copy_command, check=True, capture_output=True)

        if copy_results.returncode != 0:
            raise Exception("Copying binaries failed with exit code: {}.\n" \
            "Stack trace: {}".format(copy_results.returncode, copy_results.stderr))
        
        run_results = subprocess.run(["./cpp_qmc"], cwd=self.build_dir, check=True, capture_output=True)

        if run_results.returncode != 0:
            raise Exception("Cpp run failed with exit code: {}.\n" \
            "Stack trace: {}".format(run_results.returncode, run_results.stderr))
        
        return 0

