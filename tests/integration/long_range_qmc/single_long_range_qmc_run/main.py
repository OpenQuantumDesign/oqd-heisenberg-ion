import os

from heisenberg_ion.interface.orchestrator import Orchestrator

input_file = os.path.join(os.getcwd(), "tests/integration/long_range_qmc/single_long_range_qmc_run/wrapper_inputs.txt")
print(input_file)
Orchestrator(input_file=input_file)
