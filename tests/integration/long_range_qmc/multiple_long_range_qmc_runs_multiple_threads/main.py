import os
import time

from heisenberg_ion.interface.orchestrator import Orchestrator

a = time.time()
input_file = os.path.join(
    os.getcwd(), "tests/integration/long_range_qmc/multiple_long_range_qmc_runs_multiple_threads/wrapper_inputs_2.txt"
)
print(input_file)
Orchestrator(input_file=input_file)
b = time.time()
print(b - a)
