from heisenberg_ion.interface.workflow import Workflow
import os

input_file = os.path.join(os.getcwd(),"tests/integration/long_range_qmc/single_long_range_qmc_run/wrapper_inputs.txt")
print(input_file)
Workflow(input_file=input_file)