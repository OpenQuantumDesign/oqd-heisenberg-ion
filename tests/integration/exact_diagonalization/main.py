import os

from heisenberg_ion.interface.workflow import Workflow

input_file = os.path.join(os.getcwd(), "tests/integration/exact_diagonalization/wrapper_inputs.txt")
print(input_file)
Workflow(input_file=input_file)
