from heisenberg_ion.interface.workflow import Workflow
import os

input_file = os.path.join(os.getcwd(),"tests/integration/wrapper_inputs.txt")
print(input_file)
Workflow(InputFile=input_file)