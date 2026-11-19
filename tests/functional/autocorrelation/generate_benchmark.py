from oqd_heisenberg_ion.common.driver.factory import DriverFactory
from oqd_heisenberg_ion.common.inputs.input_reader import InputReader
from oqd_heisenberg_ion.common.preprocessor.factory import PreprocessorFactory

input_file = "tests/functional/autocorrelation/input_files/long_range_XY.txt"
inputs = InputReader(input_file_path=input_file)
inputs.read_inputs_from_file()

inputs.read_kwarg_inputs(output_folder_name="benchmark")

parameter_set_list = inputs.parameter_set_list
preprocessor_qmc = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
driver_inputs = preprocessor_qmc.preprocess()
driver = DriverFactory.create("long_range_qmc", preprocessor_qmc.simulation_folder, driver_inputs)
driver.simulate()