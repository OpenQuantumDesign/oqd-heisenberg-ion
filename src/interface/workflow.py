from common.inputs.input_file_reader import InputFileReader
from common.preprocessor.factory import PreprocessorFactory
from common.driver.factory import DriverFactory

class Workflow:

    def __init__(self, **kwargs):

        if "InputFile" in kwargs:
            self.build_from_file(kwargs)
        else:
            self.build_from_parameters(**kwargs)


    def build_from_file(self, input_file_path):

        file_inputs = InputFileReader(input_file_path)

        preprocessor = PreprocessorFactory.create(file_inputs.simulator, file_inputs.parameter_set_list)

        driver = DriverFactory.create(preprocessor.simulation_folder, preprocessor.driver_inputs)
        driver.simulate()

        return 0


    def build_from_parameters(self, **kwargs):

        simulator = kwargs["Simulator"]
        parameter_set_list = [kwargs]

        preprocessor = PreprocessorFactory.create(simulator, parameter_set_list)

        driver = DriverFactory.create(preprocessor.simulation_folder, preprocessor.driver_inputs)
        driver.simulate()

        return 0