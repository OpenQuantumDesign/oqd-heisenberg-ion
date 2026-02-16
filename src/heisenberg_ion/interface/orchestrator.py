from ..common.driver.factory import DriverFactory
from ..common.inputs.input_file_reader import InputFileReader
from ..common.preprocessor.factory import PreprocessorFactory


class Orchestrator:
    """
    Responsible for managing the entire workflow of the package. Sequentially reads the inputs, calls the preprocessor, driver and simulator
    """

    def __init__(self, **kwargs):
        """
        constructor for the Orchestrator class. Determines whether inputs should be read from file or are provided as key word arguments

        Args:
            **kwargs (dict): key word arguments. An input file is expected if one of the keys is 'input_file'
        """

        if "input_file" in kwargs:
            self.build_from_file(**kwargs)
        else:
            self.build_from_parameters(**kwargs)

    def build_from_file(self, input_file, **kwargs):
        """
        used if input file needs to be read for the simulation.

        Args:
            input_file (str): path to the input file

        Returns:
            (int): 0 if the program executes to completion
        """

        file_inputs = InputFileReader(input_file, **kwargs)
        simulator = file_inputs.simulator

        preprocessor = PreprocessorFactory.create(simulator, file_inputs.parameter_set_list)

        driver = DriverFactory.create(simulator, preprocessor.simulation_folder, preprocessor.driver_inputs)
        driver.simulate()

        return 0

    def build_from_parameters(self, **kwargs):
        """
        used if the simulation inputs are provided as key word arguments

        Returns:
            (int): 0 if the program executes to completion
        """

        simulator = kwargs["simulator"]
        parameter_set_list = [kwargs]

        preprocessor = PreprocessorFactory.create(simulator, parameter_set_list)

        driver = DriverFactory.create(simulator, preprocessor.simulation_folder, preprocessor.driver_inputs)
        driver.simulate()

        return 0
