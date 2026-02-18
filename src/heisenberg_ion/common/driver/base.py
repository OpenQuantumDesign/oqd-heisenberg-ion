import os


class Driver:
    """Driver base class

    Calls the required simulator engine after preprocessing

    """

    def __init__(self, simulation_folder):
        """Driver base class constructor

        Args:
            simulation_folder (str): root directory for storing simulator outputs

        Raises:
                Exception: if the root simulation folder does not exist
        """

        self.simulation_folder = simulation_folder
        if not os.path.exists(self.simulation_folder):
            raise Exception("Unable to find the root simulation folder\n")

    def simulate(self):
        """Driver should always have a simulate method implemented by subclasses"""

        pass
