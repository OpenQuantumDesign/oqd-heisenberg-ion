import os

class Driver:

    def __init__(self, simulation_folder, driver_inputs):

        self.simulation_folder = simulation_folder
        if not os.path.exists(self.simulation_folder):
            raise Exception("Unable to find the root simulation folder\n")

    def simulate(self):

        pass