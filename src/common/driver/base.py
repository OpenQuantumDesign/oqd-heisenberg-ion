import os

class Driver:

    def __init__(self, root_folder):

        self.root_folder = root_folder
        if not os.path.exists(self.root_folder):
            raise Exception("Unable to find the root simulation folder\n")

    def simulate(self):

        pass