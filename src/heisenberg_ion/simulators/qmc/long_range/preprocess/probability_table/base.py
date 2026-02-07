import os


class ProbabilityTable:
    args = {}

    def __init__(self, system, **sampling_args):

        self.system = system
        self.system.geometry.build()

        self.sampling_parameters = sampling_args

        self.spectrum_offset = None
        self.max_diag_norm = None
        self.max_over_states = None

    def validate_system(self):

        if not self.system.geometry.interaction_range == "long_range":
            raise Exception("Probability tables are constructed for long range systems\n")

        if self.system.hamiltonian_parameters.J <= 0:
            raise Exception("J sets the energy scale. It must be positive for QMC")
        pass

    def build(self):
        pass

    def write_to_files(self, out_dir):

        self.out_dir = out_dir
        self.prob_dir = os.path.join(self.out_dir, "probability_densities/")
        os.makedirs(self.prob_dir)

        return 0
