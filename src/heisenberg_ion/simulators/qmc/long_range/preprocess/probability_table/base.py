import os


class ProbabilityTable:
    """ProbabilityTable base class to generate probability table files"""

    args = {}

    def __init__(self, system, **sampling_args):
        """
        constructor builds the geometry and sets the base member variables

        Args:
            system (System): object representing the system to be simulated
            **sampling_args (dict): key word arguments specifying the sampling parameters
        """

        self.system = system
        self.system.geometry.build()

        self.sampling_parameters = sampling_args

        self.spectrum_offset = None
        self.max_diag_norm = None
        self.max_over_states = None

    def validate_system(self):
        """
        validates the system associated with the ProbabilityTable object

        Raises:
            Exception: if geometry is not long range, since probability tables are only required for long range geometries
            Exception: if the energy scale J is less than zero since that can result in negative probabilities
        """

        if not self.system.geometry.interaction_range == "long_range":
            raise Exception("Probability tables are constructed for long range systems\n")

        if self.system.hamiltonian_parameters.J <= 0:
            raise Exception("J sets the energy scale. It must be positive for long-range QMC")
        pass

    def build(self):
        """
        each ProbabilityTable subclass must implement a build method to populate tables
        """
        pass

    def write_to_files(self, out_dir):
        """
        makes the directories required for writing probability tables. Writing logic owned by subclasses

        Args:
            out_dir (str): directory path for probability tables
        """

        self.out_dir = out_dir
        self.prob_dir = os.path.join(self.out_dir, "probability_densities/")
        os.makedirs(self.prob_dir)
