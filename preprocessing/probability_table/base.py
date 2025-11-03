class ProbabilityTable:

    registry = {}

    def __init__(self, **sampling_args):

        self.sampling_parameters = sampling_args

        self.spectrum_offset = None
        self.max_diag_norm = None
        self.max_over_states = None

    
    def build(self):
        pass

    def write_to_files(self, out_dir):

        self.out_dir = out_dir
        return 0

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def create(cls, name, **kwargs):

        if name not in cls.registry:
            raise ValueError(f"ProbabilityTable implementation not found for name: {name}")
        else:
            return cls.registry[name](**kwargs)
        
# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically
ProbabilityTable.register = classmethod(ProbabilityTable.register)
ProbabilityTable.create = classmethod(ProbabilityTable.create)