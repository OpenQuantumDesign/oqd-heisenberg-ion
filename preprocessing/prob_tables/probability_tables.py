class ProbabilityTable:

    def __init__(self):

        self.system = None

        self.spectrum_offset = None
        self.max_diag_norm = None
        self.max_over_states = None

        return 0
    
    def build(self):
        pass

    def register(cls, name, subclass):

        cls.registry[name] = subclass

    def create(cls, name, **kwargs):

        if name not in cls.registry:
            raise ValueError(f"ProbabilityTable implementation not found for name: {name}")
        else:
            return cls.registry[name](**kwargs)
        
# Make register and create classmethods so subclasses can be added to the registry and instantiated agnostically by ProbabilityTable
ProbabilityTable.register = classmethod(ProbabilityTable.register)
ProbabilityTable.create = classmethod(ProbabilityTable.create)