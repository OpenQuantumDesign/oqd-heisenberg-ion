from input_utils import standardize_string_format

input_schema = {
        "HamiltonianType" : {
            "DataType": "categorical",
            "Categories" : ["HeisenbergAFM", "XY", "HeisenbergFM", "XXZ", "XXZh"],
            "ParameterType" : "System"
            },
        "LoopType" : {
            "DataType": "categorical", 
            "Categories": ["Determinisitic", "HeatBath", "DirectedLoop"],
            "ParameterType" : "Sampling"
            },
        "BoundaryCondition": {
            "DataType": "categorical",
            "Categories": ["Periodic", "Open"],
            "ParameterType" : "System"
            }, 
        "InteractionType" : {
            "DataType": "categorical", 
            "Categories": ["PowerLaw", "InputMatrix"],
            "ParameterType" : "System"
            },
        "SpatialDimension" : {
            "DataType": "Categorical",
            "Categories": ["1d", "2d"],
            "ParameterType" : "System"
            }, 
        "LatticeType" : {
            "DataType" : "Categorical", 
            "Categories": ["Rectangular", "Triangular"], 
            "ParameterType" : "System"
            },

        "N" : {"DataType" : int, "ParameterType" : "System"}, 
        "Delta" : {"DataType" : float, "ParameterType" : "System"}, 
        "alpha" : {"DataType" : float, "ParameterType" : "System"}, 
        "h" : {"DataType" : float, "ParameterType" : "System"},
        "J" : {"DataType" : float, "ParameterType" : "System"},
        "T" : {"DataType" : float, "ParameterType" : "Sampling"},

        "EquilibrationSteps" : {"DataType" : int, "ParameterType" : "Simulation"}, 
        "SimulationSteps" : {"DataType" : int, "ParameterType" : "Simulation"}, 
        "OperatorListUpdateMultiplier" : {"DataType" : int, "ParameterType" : "Simulation"}, 
        "Gamma" : {"DataType" : float, "ParameterType" : "Sampling"}, 
        "Ksi" : {"DataType" : float, "ParameterType" : "Sampling"}, 
        "DistanceDependantOffset" : {"DataType": bool, "ParameterType" : "Sampling"}, 

        "RootFolder" : {"DataType" : str, "ParameterType" : "Misc"}, 
        "BinFolder" : {"DataType" : str, "ParameterType" : "Misc"},
        "SouceFolder" : {"DataType" : str, "ParameterType" : "Misc"},
        "UUID" : {"DataType" : str, "ParameterType" : "Misc"}, 
        "NumberOfThreads" : {"DataType" : int, "ParameterType" : "Misc"}, 
        "Track Spin Configurations" : {"DataType" : bool, "ParameterType" : "Simulation"}, 
        "WriteFinalSpinConfiguration": {"DataType" : bool, "ParameterType" : "Simulation"},

        "InitialConfigurationIndex" : {"DataType" : str, "ParameterType" : "Simulation"}, 
        "InitialConfigurationFilePath" : {"DataType" : str, "ParameterType" : "Simulation"},
        "InitialOperatorListSize" : {"DataType": int, "ParameterType" : "Simulation"}
        }

class InputParser:

    def __init__(self, schema, **config_settings):

        self.dtype_parsers = {str: self.extract_string, 
                    int: self.extract_integer, 
                    float: self.extract_float, 
                    bool: self.extract_bool,
                    "categorical": self.extract_categorical}
        
        self.input_schema = schema

        self.simulation_config = {}

        self.build_input_config(config_settings)


    def extract_integer(self, config_settings, key):

        return int(config_settings[key])
    
    def extract_float(self, config_settings, key):

        return float(config_settings[key])
    
    def extract_string(self, config_settings, key):

        return standardize_string_format(config_settings[key])
    
    def extract_bool(self, config_settings, key):

        if config_settings[key].capitalize() == "True":
            return True
        elif config_settings[key].capitalize() == "False":
            return False
        else:
            raise ValueError("Unrecognized entry for key {}, value provided: {}. "
                            "Allowed values are 'True' or 'False'\n".format(key, config_settings[key]))
    
    def extract_categorical(self, config_settings, key):

        val = config_settings[key]
        val_std = standardize_string_format(val)

        allowed_vals = self.input_schema[key]["Categories"]
        
        if val_std in allowed_vals:
            return val_std
        else:
            raise ValueError("Unrecognized input for key: {}. Value provided: {}. " \
            "Allowed values: {}".format(key, val, allowed_vals))

    def build_input_config(self, config_settings):

        for key, item in config_settings.item():
            
            dtype = self.input_schema[key]["DataType"]
            param_type = self.input_schema[key]["ParameterType"]

            parser = self.dtype_parsers[dtype]
            val = parser(config_settings, key)

            self.simulation_config[param_type][key] = val

        return 0

        
            