from input_utils import standardize_string_format

input_schema = {

        "HamiltonianType" : {
            "dtype": "categorical",
            "Categories" : ["HeisenbergAFM", "XY", "HeisenbergFM", "XXZ", "XXZh"]
            },
        "LoopType" : {
            "dype": "categorical", 
            "Categories": ["Determinisitic", "HeatBath", "DirectedLoop"]
            },
        "BoundaryCondition": {
            "dtype": "categorical",
            "Categories": ["Periodic", "Open"]
            }, 
        "InteractionType" : {
            "dtype": "categorical", 
            "Categories": ["PowerLaw", "InputMatrix"]
            },
        "SpatialDimension" : {
            "dtype": "Categorical",
            "Categories": ["1", "2"]
            }, 
        "LatticeType" : {
            "dtype" : "Categorical", 
            "Categories": ["Rectangular", "Triangular"]
            },

        "N" : {"dtype" : int}, 
        "Delta" : {"dtype" : float}, 
        "alpha" : {"dtype" : float}, 
        "h" : {"dtype" : float}, 
        "J" : {"dtype" : float},
        "T" : {"dtype" : float},

        "EquilibrationSteps" : {"dtype" : int}, 
        "SimulationSteps" : {"dtype" : int}, 
        "OperatorListUpdateMultiplier" : {"dtype" : int}, 
        "Gamma" : {"dtype" : float}, 
        "Ksi" : {"dtype" : float}, 
        "DistanceDependantOffset" : {"dtype": bool}, 

        "RootFolder" : {"dtype" : str}, 
        "UUID" : {"dtype" : str}, 
        "NumberOfThreads" : {"dtype" : int}, 
        "Track Spin Configurations" : {"dtype" : bool}, 
        "WriteFinalSpinConfiguration": {"dtype" : bool}, 

        "InitialConfigurationIndex" : {"dtype" : str}, 
        "InitialConfigurationFilePath" : {"dype" : str},
        "InitialOperatorListSize" : {"dtype": int}

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
            
            dtype = self.input_schema[key]["dtype"]
            parser = self.dtype_parsers[dtype] 
            val = parser(config_settings, key)

            self.simulation_config[key] = val

        return 0

        
            