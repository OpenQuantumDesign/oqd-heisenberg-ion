from .utils import convert_to_snake_case
from .utils import convert_to_pascal

input_schema = {
        "Hamiltonian Name" : {
            "DataType": "categorical",
            "Categories" : ["FMHeisenbergAFMZ", "XY", "FMHeisenbergFMZ", "XXZ", "XXZh", "AFMHeisenbergFMZ"],
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": False
            },
        "Loop Type" : {
            "DataType": "categorical", 
            "Categories": ["Determinisitic", "Heatbath", "DirectedLoop"],
            "ParameterType" : "Sampling",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "Boundary": {
            "DataType": "categorical",
            "Categories": ["Periodic", "Open"],
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            }, 
        "Interaction Range": {
            "DataType": "categorical",
            "Categories": ["NearestNeighbors", "LongRange"],
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "Interaction Type" : {
            "DataType": "categorical", 
            "Categories": ["PowerLaw", "InputMatrix"],
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "Spatial Dimension" : {
            "DataType": "categorical",
            "Categories": ["1d", "2d"],
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": False
            }, 
        "Lattice Type" : {
            "DataType" : "categorical", 
            "Categories": ["Rectangular", "Triangular"], 
            "ParameterType" : "System",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "Simulator" : {
            "DataType" : "categorical", 
            "Categories" : ["LongRangeQMC", "NearestNeighborQMC", "ExactDiagonalization"],
            "ParameterType" : "Simulation",
            "ConvertKeyCase": True,
            "ConvertValCase": True
        },

        "N" : {"DataType" : int, "ParameterType" : "System", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "Delta" : {"DataType" : float, "ParameterType" : "System", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "alpha" : {"DataType" : float, "ParameterType" : "System", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "h" : {"DataType" : float, "ParameterType" : "System", "ConvertKeyCase": False, "ConvertValCase": False},
        "J" : {"DataType" : float, "ParameterType" : "System", "ConvertKeyCase": False, "ConvertValCase": False},
        "T" : {"DataType" : float, "ParameterType" : "Sampling", "ConvertKeyCase": False, "ConvertValCase": False},

        "Equilibration Steps" : {"DataType" : int, "ParameterType" : "Simulation", 
                                "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Simulation Steps" : {"DataType" : int, "ParameterType" : "Simulation", 
                             "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Operator List Update Multiplier" : {"DataType" : float, "ParameterType" : "Simulation", 
                                          "ConvertKeyCase": True, "ConvertValCase": False}, 
        "gamma" : {"DataType" : float, "ParameterType" : "Sampling", 
                   "ConvertKeyCase": False, "ConvertValCase": False}, 
        "ksi" : {"DataType" : float, "ParameterType" : "Sampling", 
                 "ConvertKeyCase": False, "ConvertValCase": False}, 
        "Distance Dependant Offset" : {"DataType": bool, "ParameterType" : "Sampling", 
                                     "ConvertKeyCase": True, "ConvertValCase": False}, 

        "Root Folder" : {"DataType" : str, "ParameterType" : "Misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Bin Folder" : {"DataType" : str, "ParameterType" : "Misc", "ConvertKeyCase": True, "ConvertValCase": False},
        "Cpp Source Folder" : {"DataType" : str, "ParameterType" : "Misc", "ConvertKeyCase": True, "ConvertValCase": False},
        "UUID" : {"DataType" : str, "ParameterType" : "Misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Number Of Threads" : {"DataType" : int, "ParameterType" : "Misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Track Spin Configurations" : {"DataType" : bool, "ParameterType" : "Simulation", "ConvertKeyCase": True, 
                                     "ConvertValCase": False}, 
        "Write Final Spin Configuration": {"DataType" : bool, "ParameterType" : "Simulation", "ConvertKeyCase": True, 
                                        "ConvertValCase": False},

        "Initial Configuration Index" : {"DataType" : int, "ParameterType" : "Simulation", 
                                       "ConvertKeyCase": True, "ConvertValCase": False}, 
        "Initial Configuration File Path" : {"DataType" : str, "ParameterType" : "Simulation", 
                                          "ConvertKeyCase": True, "ConvertValCase": False},
        "Initial Operator List Size" : {"DataType": int, "ParameterType" : "Simulation", 
                                     "ConvertKeyCase": True, "ConvertValCase": False}
        }


class InputParser:

    def __init__(self, **config_settings):

        self.dtype_parsers = {str: self.extract_string, 
                    int: self.extract_integer, 
                    float: self.extract_float, 
                    bool: self.extract_bool,
                    "categorical": self.extract_categorical}
        
        self.input_schema = input_schema

        self.simulation_config = {}

        self.build_input_config(config_settings)


    def extract_integer(self, config_settings, key):

        return int(config_settings[key])
    
    def extract_float(self, config_settings, key):

        return float(config_settings[key])
    
    def extract_string(self, config_settings, key):

        return (config_settings[key])
    
    def extract_bool(self, config_settings, key):
        
        if type(config_settings[key]) == bool:
            return config_settings[key]
        else:
            if config_settings[key].capitalize() == "True":
                return True
            elif config_settings[key].capitalize() == "False":
                return False
            else:
                raise ValueError("Unrecognized entry for key {}, value provided: {}. "
                                "Allowed values are 'True' or 'False'\n".format(key, config_settings[key]))
    
    def extract_categorical(self, config_settings, key):

        val = config_settings[key]
        val_std = convert_to_pascal(val)

        allowed_vals = self.input_schema[key]["Categories"]
        
        if val_std in allowed_vals:
            return val
        else:
            raise ValueError("Unrecognized input for key: {}. Value provided: {}. " \
            "Allowed values: {}".format(key, val, allowed_vals))


    def build_input_config(self, config_settings):

        for key, item in config_settings.items():
            
            dtype = self.input_schema[key]["DataType"]
            param_type = self.input_schema[key]["ParameterType"]

            parser = self.dtype_parsers[dtype]
            val = parser(config_settings, key)

            if not param_type in self.simulation_config:
                self.simulation_config[param_type] = {}

            std_key = convert_to_snake_case(key) if input_schema[key]["ConvertKeyCase"] else key
            std_val = convert_to_snake_case(val) if input_schema[key]["ConvertValCase"] else val

            self.simulation_config[param_type][std_key] = std_val

        return 0