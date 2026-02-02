from .utils import convert_to_snake_case
from .utils import convert_to_pascal

input_schema = {
        "hamiltonian_name" : {
            "DataType": "categorical",
            "Categories" : ["fm_heisenberg_afm_Z", "XY", "fm_heisenberg_fm_Z", "XXZ", "XXZh", "afm_heisenberg_fm_Z"],
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": False
            },
        "loop_type" : {
            "DataType": "categorical", 
            "Categories": ["deterministic", "heatbath", "directed_loop"],
            "ParameterType" : "sampling",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "boundary": {
            "DataType": "categorical",
            "Categories": ["periodic", "open"],
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            }, 
        "interaction_range": {
            "DataType": "categorical",
            "Categories": ["nearest_neighbors", "long_range"],
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "interaction_type" : {
            "DataType": "categorical", 
            "Categories": ["power_law", "input_matrix"],
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "spatial_dimension" : {
            "DataType": "categorical",
            "Categories": ["1d", "2d"],
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": False
            }, 
        "lattice_type" : {
            "DataType" : "categorical", 
            "Categories": ["rectangular", "triangular"], 
            "ParameterType" : "system",
            "ConvertKeyCase": True,
            "ConvertValCase": True
            },
        "simulator" : {
            "DataType" : "categorical", 
            "Categories" : ["long_range_qmc", "nearest_neighbor_qmc", "exact_diagonalization"],
            "ParameterType" : "simulation",
            "ConvertKeyCase": True,
            "ConvertValCase": True
        },

        "N" : {"DataType" : int, "ParameterType" : "system", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "Delta" : {"DataType" : float, "ParameterType" : "system", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "alpha" : {"DataType" : float, "ParameterType" : "system", "ConvertKeyCase": False, "ConvertValCase": False}, 
        "h" : {"DataType" : float, "ParameterType" : "system", "ConvertKeyCase": False, "ConvertValCase": False},
        "J" : {"DataType" : float, "ParameterType" : "system", "ConvertKeyCase": False, "ConvertValCase": False},
        "T" : {"DataType" : float, "ParameterType" : "sampling", "ConvertKeyCase": False, "ConvertValCase": False},

        "equilibration_steps" : {"DataType" : int, "ParameterType" : "simulation", 
                                "ConvertKeyCase": True, "ConvertValCase": False}, 
        "simulation_steps" : {"DataType" : int, "ParameterType" : "simulation", 
                             "ConvertKeyCase": True, "ConvertValCase": False}, 
        "operator_list_update_multiplier" : {"DataType" : float, "ParameterType" : "simulation", 
                                          "ConvertKeyCase": True, "ConvertValCase": False}, 
        "gamma" : {"DataType" : float, "ParameterType" : "sampling", 
                   "ConvertKeyCase": False, "ConvertValCase": False}, 
        "ksi" : {"DataType" : float, "ParameterType" : "sampling", 
                 "ConvertKeyCase": False, "ConvertValCase": False}, 
        "distance_dependent_offset" : {"DataType": bool, "ParameterType" : "sampling", 
                                     "ConvertKeyCase": True, "ConvertValCase": False}, 

        "root_folder" : {"DataType" : str, "ParameterType" : "misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "bin_folder" : {"DataType" : str, "ParameterType" : "misc", "ConvertKeyCase": True, "ConvertValCase": False},
        "cpp_source_folder" : {"DataType" : str, "ParameterType" : "misc", "ConvertKeyCase": True, "ConvertValCase": False},
        "uuid" : {"DataType" : str, "ParameterType" : "misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "number_of_threads" : {"DataType" : int, "ParameterType" : "misc", "ConvertKeyCase": True, "ConvertValCase": False}, 
        "track_spin_configurations" : {"DataType" : bool, "ParameterType" : "simulation", "ConvertKeyCase": True, 
                                     "ConvertValCase": False}, 
        "write_final_spin_configuration": {"DataType" : bool, "ParameterType" : "simulation", "ConvertKeyCase": True, 
                                        "ConvertValCase": False},

        "initial_configuration_index" : {"DataType" : int, "ParameterType" : "simulation", 
                                       "ConvertKeyCase": True, "ConvertValCase": False}, 
        "initial_configuration_file_path" : {"DataType" : str, "ParameterType" : "simulation", 
                                          "ConvertKeyCase": True, "ConvertValCase": False},
        "initial_operator_list_size" : {"DataType": int, "ParameterType" : "simulation", 
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

        allowed_vals = self.input_schema[key]["Categories"]
        
        if val in allowed_vals:
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

            std_key = convert_to_snake_case(key) if self.input_schema[key]["ConvertKeyCase"] else key
            std_val = convert_to_snake_case(val) if self.input_schema[key]["ConvertValCase"] else val

            self.simulation_config[param_type][std_key] = std_val

        return 0