from preprocessing.probability_table.base import ProbabilityTable
from preprocessing.system.base import System

class ProbabilityTableBuilder:

    def __init__(self, system, **sampling_args):
        
        sampling_type = sampling_args["sampling_type"]

        if sampling_args["write_tables"]:
            out_dir = sampling_args["probability_table_folder"]
            self.prob_table.write_to_files(out_dir)

        self.prob_table = ProbabilityTable.create(sampling_type, **sampling_args)

        