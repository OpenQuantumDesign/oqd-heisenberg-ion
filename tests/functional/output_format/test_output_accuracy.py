import os

import h5py
import numpy as np

from oqd_heisenberg_ion.interface.orchestrator import Orchestrator

simulator='long_range_qmc'

# System Parameters
hamiltonian_name='XY'
boundary='periodic'
interaction_range='long_range'
interaction_type='power_law'
alpha=0.0
spatial_dimension='1d'
lattice_type='rectangular'
N=5
J=1.0

# Simulation Parameters
T = 0.05
equilibration_steps=5000
simulation_steps=5000
operator_list_update_multiplier=1.73
loop_type='deterministic'

# Misc.
test_id = np.random.randint(1000, 1e6, 1)[0]
root_folder='./tests/functional/output_format/results'
output_folder_name="test_id_{}".format(test_id)
track_spin_configurations='True'
write_final_spin_configuration='True'
number_of_threads=1

# Initial Configuration Parameters
initial_configuration_index=1
initial_operator_list_size=-1

# Seeds
initial_config_seed=1
disconnected_spin_flip_seed=4
off_diagonal_update_seed=123
metropolis_insert_seed=6
metropolis_remove_seed=8
diagonal_update_seed=2
exit_leg_seed=3

orchestrator = Orchestrator(simulator=simulator, 
                            hamiltonian_name=hamiltonian_name, 
                            boundary=boundary, 
                            interaction_range=interaction_range, 
                            interaction_type=interaction_type, 
                            alpha=alpha, 
                            spatial_dimension=spatial_dimension, 
                            lattice_type=lattice_type,
                            N=N, 
                            J=J, 
                            T=T, 
                            equilibration_steps=equilibration_steps, 
                            simulation_steps=simulation_steps, 
                            operator_list_update_multiplier=operator_list_update_multiplier, 
                            loop_type=loop_type,
                            root_folder=root_folder, 
                            output_folder_name=output_folder_name, 
                            track_spin_configurations=track_spin_configurations, 
                            write_final_spin_configuration=write_final_spin_configuration,
                            number_of_threads=number_of_threads,
                            initial_configuration_index=initial_configuration_index, 
                            initial_operator_list_size=initial_operator_list_size, 
                            initial_config_seed=initial_config_seed, 
                            disconnected_spin_flip_seed=disconnected_spin_flip_seed, 
                            off_diagonal_update_seed=off_diagonal_update_seed, 
                            metropolis_insert_seed=metropolis_insert_seed, 
                            metropolis_remove_seed=metropolis_remove_seed,
                            diagonal_update_seed=diagonal_update_seed,
                            exit_leg_seed=exit_leg_seed)
orchestrator.simulate()

bechmark_shot_data_file = "/Users/shaeermoeed/Github/oqd_heisenberg_ion/tests/functional/output_format/results/benchmark_shots/spin_configurations_diagnostics.csv"
benchmark_shot_data = np.loadtxt(bechmark_shot_data_file, delimiter=",", skiprows=2, dtype=int)

num_param_sets = orchestrator.preprocessor.num_parameter_sets
run_dir = orchestrator.preprocessor.processed_configs[0]["misc"]["run_folder"]
shot_data_file = os.path.join(run_dir, 'qmc_output/spin_configurations.h5')

with h5py.File(shot_data_file, 'r') as f:
    dataset = f['shot_data']
    shot_data_array = dataset[:]

print(test_id)
print(shot_data_array)

check_arr_equal = np.all(np.equal(shot_data_array, benchmark_shot_data))

if not check_arr_equal:
    raise Exception("Shot data not equal to the benchmark shot data")
print("Shot data equal to the benchmark shot data")