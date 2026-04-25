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
N=[3,125]
J=1.0

# Simulation Parameters
T = 0.05
equilibration_steps=1.0*(10**6)
simulation_steps=1.0*(10**6)
operator_list_update_multiplier=[1.5, 1.13]
loop_type='deterministic'

# Misc.
root_folder='./results'
output_folder_name='SSE_MFT_Limit_2'
track_spin_configurations='True'
write_final_spin_configuration='False'
number_of_threads=8

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

uuid_list = []
for i in range(len(N)):

    alpha_str = str(round(alpha,2))
    N_str = str(N[i])

    uuid_str = 'N_{}_hamiltonian_type_0_Delta_0.0_h_0.0_alpha_{}_gamma_0.0_ksi_0.0_J_1.0_dist_dep_offset_0_boundary_1_T_{}_{}_input_config_{}_initial_input_config_{}'.format(N_str, alpha_str, T, loop_type, initial_configuration_index, initial_configuration_index)
    uuid_list.append(uuid_str)

uuid = uuid_list

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
                            uuid=uuid, 
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