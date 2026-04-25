from oqd_heisenberg_ion.interface.orchestrator import Orchestrator

N_list = [5, 7, 9, 11, 13]
theta_list = [-0.02, 0.0, 0.02]
N_inputs = []
theta_inputs = []
uuid_inputs = []
for i in range(len(N_list)):
    for j in range(len(theta_list)):
        N_inputs.append(N_list[i])
        theta_inputs.append(theta_list[j])
        uuid = "hamiltonian_XY_N_{}_theta_{}_boundary_1".format(N_list[i], theta_list[j])
        uuid_inputs.append(uuid)

simulator = 'exact_diagonalization'

# System Parameters
hamiltonian_name = 'XY'
boundary = 'periodic'
interaction_range = 'nearest_neighbor'
spatial_dimension = '1d'
lattice_type = 'rectangular'
N =N_inputs
J = 1.0
theta=theta_inputs
root_folder = './Results'
output_folder_name = 'ED'
uuid=uuid_inputs

orchestrator = Orchestrator(simulator=simulator, 
                            hamiltonian_name=hamiltonian_name, 
                            boundary=boundary, 
                            interaction_range=interaction_range,  
                            spatial_dimension=spatial_dimension, 
                            lattice_type=lattice_type,
                            N=N, 
                            J=J,
                            theta=theta,
                            root_folder=root_folder, 
                            output_folder_name=output_folder_name, 
                            uuid=uuid)
orchestrator.simulate()

