import numpy as np
import os
import SSE_simulator as SSE_sim
import SSE_initializer as SSE_init
import exact_diagonalization as ed
import matplotlib.pyplot as plt
import datetime
import statistical_analysis as st

# Some test parameters
gamma = 0.0
ksi = 0.0
N = 8
J = 1.0
T = 0.5
beta = J/T
#Delta_list = np.linspace(0.0,0.5,num=2,endpoint=True)
Delta_list = [1.0]
#h_list = np.linspace(0.0,0.5,num=4,endpoint=True)
#h_list = np.linspace(0.0,2.0,num=5,endpoint=True)
h_list = [0.0, 1.0]
h_B_list = [h/(J*(N-1)) for h in h_list]
num_bonds = int(N * (N-1)/2)
alpha = 1.0
M_init = 10
#mc_steps=5000
#equilibration_steps=1000
mc_steps=800
equilibration_steps=500
N_l_init = 10

E_MC_List = np.zeros((len(Delta_list), len(h_list)))
E_err_List = np.zeros((len(Delta_list), len(h_list)))
E_ED_List = np.zeros((len(Delta_list), len(h_list)))

M_z_MC_List = np.zeros((len(Delta_list), len(h_list)))
M_z_err_List = np.zeros((len(Delta_list), len(h_list)))
M_z_ED_list = np.zeros((len(Delta_list), len(h_list)))

out_dir = "Results"
date_time = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = os.path.join(out_dir,"energies_table_N_{}_alpha_{}_{}.csv".format(N,alpha, date_time))

with open(out_file, 'w') as f:
    f.write("gamma={},N={},T={},J={},alpha={},mc steps={},eq steps={},Delta pts={},h pts={}\n".format(gamma, N, T, J, alpha, mc_steps, equilibration_steps, len(Delta_list), len(h_list)))
    f.write("Delta, h, E0 ED, E0 SSE, E0 SSE Error, M_z ED, M_z SSE, M_z SSE Error\n")
    for i1 in range(len(Delta_list)):
        for i2 in range(len(h_list)):

            Delta = np.round(Delta_list[i1],1)
            h = h_list[i2]
            h_B = h_B_list[i2]
            
            if N == 2:
                H_mat = ed.construct_Hamiltonian_2(Delta, h=h)
                evals, evecs = np.linalg.eigh(H_mat)
            elif N == 3:
                H_mat = ed.construct_Hamiltonian_3(Delta, alpha, J=J, h=h)
                evals, evecs = np.linalg.eigh(H_mat)
            
            comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha)
            state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)

            #diag_update_prob_table, max_norm_diag_probs, diag_norm = SSE_init.generate_diag_update_table(N, num_bonds, epsilon, Delta, alpha, h_B)

            #vertex_weights, spectrum_offset = SSE_init.generate_vertex_weights(epsilon, h_B, Delta, alpha, N, num_bonds)

            #loop_update_prob_table = SSE_init.generate_heat_bath_prob_table(num_bonds, vertex_weights)
            sites = SSE_init.geometry(N, num_bonds)

            diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, loop_update_table = SSE_init.compute_prob_tables_directed_loops(num_bonds, sites, alpha, gamma, h_B, Delta, ksi)
            #diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, loop_update_table = SSE_init.compute_prob_tables_heat_bath(num_bonds, sites, alpha, gamma, h_B, Delta)
            if h <= 0.3:
                starting_spin_config_1 = np.zeros(N, dtype=int)
                starting_spin_config_2 = np.zeros(N, dtype=int)
                starting_spin_config_3 = np.zeros(N, dtype=int)
                starting_spin_config_4 = np.zeros(N, dtype=int)
                for i in range(N):
                    starting_spin_config_1[i] = 1
                    starting_spin_config_2[i] = -1
                    starting_spin_config_3[i] = (-1)**i
                    starting_spin_config_4[i] = np.random.choice([1,-1])
            else:
                starting_spin_config_1 = -np.ones(N)
                starting_spin_config_2 = -np.ones(N)
                starting_spin_config_3 = -np.ones(N)
                starting_spin_config_4 = -np.ones(N)
            energy_array_1, magnetization_array_1 = SSE_sim.simulate_XXZ(N, M_init, num_bonds, equilibration_steps, mc_steps, sites, beta, 
                diag_prob_table, max_over_states, max_diag_norm, loop_update_table, starting_spin_config_1)
            
            np.savetxt("Results/step_energy_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_1_timestamp_{}".format(N,alpha,Delta,h,T, date_time), energy_array_1, header="E")
            np.savetxt("Results/step_magnetization_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_1_timestamp_{}".format(N,alpha,Delta,h,T, date_time), magnetization_array_1, header="M_z")
            
            energy_array_2, magnetization_array_2 = SSE_sim.simulate_XXZ(N, M_init, num_bonds, equilibration_steps, mc_steps, sites, beta, 
                diag_prob_table, max_over_states, max_diag_norm, loop_update_table, starting_spin_config_2)
            
            np.savetxt("Results/step_energy_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_2_timestamp_{}".format(N,alpha,Delta,h,T, date_time), energy_array_2, header="E")
            np.savetxt("Results/step_magnetization_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_2_timestamp_{}".format(N,alpha,Delta,h,T, date_time), magnetization_array_2, header="M_z")
            
            energy_array_3, magnetization_array_3 = SSE_sim.simulate_XXZ(N, M_init, num_bonds, equilibration_steps, mc_steps, sites, beta, 
                diag_prob_table, max_over_states, max_diag_norm, loop_update_table, starting_spin_config_3)
            
            np.savetxt("Results/step_energy_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_3_timestamp_{}".format(N,alpha,Delta,h,T, date_time), energy_array_3, header="E")
            np.savetxt("Results/step_magnetization_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_3_timestamp_{}".format(N,alpha,Delta,h,T, date_time), magnetization_array_3, header="M_z")
            
            energy_array_4, magnetization_array_4 = SSE_sim.simulate_XXZ(N, M_init, num_bonds, equilibration_steps, mc_steps, sites, beta, 
                diag_prob_table, max_over_states, max_diag_norm, loop_update_table, starting_spin_config_4)
            
            np.savetxt("Results/step_energy_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_4_timestamp_{}".format(N,alpha,Delta,h,T, date_time), energy_array_4, header="E")
            np.savetxt("Results/step_magnetization_N_{}_alpha_{}_Delta_{}_h_{}_T_{}_start_4_timestamp_{}".format(N,alpha,Delta,h,T, date_time), magnetization_array_4, header="M_z")

            energy_array = np.zeros(4*mc_steps)
            energy_array[0:mc_steps] = energy_array_1
            energy_array[mc_steps:2*mc_steps] = energy_array_2
            energy_array[2*mc_steps:3*mc_steps] = energy_array_3
            energy_array[3*mc_steps:4*mc_steps] = energy_array_4

            magnetization_array = np.zeros(4*mc_steps)
            magnetization_array[0:mc_steps] = magnetization_array_1
            magnetization_array[mc_steps:2*mc_steps] = magnetization_array_2
            magnetization_array[2*mc_steps:3*mc_steps] = magnetization_array_3
            magnetization_array[3*mc_steps:4*mc_steps] = magnetization_array_4

            energy_mean, energy_error = st.statistics_binning(energy_array)
            magnetization_mean, magnetization_error = st.statistics_binning(magnetization_array)

            energy_array = energy_array + J*spectrum_offset

            ED_energy = 0.0
            ED_magnetization = 0.0
            partition_func = 0.0
            for i in range(len(evals)):
                ED_energy += evals[i]*np.exp(-beta * evals[i])
                ED_magnetization +=  mag_z[i] * np.exp(-beta * evals[i])
                partition_func += np.exp(-beta * evals[i])
            ED_energy /= partition_func
            ED_magnetization /= partition_func

            E_ED_List[i1,i2] = ED_energy
            E_MC_List[i1,i2] = energy_mean + J*spectrum_offset
            E_err_List[i1,i2] = energy_error

            M_z_ED_list[i1,i2] = ED_magnetization
            M_z_MC_List[i1,i2] = magnetization_mean
            M_z_err_List[i1,i2] = magnetization_error

            print("h=",h)
            print("offset=",spectrum_offset)

            f.write(str(Delta) + "," + str(h) + "," + str(ED_energy) + "," + str(energy_mean + J*spectrum_offset) + "," + str(energy_error) + "," + str(ED_magnetization) + "," + str(magnetization_mean) + "," + str(magnetization_error) + "\n")


with open(out_file, 'r') as f:
    line = f.readline()
    line_data = line.split(",")
    for label in line_data:
        if label.strip().split("=")[0] == "Delta pts":
            num_delta_pts = int(label.strip().split("=")[1].strip())
        if label.strip().split("=")[0] == "h pts":
            num_h_pts = int(label.strip().split("=")[1].strip())

Delta_arr = np.zeros(num_delta_pts)
h_arr = np.zeros(num_h_pts)

E_ED_arr = np.zeros((num_delta_pts, num_h_pts))
E_MC_arr = np.zeros((num_delta_pts, num_h_pts))
E_err_arr = np.zeros((num_delta_pts, num_h_pts))

M_ED_arr = np.zeros((num_delta_pts, num_h_pts))
M_MC_arr = np.zeros((num_delta_pts, num_h_pts))
M_err_arr = np.zeros((num_delta_pts, num_h_pts))

data_from_file = np.loadtxt(out_file, skiprows=2, delimiter=",")

for i in range(num_delta_pts):
    Delta_arr[i] = data_from_file[i*num_h_pts,0]
    for j in range(num_h_pts):
        h_arr[j] = data_from_file[i*num_h_pts + j,1]
        E_ED_arr[i,j] = data_from_file[i*num_h_pts + j, 2]
        E_MC_arr[i,j] = data_from_file[i*num_h_pts + j, 3]
        E_err_arr[i,j] = data_from_file[i*num_h_pts + j, 4]
        M_ED_arr[i,j] = data_from_file[i*num_h_pts + j, 5]
        M_MC_arr[i,j] = data_from_file[i*num_h_pts + j, 6]
        M_err_arr[i,j] = data_from_file[i*num_h_pts + j, 7]


plt.figure()
for i in range(len(Delta_arr)):
    plt.plot(h_arr, E_ED_arr[i,:], label=r'$\Delta$={}'.format(Delta_arr[i]), color='C{}'.format(i))
    plt.scatter(h_arr, E_MC_arr[i,:], marker='o', color='C{}'.format(i))
    plt.errorbar(h_arr, E_MC_arr[i,:], E_err_arr[i,:], fmt="none", capsize=5, color='C{}'.format(i))

plt.xlabel(r"$h$")
plt.ylabel(r"$E$")
plt.title(r"N={}, $\alpha$={}".format(N,alpha))
plt.legend()
plt.savefig("Results/Energy_N_{}_alpha_{}.png".format(N,alpha))

plt.figure()
for i in range(len(Delta_arr)):
    plt.plot(h_arr, M_ED_arr[i,:], label=r'$\Delta$={}'.format(Delta_arr[i]), color='C{}'.format(i))
    plt.scatter(h_arr, M_MC_arr[i,:], marker='o', color='C{}'.format(i))
    plt.errorbar(h_arr, M_MC_arr[i,:], M_err_arr[i,:], fmt="none", capsize=5, color='C{}'.format(i))

plt.xlabel(r"$h$")
plt.ylabel(r"$M_z$")
plt.title(r"N={}, $\alpha$={}".format(N,alpha))
plt.legend()
plt.savefig("Results/Magnetization_N_{}_alpha_{}.png".format(N,alpha))

'''
plt.figure()
plt.plot(Delta_arr, E_ED_arr[:,0], label='ED', color = 'C3')
plt.scatter(Delta_arr, E_MC_arr[:,0], marker='o', label='SSE', color='C0')
plt.errorbar(Delta_arr, E_MC_arr[:,0], E_err_arr[:,0], fmt="none", capsize=5, color='C0')

plt.xlabel(r"$\Delta$")
plt.ylabel(r"$E$")
plt.title(r"N={}, $\alpha$={}".format(N,alpha))
plt.legend()
plt.savefig("Results/Energy_N_{}_alpha_{}.png".format(N,alpha))

plt.figure()
plt.plot(Delta_arr, M_ED_arr[:,0], label='ED', color = 'C3')
plt.scatter(Delta_arr, M_MC_arr[:,0], marker='o', label='SSE', color='C0')
plt.errorbar(Delta_arr, M_MC_arr[:,0], M_err_arr[:,0], fmt="none", capsize=5, color='C0')

plt.xlabel(r"$\Delta$")
plt.ylabel(r"$M_z$")
plt.title(r"N={}, $\alpha$={}".format(N,alpha))
plt.legend()
plt.savefig("Results/Magnetization_N_{}_alpha_{}.png".format(N,alpha))
'''

