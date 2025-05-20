import numpy as np
import os
import SSE_simulator as SSE_sim
import SSE_initializer as SSE_init
import exact_diagonalization as ed
import matplotlib.pyplot as plt

# Some test parameters
gamma = 0.1
N = 3
J = 1.0
T = 0.1
beta = J/T
Delta_list = np.linspace(-2.0,3.0,num=5,endpoint=True)
#Delta_list = np.linspace(0.1,2.0,num=3,endpoint=True)
h_list = np.linspace(0.0,4.0,num=5,endpoint=True)
#h_list = np.linspace(0.0,2.0,num=5,endpoint=True)
h_B_list = h_list/(J*(N-1))
#Delta_list = [0.8]
#h_list = [0.15, 0.25, 0.3]
#h_B_list = h_list
num_bonds = int(N * (N-1)/2)
alpha = 10.0
M_init = 10
#mc_steps=5000
#equilibration_steps=1000
mc_steps=10000
equilibration_steps=1000
N_l_init = 10

E_MC_List = np.zeros((len(Delta_list), len(h_list)))
E_err_List = np.zeros((len(Delta_list), len(h_list)))
E_ED_List = np.zeros((len(Delta_list), len(h_list)))

out_dir = "Results"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = os.path.join(out_dir,"energies_table_N_{}_alpha_{}.csv".format(N,alpha))

with open(out_file, 'w') as f:
    f.write("gamma={},N={},T={},J={},alpha={},mc steps={},eq steps={},Delta pts={},h pts={}\n".format(gamma, N, T, J, alpha, mc_steps, equilibration_steps, len(Delta_list), len(h_list)))
    f.write("Delta, h, E0 ED, E0 SSE, E0 SSE Error\n")
    for i1 in range(len(Delta_list)):
        for i2 in range(len(h_list)):

            Delta = Delta_list[i1]
            h = h_list[i2]
            h_B = h_B_list[i2]

            #diag_update_prob_table, max_norm_diag_probs, diag_norm = SSE_init.generate_diag_update_table(N, num_bonds, epsilon, Delta, alpha, h_B)

            #vertex_weights, spectrum_offset = SSE_init.generate_vertex_weights(epsilon, h_B, Delta, alpha, N, num_bonds)

            #loop_update_prob_table = SSE_init.generate_heat_bath_prob_table(num_bonds, vertex_weights)
            sites = SSE_init.geometry(N, num_bonds)

            diag_prob_table, max_over_states, max_diag_norm, vertex_weights, spectrum_offset, directed_loop_prob_table = SSE_init.compute_prob_tables_directed_loops(num_bonds, sites, alpha, gamma, h_B, Delta)
            #loop_update_prob_table = SSE_init.directed_loop_prob_table(num_bonds, vertex_weights, sites, alpha, epsilon, h_B, Delta)

            energy_array, energy_mean, energy_error = SSE_sim.simulate_XXZ(N, M_init, num_bonds, equilibration_steps, mc_steps, sites, beta, 
                diag_prob_table, max_over_states, max_diag_norm, directed_loop_prob_table)

            energy_array = energy_array + J*spectrum_offset

            if N == 2:
                H_mat = ed.construct_Hamiltonian_2(Delta, h=h)
            elif N == 3:
                H_mat = ed.construct_Hamiltonian_3(Delta, alpha, J=J, h=h)

            evals, evecs = np.linalg.eigh(H_mat)

            ED_energy = 0.0
            partition_func = 0.0
            for i in range(len(evals)):
                ED_energy += evals[i]*np.exp(-beta * evals[i])
                partition_func += np.exp(-beta * evals[i])
            ED_energy /= partition_func

            E_ED_List[i1,i2] = ED_energy
            E_MC_List[i1,i2] = energy_mean + J*spectrum_offset
            E_err_List[i1,i2] = energy_error

            f.write(str(Delta) + "," + str(h) + "," + str(ED_energy) + "," + str(energy_mean + J*spectrum_offset) + "," + str(energy_error) + "\n")


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
data_from_file = np.loadtxt(out_file, skiprows=2, delimiter=",")

for i in range(num_delta_pts):
    Delta_arr[i] = data_from_file[i*num_h_pts,0]
    for j in range(num_h_pts):
        h_arr[j] = data_from_file[i*num_h_pts + j,1]
        E_ED_arr[i,j] = data_from_file[i*num_h_pts + j, 2]
        E_MC_arr[i,j] = data_from_file[i*num_h_pts + j, 3]
        E_err_arr[i,j] = data_from_file[i*num_h_pts + j, 4]

plt.figure()
for i in range(len(Delta_list)):
    plt.plot(h_arr, E_ED_arr[i,:], label=r'$\Delta$={}'.format(Delta_arr[i]))
    plt.scatter(h_arr, E_MC_arr[i,:], marker='o')
    plt.errorbar(h_arr, E_MC_arr[i,:], E_err_arr[i,:], fmt="none", capsize=5)

plt.xlabel(r"$h$")
plt.ylabel(r"$E$")
plt.title(r"N={}, $\alpha$={}".format(N,alpha))
plt.legend()
plt.savefig("Results/N_{}_alpha_{}.png".format(N,alpha))