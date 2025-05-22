import numpy as np
import os
import deterministic_NN_SSE_simulator as SSE_sim
import deterministic_NN_SSE_initializer as SSE_init
import NN_exact_diagonalization as ed_iso
import matplotlib.pyplot as plt

# Some test parameters
N = 3
J = 1.0
T_list = [0.1, 0.2, 0.3, 0.4, 0.5]
beta_list = [J/T for T in T_list]
print(beta_list)
boundary_conditions = "PBC"
interaction_type = "Ferromagnetic"
M_init = 10
#mc_steps=5000
#equilibration_steps=1000
mc_steps=100000
equilibration_steps=10000

E_MC_List = np.zeros(len(beta_list))
E_err_List = np.zeros(len(beta_list))
E_ED_List = np.zeros(len(beta_list))

out_dir = "Results"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = os.path.join(out_dir,"energies_table_N_{}_XY_{}_{}.csv".format(N, interaction_type, boundary_conditions))

with open(out_file, 'w') as f:
    f.write("N={},J={},mc steps={},eq steps={}, num_beta_pts={}, boundary_conditions = {}\n".format(N, J, mc_steps, equilibration_steps, len(beta_list), boundary_conditions))
    f.write("T, beta, E0 ED, E0 SSE, E0 SSE Error\n")
    for i1 in range(len(T_list)):

        temperature = T_list[i1]
        beta = beta_list[i1]

        sites = SSE_init.geometry(N, boundary_conditions, interaction_type)
        interaction_type_index = SSE_init.isotropic_interaction_type_map[interaction_type]

        energy_array, energy_mean, energy_error, spectrum_offset = SSE_sim.simulate_XY(N, M_init, equilibration_steps, mc_steps, sites, beta, boundary_conditions)

        energy_array = energy_array + J*spectrum_offset
        boundary_index = SSE_init.boundary_map[boundary_conditions]
        
        if N == 2:
            H_mat_2 = ed_iso.construct_Hamiltonian_2(0.0, interaction_type_index, J=1.0)
        elif N == 3:
            H_mat_2 = ed_iso.construct_Hamiltonian_3(0.0, boundary_conditions, interaction_type_index, J=1.0)

        print(H_mat_2)
        evals, evecs = np.linalg.eigh(H_mat_2)

        ED_energy = 0.0
        partition_func = 0.0
        for i in range(len(evals)):
            ED_energy += evals[i]*np.exp(-beta * evals[i])
            partition_func += np.exp(-beta * evals[i])
        ED_energy /= partition_func

        E_ED_List[i1] = ED_energy
        E_MC_List[i1] = energy_mean + J*spectrum_offset
        E_err_List[i1] = energy_error

        f.write(str(temperature) + "," + str(beta) + "," + str(ED_energy) + "," + str(energy_mean + J*spectrum_offset) + "," + str(energy_error) + "\n")

with open(out_file, 'r') as f:
    line = f.readline()
    line_data = line.split(",")
    for label in line_data:
        if label.strip().split("=")[0] == "num_beta_pts":
            num_beta_pts = int(label.strip().split("=")[1].strip())
        if label.strip().split("=")[0] == "J":
            J = float(label.strip().split("=")[1].strip())

beta_arr = np.zeros(num_beta_pts)

E_ED_arr = np.zeros(num_beta_pts)
E_MC_arr = np.zeros(num_beta_pts)
E_err_arr = np.zeros(num_beta_pts)
data_from_file = np.loadtxt(out_file, skiprows=2, delimiter=",")

for i in range(num_beta_pts):
    beta_arr[i] = data_from_file[i,1]
    E_ED_arr[i] = data_from_file[i, 2]
    E_MC_arr[i] = data_from_file[i, 3]
    E_err_arr[i] = data_from_file[i, 4]

plt.figure()
plt.plot(beta_arr, E_ED_arr, label='ED')
plt.scatter(beta_arr, E_MC_arr, marker='o', label='SSE')
plt.errorbar(beta_arr, E_MC_arr, E_err_arr, fmt="none", capsize=5)

plt.xlabel(r"$\beta$")
plt.ylabel(r"$E$")
plt.title(r"N={}, $\Delta$=1, $h$=0, {}, {}".format(N, interaction_type, boundary_conditions))
plt.legend()
plt.savefig("Results/N_{}_{}_XY_{}.png".format(N, interaction_type, boundary_conditions))
