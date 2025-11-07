import numpy as np
import os
import deterministic_NN_SSE_simulator as SSE_sim
import deterministic_NN_SSE_initializer as SSE_init
import legacy.ed.NN_exact_diagonalization as ed_iso
import matplotlib.pyplot as plt

# Some test parameters
N = 8
J = 1.0
T_list = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
beta_list = [J/T for T in T_list]
boundary_conditions = "OBC"
#interaction_type = "Anti-ferromagnetic"
interaction_type = "Anti-ferromagnetic"

M_init = 10
#mc_steps=5000
#equilibration_steps=1000
mc_steps=5000
equilibration_steps=1000

E_MC_List = np.zeros(len(beta_list))
E_err_List = np.zeros(len(beta_list))
E_ED_List = np.zeros(len(beta_list))

M_z_MC_List = np.zeros(len(beta_list))
M_z_err_List = np.zeros(len(beta_list))
M_z_ED_list = np.zeros(len(beta_list))

out_dir = "Results"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = os.path.join(out_dir,"estimators_table_N_{}_{}_isotropic_{}.csv".format(N, interaction_type, boundary_conditions))

with open(out_file, 'w') as f:
    f.write("N={},J={},mc steps={},eq steps={}, num_beta_pts={}, interaction_type={}, boundary_conditions = {}\n".format(N, J, mc_steps, equilibration_steps, len(beta_list), interaction_type, boundary_conditions))
    f.write("T, beta, E0 ED, E0 SSE, E0 SSE Error, M_z ED, M_z SSE, M_z SSE Error\n")
    for i1 in range(len(T_list)):
        
        boundary_index = SSE_init.boundary_map[boundary_conditions]
        interaction_type_index = SSE_init.isotropic_interaction_type_map[interaction_type]

        '''
        if N == 2:
            H_mat_2 = ed_iso.construct_Hamiltonian_2(1.0, interaction_type_index, J=1.0)
            evals, evecs = np.linalg.eigh(H_mat_2)
        elif N == 3:
            H_mat_2 = ed_iso.construct_Hamiltonian_3(1.0, boundary_conditions, interaction_type_index, J=1.0)
            evals, evecs = np.linalg.eigh(H_mat_2)
        '''
        comparison_file = "eigenvalues_magnetization_{}_Heisenberg_{}.csv".format(N, boundary_conditions)
        state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)

        temperature = T_list[i1]
        beta = beta_list[i1]

        sites = SSE_init.geometry(N, boundary_conditions, interaction_type)

        energy_array, energy_mean, energy_error, spectrum_offset, magnetization_array, magnetization_mean, magnetization_error, step_spin = SSE_sim.simulate_isotropic(N, M_init, equilibration_steps, mc_steps, sites, beta, interaction_type, boundary_conditions)

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

        E_ED_List[i1] = ED_energy
        E_MC_List[i1] = energy_mean + J*spectrum_offset
        E_err_List[i1] = energy_error

        M_z_ED_list[i1] = ED_magnetization
        M_z_MC_List[i1] = magnetization_mean
        M_z_err_List[i1] = magnetization_error

        np.savetxt(os.path.join(out_dir, "Isotropic_N_{}_T_{}_J_{}_mc_steps_{}_eq_steps_{}_interaction_type_{}_boundary_conditions_{}.csv".format(N, temperature, J, mc_steps, equilibration_steps, interaction_type, boundary_conditions)), step_spin, delimiter=",", fmt="%d")

        f.write(str(temperature) + "," + str(beta) + "," + str(ED_energy) + "," + str(energy_mean + J*spectrum_offset) + "," + str(energy_error) + "," + str(ED_magnetization) + "," + str(magnetization_mean) + "," + str(magnetization_error) + "\n")


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
M_z_ED_arr = np.zeros(num_beta_pts)
M_z_MC_arr= np.zeros(num_beta_pts)
M_z_err_arr = np.zeros(num_beta_pts)
T_arr, beta_arr, E_ED_arr, E_MC_arr, E_err_arr, M_z_ED_arr, M_z_MC_arr, M_z_err_arr = np.loadtxt(out_file, skiprows=2, delimiter=",", unpack=True)

plt.figure()
plt.plot(beta_arr, E_ED_arr, label='ED', color='C3')
plt.scatter(beta_arr, E_MC_arr, marker='o', label='SSE', color='C0')
plt.errorbar(beta_arr, E_MC_arr, E_err_arr, fmt="none", capsize=5, color='C0')

plt.xlabel(r"$\beta$")
plt.ylabel(r"$E$")
plt.title(r"N={}, $\Delta$=1, $h$=0, {}, {}".format(N, interaction_type, boundary_conditions))
plt.legend()
plt.tight_layout()
plt.savefig("Results/Energy_N_{}_{}_isotropic_{}.png".format(N, interaction_type, boundary_conditions))

plt.figure()
plt.plot(beta_arr, M_z_ED_arr, label='ED', color='C3')
plt.scatter(beta_arr, M_z_MC_arr, marker='o', label='SSE', color='C0')
plt.errorbar(beta_arr, M_z_MC_arr, M_z_err_arr, fmt="none", capsize=5, color='C0')

plt.xlabel(r"$\beta$")
plt.ylabel(r"$M_z$")
plt.title(r"N={}, $\Delta$=1, $h$=0, {}, {}".format(N, interaction_type, boundary_conditions))
plt.legend()
plt.tight_layout()
plt.savefig("Results/Magnetization_N_{}_{}_isotropic_{}.png".format(N, interaction_type, boundary_conditions))