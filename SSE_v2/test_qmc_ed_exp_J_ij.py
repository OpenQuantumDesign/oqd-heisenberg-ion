import numpy as np
import statistical_analysis as stats
import matplotlib.pyplot as plt

N = 11
T = 0.005
beta = 1.0/T
energy_ED_list = []
energy_SSE_list = []
energy_SSE_err_list = []

alpha_list = ["exp_mu_1.1", "exp_mu_1.2", "exp_mu_1.3", "exp_mu_1.4", "exp_mu_1.5", "exp_mu_1.6", "exp_mu_1.7", "exp_mu_1.8", "exp_mu_1.9", "exp_mu_2.0", "exp_mu_1.41"]
mu_list = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
eq_drop = 0

Delta = 0.0
hamiltonian_type = int(Delta)
h = 0.0
gamma = 0.0
loop_type = "deterministic"
dist_dep_offset = 0
boundary=0
init_start_config = 1
start_config = init_start_config

for i in range(len(alpha_list)):

    alpha = alpha_list[i]
    
    if alpha[0:3] == "exp":
        J = "exp"
    else:
        J = 1.0
    
    auto_corr_drop = 1

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_{}_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, J, dist_dep_offset, boundary, T, loop_type, start_config, init_start_config)
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)

    energy_array = energy_arr_1/N
    stiffness_array = stiffness_arr_1

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)

    if alpha == alpha_list[-1]:
        energy_SSE_exp_J_ij = energy[0]
        energy_SSE_err_exp_J_ij = energy[1]
    else:
        energy_SSE_list.append(energy[0])
        energy_SSE_err_list.append(energy[1])

for i in range(len(alpha_list)):

    alpha = alpha_list[i]

    if alpha[0:3] == "exp":
        J = "exp"
    else:
        J = 1.0

    comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_0.0_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha)
    state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)

    if alpha == alpha_list[-1]:
        energy_ED_exp_J_ij = evals[0]/N
    else:
        energy_ED_list.append(evals[0]/N)

print(energy_ED_list)
print(energy_SSE_list)

print(energy_ED_exp_J_ij)
print(energy_SSE_exp_J_ij + energy_SSE_err_exp_J_ij)
print(energy_SSE_exp_J_ij - energy_SSE_err_exp_J_ij)

plt.figure()
plt.scatter(mu_list, energy_SSE_list, label='SSE', color='C0')
plt.errorbar(mu_list, energy_SSE_list, energy_SSE_err_list, fmt='None', capsize=5)
plt.plot(mu_list, energy_ED_list, label='ED', color='C3')
plt.plot(mu_list, [energy_SSE_err_exp_J_ij + energy_SSE_exp_J_ij]*len(alpha_list[:-1]), linestyle='dashed', color='C0')
plt.plot(mu_list, [-energy_SSE_err_exp_J_ij + energy_SSE_exp_J_ij]*len(alpha_list[:-1]), linestyle='dashed', color='C0')
plt.plot(mu_list, [energy_ED_exp_J_ij]*len(alpha_list[:-1]), linestyle='dashed', color='C3')
plt.ylabel('E/N')
plt.xlabel(r"$\lambda$")
#plt.plot(alpha_list, [0.0]*len(alpha_list), color='C3')
plt.legend()
plt.savefig("Exp_J_ij_Energies.png")