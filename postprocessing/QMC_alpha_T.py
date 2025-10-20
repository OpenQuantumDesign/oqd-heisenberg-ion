import os
import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats
import spin_wave_analysis as lsw

N=51
h = 0.0
alpha_list = []
alpha_list_2 = [2.0 + i*0.2 for i in range(6)]
alpha_list += alpha_list_2
print(alpha_list)
J = 1.0
T_list = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
beta_list = []

stiffness_SSE_1_list = []
stiffness_SSE_1_err_list = []

start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0

find_alpha_val = 2.8
stiffness_find_alpha = []
stiffness_err_find_alpha = []

fig, ax1 = plt.subplots()
for T_1 in T_list:

    beta_1 = J/T_1

    for i in range(len(alpha_list)):

        alpha = round(alpha_list[i], 1)
        
        norm = 0.0
        for r in range(1, int((N+1)/2)):
            norm += 1.0/(r**(alpha-2))

        folder = "SSE"

        init_start_config = 1
        start_config = init_start_config
        boundary=1
        file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(folder, N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T_1, loop_type, start_config, init_start_config)
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
        
        drop_1 = int(step_number_1[-1]/2)
        auto_corr_drop = 2

        print(step_number_1[-1])

        if hamiltonian_type == 0:
            stiffness_array = np.array(stiffness_arr_1)/norm

        mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
        stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

        stiffness_SSE_1_list.append(stiffness[0])
        stiffness_SSE_1_err_list.append(stiffness[1])

        if alpha == find_alpha_val:
            stiffness_find_alpha.append(stiffness[0])
            stiffness_err_find_alpha.append(stiffness[1])

    beta_1 = round(beta_1, 3)
    beta_list.append(beta_1)
    ax1.scatter(alpha_list, np.array(stiffness_SSE_1_list), label=r'$\beta={}$'.format((beta_1)))
    ax1.errorbar(alpha_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt='None', capsize=5)

    stiffness_SSE_1_list = []
    stiffness_SSE_1_err_list = []

#ax1.plot(alpha_list, [1/np.pi] * len(alpha_list), color='green', linestyle='dashed')
ax1.legend(prop={'size': 12})

ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$\rho_s$", fontsize=16)

plt.savefig("rho_s_T.png", bbox_inches='tight')

plt.figure()
plt.scatter(T_list, stiffness_find_alpha, label=r'$\alpha$={}'.format(find_alpha_val))
plt.errorbar(T_list, stiffness_find_alpha, stiffness_err_find_alpha, fmt='None', capsize=5)
plt.xlabel(r'$T$')
plt.ylabel(r'$\rho_s$')
plt.savefig("rho_alpha_{}.png".format(find_alpha_val))