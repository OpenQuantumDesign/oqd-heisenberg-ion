import os
import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats

N=10
h = 0.0
#alpha_list=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
#alpha_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
alpha_list=[1.0, 2.0, 3.0, 4.0, 5.0, 10.0]
J = 1.0
T = 0.1
beta = J/T
energy_ED_list = []
mag_z_ED_list = []
stiffness_ED_list = []
energy_SSE_list = []
energy_SSE_err_list = []
mag_z_SSE_list = []
mag_z_SSE_err_list = []
stiffness_SSE_list = []
stiffness_SSE_err_list = []
theta = 0.001
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
boundary = 1

for i in range(len(alpha_list)):

    alpha = alpha_list[i]

    '''
    hamiltonian_type = 1
    if hamiltonian_type == -1:
        start_config = 2
    else:
        start_config = 1
    '''
    '''
    init_start_config = 0
    start_config = init_start_config
    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, loop_type, start_config, init_start_config)
    if alpha <= 3.0 and hamiltonian_type == 0:
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
    else:
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
    '''
    init_start_config = 1
    start_config = init_start_config
    file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T, loop_type, start_config, init_start_config)
    step_number_2, energy_arr_2, magnetization_arr_2, stiffness_arr_2, s_k_pi = np.loadtxt(file_2, delimiter=",", skiprows=2, unpack=True)

    '''
    init_start_config = 2
    start_config = init_start_config
    file_3 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, loop_type, start_config, init_start_config)
    step_number_3, energy_arr_3, magnetization_arr_3, stiffness_arr_3  = np.loadtxt(file_3, delimiter=",", skiprows=2, unpack=True)

    init_start_config = 3
    start_config = init_start_config
    file_4 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_boundary_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, loop_type, start_config, init_start_config)
    step_number_4, energy_arr_4, magnetization_arr_4, stiffness_arr_4 = np.loadtxt(file_4, delimiter=",", skiprows=2, unpack=True)
    '''
    drop_1 = 0
    auto_corr_drop = 2
    energy_array = stats.combine_different_runs([], energy_arr_2, 0, drop_1)
    #energy_array = stats.combine_different_runs(energy_array, energy_arr_3, 0, drop_1)
    #energy_array = stats.combine_different_runs(energy_array, energy_arr_4, 0, drop_1)
    energy_array = energy_arr_2

    stiffness_array = stats.combine_different_runs([], stiffness_arr_2, 0, drop_1)
    #stiffness_array = stats.combine_different_runs(stiffness_array, stiffness_arr_3, 0, drop_1)
    #stiffness_array = stats.combine_different_runs(stiffness_array, stiffness_arr_4, 0, drop_1)
    stiffness_array = stiffness_arr_2

    '''
    if hamiltonian_type == 0:
        energy_array = energy_arr_1
        stiffness_array = stiffness_arr_1
    elif hamiltonian_type == 1:
        energy_array = energy_arr_2
        stiffness_array = stiffness_arr_2
    elif hamiltonian_type == -1:
        energy_array = energy_arr_3
        stiffness_array = stiffness_arr_3
    '''

    comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_0.0_Heisenberg_PBC.csv".format(N, Delta, h, J, J, alpha)
    state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)
    comp_file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(N, Delta, h, J, J, alpha, theta)
    state_index_theta, evals_theta, mag_z_theta, mag_x_theta = np.loadtxt(comp_file_1, skiprows=1, delimiter=",", unpack=True)
    comp_file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(N, Delta, h, J, J, alpha, -theta)
    state_index_minus_theta, evals_minus_theta, mag_z_minus_theta, mag_x_minus_theta = np.loadtxt(comp_file_2, skiprows=1, delimiter=",", unpack=True)

    ED_magnetization = 0.0
    energy_theta = 0.0
    energy_zero = 0.0
    energy_minus_theta = 0.0
    partition_theta = 0.0
    partition_zero = 0.0
    partition_minus_theta = 0.0
    for i in range(len(evals)):
        energy_theta += np.exp(-beta * evals_theta[i]) * evals_theta[i]
        partition_theta += np.exp(-beta * evals_theta[i])
        energy_zero += np.exp(-beta * evals[i]) * evals[i]
        partition_zero += np.exp(-beta * evals[i])
        energy_minus_theta += np.exp(-beta * evals_minus_theta[i]) * evals_minus_theta[i]
        partition_minus_theta += np.exp(-beta * evals_minus_theta[i])
        ED_magnetization += mag_z[i] * np.exp(-beta * evals[i])
            
    energy_theta /= partition_theta
    energy_zero /= partition_zero
    energy_minus_theta /= partition_minus_theta
    ED_magnetization /= partition_zero
    free_energy_theta = (-1.0/beta) * np.log(partition_theta)
    free_energy_zero = (-1.0/beta) * np.log(partition_zero)
    free_energy_minus_theta = (-1.0/beta) * np.log(partition_minus_theta)
    second_derivative = (free_energy_theta/N + free_energy_minus_theta/N - 2.0*free_energy_zero/N)/(theta**2)
    #second_derivative = (evals_theta[0]/N + evals_minus_theta[0]/N - 2.0*evals[0]/N)/(theta**2)
    spin_stiffness = second_derivative
    stiffness_ED_list.append(spin_stiffness)

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_2, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

    stiffness_SSE_list.append(stiffness[0])
    stiffness_SSE_err_list.append(stiffness[1])

    energy_ED_list.append(energy_zero)
    mag_z_ED_list.append(ED_magnetization)
    energy_SSE_list.append(energy[0])
    energy_SSE_err_list.append(energy[1])
    mag_z_SSE_list.append(mag_z_exp[0])
    mag_z_SSE_err_list.append(mag_z_exp[1])

    print(step_number_2[-1])

print(stiffness_ED_list)
print(stiffness_SSE_list)
print(stiffness_SSE_err_list)

plt.figure()
plt.plot(alpha_list, [0.0]*len(alpha_list), color="C3")
plt.scatter(alpha_list, np.array(energy_ED_list)/N-np.array(energy_SSE_list)/N, label="ED-SSE", color='C0')
plt.errorbar(alpha_list, np.array(energy_ED_list)/N - np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel("E/N")
plt.title("E/N Comparison vs alpha")
plt.tight_layout()
plt.savefig("Energy per site comparison difference.png")

plt.figure()
plt.scatter(alpha_list, np.array(stiffness_SSE_list) - np.array(stiffness_ED_list), color='C0', label='SSE-ED')
plt.plot(alpha_list, [0.0]*len(alpha_list), color="C3")
plt.errorbar(alpha_list, np.array(stiffness_SSE_list) - np.array(stiffness_ED_list), np.array(stiffness_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\rho_s$")
plt.title("Superfluid Density vs alpha")
plt.tight_layout()
plt.savefig("Superfluid_Density_Comparison difference.png")

plt.figure()
plt.plot(alpha_list, np.array(energy_ED_list)/N, color="C3", label='ED')
plt.scatter(alpha_list, np.array(energy_SSE_list)/N, label="SSE", color='C0')
plt.errorbar(alpha_list, np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel("E/N")
plt.title("E/N Comparison vs alpha")
plt.tight_layout()
plt.savefig("Energy_per_site_comparison.png")

plt.figure()
plt.plot(alpha_list, np.array(mag_z_ED_list), color="C3", label='ED')
plt.scatter(alpha_list, np.array(mag_z_SSE_list), color='C0', label='SSE')
plt.errorbar(alpha_list, np.array(mag_z_SSE_list), np.array(mag_z_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel("M_z Difference")
plt.title("Magnetization Comparison vs alpha")
plt.tight_layout()
plt.savefig("Magnetization_comparison.png")

plt.figure()
plt.scatter(alpha_list, np.array(stiffness_SSE_list), color='C0', label='SSE')
plt.plot(alpha_list, np.array(stiffness_ED_list), color="C3", label='ED')
plt.errorbar(alpha_list, np.array(stiffness_SSE_list), np.array(stiffness_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\rho_s$")
plt.title(r"Superfluid Density ($\Delta$={})".format(Delta))
plt.tight_layout()
plt.savefig("Superfluid_Density_Comparison.png")