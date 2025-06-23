import os
import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats

def combine_different_runs(data_1, data_2, drop_samples_1, drop_samples_2):

    data_3 = []
    data_1 = data_1[drop_samples_1:]
    data_2 = data_2[drop_samples_2:]

    for i in range(len(data_1)):
        data_3.append(data_1[i])

    for i in range(len(data_2)):
        data_3.append(data_2[i])

    return data_3

N=9
Delta_list = [-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
Delta_list = [-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-14.0,-13.0,-12.0, -11.0]
Delta_list = [-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-14.0,-13.0,-12.0, -11.0,-10.0,-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0]
#Delta_list = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
#Delta_list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
h = 0.0
alpha=1.0
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
theta = 0.1
start_config = 1
auto_corr_drop = 20
eq_drop=10000
loop_type = "directed_loops"
dist_dep_offset = 0

#file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}/Cluster Histogram.csv".format(N, -20.0, h, alpha, 0.1, dist_dep_offset, loop_type, start_config)
#sites, cluster_probs = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
#plt.bar(sites, cluster_probs)
#plt.show()

for Delta in Delta_list:

    # Best run for N=10
    '''
    if Delta <= -15.0:
        gamma = 0.0
        eq_drop=1800000
        auto_corr_drop = 10
        start_config = 1
    elif Delta <= -9.0:
        gamma = 0.0
        eq_drop=1000000
        auto_corr_drop = 10
        start_config = 1
    else:
        gamma = 0.0
        eq_drop=1000000
        auto_corr_drop = 10
        start_config = 1
    '''

    gamma = 0.1
    if Delta <= -11.0:
        #start_config = 1
        eq_drop = 0
        gamma = 0.0
        dist_dep_offset = 1
    elif Delta <0.0:
        start_config = 2
        gamma = 0.1
        eq_drop = 10000
        dist_dep_offset = 0
    else:
        start_config=1
        gamma = 0.1
        eq_drop = 10000
        dist_dep_offset = 0
    auto_corr_drop = 1

    gamma=0.0
    eq_drop=0
    dist_dep_offset=1

    init_start_config = 0
    start_config = init_start_config if Delta == -1 else -2
    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)

    init_start_config = 1
    start_config = init_start_config if Delta == -1 else -2
    file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_2, energy_arr_2, magnetization_arr_2, stiffness_arr_2 = np.loadtxt(file_2, delimiter=",", skiprows=2, unpack=True)

    init_start_config = 2
    start_config = init_start_config if Delta == -1 else -2
    file_3 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_3, energy_arr_3, magnetization_arr_3, stiffness_arr_3 = np.loadtxt(file_3, delimiter=",", skiprows=2, unpack=True)

    init_start_config = 3
    start_config = init_start_config if Delta == -1 else -2
    file_4 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_4, energy_arr_4, magnetization_arr_4, stiffness_arr_4 = np.loadtxt(file_4, delimiter=",", skiprows=2, unpack=True)

    init_start_config = 4
    start_config = init_start_config if Delta == -1 else -2
    file_5 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_5, energy_arr_5, magnetization_arr_5, stiffness_arr_5 = np.loadtxt(file_5, delimiter=",", skiprows=2, unpack=True)

    drop_1 = 50000
    energy_array = combine_different_runs(energy_arr_1, energy_arr_2, drop_1, drop_1)
    energy_array = combine_different_runs(energy_array, energy_arr_3, 0, drop_1)
    energy_array = combine_different_runs(energy_array, energy_arr_4, 0, drop_1)
    #energy_array = combine_different_runs(energy_array, energy_arr_5, 0, drop_1)

    stiffness_array = combine_different_runs(stiffness_arr_1, stiffness_arr_2, drop_1, drop_1)
    stiffness_array = combine_different_runs(stiffness_array, stiffness_arr_3, 0, drop_1)
    stiffness_array = combine_different_runs(stiffness_array, stiffness_arr_4, 0, drop_1)
    #stiffness_array = combine_different_runs(stiffness_array, stiffness_arr_5, 0, drop_1)

    energy_arr_1 = energy_array
    stiffness_arr_1 = stiffness_array

    print(step_number_1[-1])

    #file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_directed_loops_input_config_1/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma)
    #step_number_2, energy_arr_2, magnetization_arr_2 = np.loadtxt(file_2, delimiter=",", skiprows=2, unpack=True)

    comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_0.0_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha)
    state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)
    comp_file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha, theta)
    state_index_theta, evals_theta, mag_z_theta, mag_x_theta = np.loadtxt(comp_file_1, skiprows=1, delimiter=",", unpack=True)
    comp_file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha, -theta)
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
    second_derivative = (energy_theta + energy_minus_theta - 2.0*energy_zero)/(theta**2)
    #second_derivative = (evals_theta[0] - evals[0])/(theta**2)
    spin_stiffness = (3.0/(2.0 * N))*second_derivative
    stiffness_ED_list.append(spin_stiffness)

    energy = stats.statistics_binning(energy_arr_1, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_arr_1, auto_corr_drop, eq_drop)

    stiffness_SSE_list.append(stiffness[0])
    stiffness_SSE_err_list.append(stiffness[1])

    #energy = combine_dataset_stats(energy_arr_1, energy_arr_2)
    #mag_z_exp = combine_dataset_stats(magnetization_arr_1, magnetization_arr_2)
    '''
    print(energy)
    print(mag_z_exp)
    print(energy_zero)
    print(ED_magnetization)
    print(stiffness)
    '''

    energy_ED_list.append(energy_zero)
    mag_z_ED_list.append(ED_magnetization)
    energy_SSE_list.append(energy[0])
    energy_SSE_err_list.append(energy[1])
    mag_z_SSE_list.append(mag_z_exp[0])
    mag_z_SSE_err_list.append(mag_z_exp[1])

#print(stiffness_SSE_list)
#print(stiffness_ED_list)

plt.figure()
plt.plot(Delta_list, [0.0]*len(Delta_list), color="C3")
plt.scatter(Delta_list, np.array(energy_ED_list)/N-np.array(energy_SSE_list)/N, label="ED-SSE", color='C0')
plt.errorbar(Delta_list, np.array(energy_ED_list)/N - np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("E/N Difference")
plt.title("E/N Comparison")
plt.tight_layout()
plt.savefig("Energy per site comparison difference.png")

plt.figure()
plt.plot(Delta_list, [0.0]*len(Delta_list), color="C3")
plt.scatter(Delta_list, np.array(mag_z_ED_list)-np.array(mag_z_SSE_list), color='C0', label='ED-SSE')
plt.errorbar(Delta_list, np.array(mag_z_ED_list)-np.array(mag_z_SSE_list), np.array(mag_z_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("M_z Difference")
plt.title("Magnetization Comparison")
plt.tight_layout()
plt.savefig("Magnetization_comparison difference.png")

plt.figure()
plt.scatter(Delta_list, np.array(stiffness_SSE_list) - np.array(stiffness_ED_list), color='C0', label='SSE-ED')
plt.plot(Delta_list, [0.0]*len(Delta_list), color="C3")
plt.errorbar(Delta_list, np.array(stiffness_SSE_list) - np.array(stiffness_ED_list), np.array(stiffness_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel(r"$\rho_s$")
plt.title("Superfluid Density")
plt.tight_layout()
plt.savefig("Superfluid_Density_Comparison difference.png")

plt.figure()
plt.plot(Delta_list, np.array(energy_ED_list)/N, color="C3", label='ED')
plt.scatter(Delta_list, np.array(energy_SSE_list)/N, label="SSE", color='C0')
plt.errorbar(Delta_list, np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("E/N")
plt.title("E/N Comparison")
plt.tight_layout()
plt.savefig("Energy per site comparison.png")

plt.figure()
plt.plot(Delta_list, np.array(mag_z_ED_list), color="C3", label='ED')
plt.scatter(Delta_list, np.array(mag_z_SSE_list), color='C0', label='SSE')
plt.errorbar(Delta_list, np.array(mag_z_SSE_list), np.array(mag_z_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("M_z Difference")
plt.title("Magnetization Comparison")
plt.tight_layout()
plt.savefig("Magnetization_comparison.png")

plt.figure()
plt.scatter(Delta_list, np.array(stiffness_SSE_list), color='C0', label='SSE')
plt.plot(Delta_list, np.array(stiffness_ED_list), color="C3", label='ED')
plt.errorbar(Delta_list, np.array(stiffness_SSE_list), np.array(stiffness_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel(r"$\rho_s$")
plt.title("Superfluid Density")
plt.tight_layout()
plt.savefig("Superfluid_Density_Comparison.png")
