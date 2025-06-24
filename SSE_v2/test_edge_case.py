import os
import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats

N=3
h = 0.0
alpha=9.0
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
loop_type = "deterministic"
dist_dep_offset = 0
hamiltonian_types = [-1,0,1]
Delta_list = [-1.0, 0.0, 1.0]
gamma = 0.0
auto_corr_drop = 1
eq_drop = 0

for i in range(len(hamiltonian_types)):

    Delta = Delta_list[i]
    hamiltonian_type = hamiltonian_types[i]
    if hamiltonian_type == -1:
        start_config = 3
    else:
        start_config = 1
    init_start_config = start_config

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)

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

    energy_ED_list.append(energy_zero)
    mag_z_ED_list.append(ED_magnetization)
    energy_SSE_list.append(energy[0])
    energy_SSE_err_list.append(energy[1])
    mag_z_SSE_list.append(mag_z_exp[0])
    mag_z_SSE_err_list.append(mag_z_exp[1])

    print(step_number_1[-1])

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