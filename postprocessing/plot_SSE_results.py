import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats

N=100
Delta_list = [-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
h = 0.0
alpha=1.0
J = 1.0
T = 0.1
beta = J/T
energy_ED_list = []
mag_z_ED_list = []
energy_SSE_list = []
energy_SSE_err_list = []
mag_z_SSE_list = []
mag_z_SSE_err_list = []
stiffness_SSE_list = []
stiffness_SSE_err_list = []
auto_corr_drop=5

for Delta in Delta_list:

    if Delta < 0.0:
        gamma = 0.1
    else:
        gamma = 0.0

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_directed_loops_input_config_-1/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma)
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)

    print("\n")
    print(file_1)
    print("\n")
    energy = stats.statistics_binning(energy_arr_1, auto_corr_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop)
    stiffness = stats.statistics_binning(stiffness_arr_1, auto_corr_drop)

    energy_SSE_list.append(energy[0])
    energy_SSE_err_list.append(energy[1])

    mag_z_SSE_list.append(mag_z_exp[0])
    mag_z_SSE_err_list.append(mag_z_exp[1])

    stiffness_SSE_list.append(stiffness[0])
    stiffness_SSE_err_list.append(stiffness[1])

print(energy_SSE_list)
plt.figure()
plt.scatter(Delta_list, np.array(energy_SSE_list)/N, label="SSE", color='C0')
plt.errorbar(Delta_list, np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("E/N")
plt.title("E/N")
plt.tight_layout()
plt.savefig("Energy per site.png")

plt.figure()
plt.scatter(Delta_list, np.array(mag_z_SSE_list), color='C0', label='SSE')
plt.errorbar(Delta_list, np.array(mag_z_SSE_list), np.array(mag_z_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel("M_z ")
plt.title("Magnetization")
plt.tight_layout()
plt.savefig("Magnetization.png")

plt.figure()
plt.scatter(Delta_list, np.array(stiffness_SSE_list), color='C0', label='SSE')
plt.errorbar(Delta_list, np.array(stiffness_SSE_list), np.array(stiffness_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\Delta$")
plt.ylabel(r"$\rho_s$")
plt.title("Spin Stiffness")
plt.tight_layout()
plt.savefig("Spin Stiffness.png")
