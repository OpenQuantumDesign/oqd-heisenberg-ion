import matplotlib.pyplot as plt
import numpy as np

from heisenberg_ion.common.postprocess import utils as stats

N_list = [3, 5, 7, 9, 13, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125]
h = 0.0
alpha = 0.0
J = 1.0
T = 0.05
beta = J / T
energy_SSE_1_list = []
energy_SSE_1_err_list = []
energy_SSE_2_list = []
energy_SSE_2_err_list = []
stiffness_SSE_1_list = []
stiffness_SSE_1_err_list = []
stiffness_SSE_2_list = []
stiffness_SSE_2_err_list = []
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
one_over_N_squared = []
exact_results = []

for i in range(len(N_list)):
    N = N_list[i]
    one_over_N_squared.append(1.0 / N**2)

    init_start_config = 1
    start_config = init_start_config
    boundary = 1

    gamma_0 = (N**2 - 1) / 24.0
    gamma_0_0 = N * gamma_0 / 4.0
    mft_energy = -N * (N - 1) / 8.0

    lmg_stiffness = (N**2 - 1) * (N + 1) / (24.0 * 4.0)
    lmg_energy = -(N**2 - 1) / 8.0

    file_1 = "Results/SSE_18_08_25/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
        N,
        hamiltonian_type,
        Delta,
        h,
        alpha,
        gamma,
        dist_dep_offset,
        boundary,
        T,
        loop_type,
        start_config,
        init_start_config,
    )
    qmc_data = np.loadtxt(file_1, delimiter=",", skiprows=2)
    energy_arr_1 = qmc_data[:, 1]
    stiffness_arr_1 = qmc_data[:, 3]

    auto_corr_drop = 2

    energy_array_1 = (energy_arr_1) / (mft_energy)
    energy_array_2 = (energy_arr_1) / lmg_energy
    stiffness_array_1 = (stiffness_arr_1) / gamma_0_0
    stiffness_array_2 = (stiffness_arr_1) / lmg_stiffness

    energy_1 = stats.statistics_binning(energy_array_1, auto_corr_drop, eq_drop)
    energy_2 = stats.statistics_binning(energy_array_2, auto_corr_drop, eq_drop)
    stiffness_1 = stats.statistics_binning(stiffness_array_1, auto_corr_drop, eq_drop)
    stiffness_2 = stats.statistics_binning(stiffness_array_2, auto_corr_drop, eq_drop)

    stiffness_SSE_1_list.append(stiffness_1[0])
    stiffness_SSE_1_err_list.append(stiffness_1[1])

    stiffness_SSE_2_list.append(stiffness_2[0])
    stiffness_SSE_2_err_list.append(stiffness_2[1])

    energy_SSE_1_list.append(energy_1[0])
    energy_SSE_1_err_list.append(energy_1[1])

    energy_SSE_2_list.append(energy_2[0])
    energy_SSE_2_err_list.append(energy_2[1])

N_exact_list = []
for n in range(63):
    N_n = N_list[0] + 2 * n
    N_exact_list.append(N_n)

print(energy_SSE_1_list)
print(energy_SSE_1_err_list)
print("\n")
print(stiffness_SSE_2_list)
print(stiffness_SSE_2_err_list)
print("\n")

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
fig, ax1 = plt.subplots()
left, bottom, width, height = [0.32, 0.47, 0.36, 0.36]
ax2 = fig.add_axes([left, bottom, width, height])

ax2.scatter(N_list, np.array(stiffness_SSE_2_list), color="C0", label="SSE")
ax2.errorbar(N_list, np.array(stiffness_SSE_2_list), np.array(stiffness_SSE_2_err_list), fmt="None", capsize=5)
ax2.plot(N_exact_list, [1.0] * len(N_exact_list), color="black", linestyle="dashed")
ax2.set_xlabel(r"$N$", fontsize=10)
ax2.set_ylabel(r"$\overline{\rho_s}/\rho_s$", fontsize=12)

ax1.plot(N_exact_list, [1 + 1 / N for N in N_exact_list], color="C3", label="Exact")
ax1.scatter(N_list, [1 + 1 / N for N in N_list], color="C3")
ax1.scatter(N_list, np.array(stiffness_SSE_1_list), label="SSE", color="C0")
ax1.errorbar(N_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt="None", capsize=5)
ax1.plot(N_exact_list, [1.0] * len(N_exact_list), color="black", linestyle="dashed")
ax1.legend(prop={"size": 12})
ax1.set_xlabel(r"$N$", fontsize=16)
ax1.set_ylabel(r"$\overline{\rho_s}/\rho_s^{MF}$", fontsize=16)

plt.savefig("Superfluid_Density_MFT.pdf", dpi=1200, bbox_inches="tight")

plt.figure()
fig, ax1 = plt.subplots()
left, bottom, width, height = [0.44, 0.36, 0.36, 0.36]
ax2 = fig.add_axes([left, bottom, width, height])

ax2.scatter(N_list, np.array(energy_SSE_2_list), color="C0", label="SSE")
ax2.errorbar(N_list, np.array(energy_SSE_2_list), np.array(energy_SSE_2_err_list), fmt="None", capsize=5)
ax2.plot(N_exact_list, [1.0] * len(N_exact_list), color="black", linestyle="dashed")
ax2.set_xlabel(r"$N$", fontsize=10)
ax2.set_ylabel(r"$\overline{E_0}/E_0$", fontsize=12)

ax1.plot(N_exact_list, [1 + 1 / N for N in N_exact_list], color="C3", label="Exact")
ax1.scatter(N_list, [1 + 1 / N for N in N_list], color="C3")
ax1.scatter(N_list, np.array(energy_SSE_1_list), label="SSE", color="C0")
ax1.errorbar(N_list, np.array(energy_SSE_1_list), np.array(energy_SSE_1_err_list), fmt="None", capsize=5)
ax1.plot(N_exact_list, [1.0] * len(N_exact_list), color="black", linestyle="dashed")
ax1.legend(prop={"size": 12})
ax1.set_xlabel(r"$N$", fontsize=16)
ax1.set_ylabel(r"$\overline{E_0}/E_0^{MF}$", fontsize=16)

plt.savefig("Energy_MFT.pdf", dpi=1200, bbox_inches="tight")
