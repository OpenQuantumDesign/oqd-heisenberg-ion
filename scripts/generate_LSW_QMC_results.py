import matplotlib.pyplot as plt
import numpy as np

from oqd_heisenberg_ion.common.postprocess import swt as lsw
from oqd_heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = []
alpha_list_2 = [round(0.0 + i * 0.1, 2) for i in range(81)]
alpha_list += alpha_list_2
print(alpha_list)
J = 1.0
T_1 = 0.005
beta_1 = J / T_1
energy_SSE_2_list = []
energy_SSE_2_err_list = []
energy_SSE_3_list = []
energy_SSE_3_err_list = []
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
num_bonds = int(N * (N - 1) / 2)
theta = 0.01 * (1 / (N - 1))
lsw_energy = []

for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    E0_lsw = lsw.E_0_LSW(N, alpha, 1.0)

    E0_mf = lsw.E_0_MF(N, alpha, 1.0, 0.0)

    norm_energy = 0.0
    for r in range(1, int((N + 1) / 2)):
        norm_energy += 1.0 / (r ** (alpha - 2))

    lsw_energy.append(E0_lsw / norm_energy)

    folder = "SSE"
    beta_1 = 1.0 / T_1

    init_start_config = 1
    start_config = init_start_config
    boundary = 1
    file_1 = "Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
        folder,
        N,
        hamiltonian_type,
        Delta,
        h,
        alpha,
        gamma,
        dist_dep_offset,
        boundary,
        T_1,
        loop_type,
        start_config,
        init_start_config,
    )
    qmc_data = np.loadtxt(file_1, delimiter=",", skiprows=2)
    energy_arr_1 = qmc_data[:, 1]

    auto_corr_drop = 2

    energy_array_2 = energy_arr_1 / E0_lsw

    energy_array_3 = energy_arr_1 / E0_mf

    energy_2 = stats.statistics_binning(energy_array_2, auto_corr_drop, eq_drop)
    energy_3 = stats.statistics_binning(energy_array_3, auto_corr_drop, eq_drop)

    energy_SSE_2_list.append(energy_2[0])
    energy_SSE_2_err_list.append(energy_2[1])

    energy_SSE_3_list.append(energy_3[0])
    energy_SSE_3_err_list.append(energy_3[1])

large_alpha_lsw_energy = lsw.E_0_LSW_NN(100 * N, 1.0)
print(large_alpha_lsw_energy / (100 * N))
delta_E = large_alpha_lsw_energy / (100 * N) + 1 / np.pi
ratio_E = (-1 / np.pi) / (large_alpha_lsw_energy / (100 * N))
print(delta_E)
print(100 * (delta_E * np.pi))

classical_energy_nn = -0.25
ratio_E_classical = (-1 / np.pi) / (classical_energy_nn)
ratio_rho_classical = (1 / np.pi) / (-classical_energy_nn)

mft_E0 = -J * (N * (N - 1)) / 8.0
lmg_E0 = -J * (N**2 - 1) / 8.0

print(lmg_E0 / N)
print(mft_E0 / N)
print(lsw_energy[0])

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
plt.scatter(alpha_list, np.array(energy_SSE_2_list), label=r"$LSW$", color="C0")
plt.errorbar(
    alpha_list, np.array(energy_SSE_2_list), np.array(energy_SSE_2_err_list), fmt="None", capsize=5, color="C0"
)
plt.scatter(alpha_list, np.array(energy_SSE_3_list), label=r"$Classical$", color="C3")
plt.errorbar(
    alpha_list, np.array(energy_SSE_3_list), np.array(energy_SSE_3_err_list), fmt="None", capsize=5, color="C3"
)
plt.plot(alpha_list, [ratio_E] * len(alpha_list), color="C0", linestyle="dashed")
plt.plot(alpha_list, [ratio_E_classical] * len(alpha_list), color="C3", linestyle="dashed")
plt.legend(prop={"size": 12})

plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$E_{0}^{Exact}/E_{0}^{Approx.}$", fontsize=16)
plt.savefig("Energy_LSW_Comparison.pdf", bbox_inches="tight", dpi=1200)
