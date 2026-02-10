import matplotlib.pyplot as plt
import numpy as np

from heisenberg_ion.common.postprocess import swt as lsw
from heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = []
alpha_list_2 = [round(0.0 + i * 0.1, 2) for i in range(81)]
alpha_list += alpha_list_2
# alpha_list.append(10.0)
print(alpha_list)
J = 1.0
T_1 = 0.005
beta_1 = J / T_1
energy_ED_list = []
mag_z_ED_list = []
stiffness_ED_list = []
energy_SSE_1_list = []
energy_SSE_1_err_list = []
energy_SSE_2_list = []
energy_SSE_2_err_list = []
energy_SSE_3_list = []
energy_SSE_3_err_list = []
mag_z_SSE_1_list = []
mag_z_SSE_1_err_list = []
stiffness_SSE_1_list = []
stiffness_SSE_1_err_list = []
stiffness_SSE_2_list = []
stiffness_SSE_2_err_list = []
stiffness_SSE_3_list = []
stiffness_SSE_3_err_list = []
stiffness_SSE_4_list = []
stiffness_SSE_4_err_list = []
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
num_bonds = int(N * (N - 1) / 2)
lsw_stiffness = []
lsw_stiffness_2 = []
theta = 0.01 * (1 / (N - 1))
pct_diff_list = []
lsw_energy = []
lsw_energy_2 = []
gamma_0_list = []
alpha_shifted_list = []
lsw_over_classical = []

for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    E0_lsw = lsw.E_0_LSW(N, alpha, 1.0)

    E0_mf = lsw.E_0_MF(N, alpha, 1.0, 0.0)

    norm = 0.0
    norm_energy = 0.0
    for r in range(1, int((N + 1) / 2)):
        norm += 1.0 / (r ** (alpha - 2))
        norm_energy += 1.0 / (r ** (alpha - 2))

    lsw_energy.append(E0_lsw / norm_energy)
    lsw_energy_2.append(E0_mf / norm_energy)

    rho_lsw_2 = lsw.rho_2(N, alpha, theta, J)
    lsw_stiffness_2.append(rho_lsw_2 / norm)

    rho_lsw = lsw.rho_mf(N, alpha, theta, J)
    lsw_stiffness.append(rho_lsw / norm)

    folder = "SSE"
    beta_1 = 1.0 / T_1

    init_start_config = 1
    start_config = init_start_config
    boundary = 1
    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
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
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(
        file_1, delimiter=",", skiprows=2, unpack=True
    )

    # drop_1 = int(step_number_1[-1]/2)
    auto_corr_drop = 2

    """
    if alpha >= 3:
        stiffness_SSE_3_list.append(-energy_SSE_1_list[alpha_list.index(round(alpha-2,2))])
        stiffness_SSE_3_err_list.append(energy_SSE_1_err_list[alpha_list.index(round(alpha-2,2))])
        alpha_shifted_list.append(alpha)
    """

    energy_array = energy_arr_1 / norm_energy
    stiffness_array = np.array(stiffness_arr_1) / norm
    s_k_pi_array = s_k_pi_arr_1
    energy_array_2 = energy_arr_1 / E0_lsw

    energy_array_3 = energy_arr_1 / E0_mf

    stiffness_array_4 = stiffness_arr_1 / rho_lsw_2
    stiffness_4 = stats.statistics_binning(stiffness_array_4, auto_corr_drop, eq_drop)
    stiffness_SSE_4_list.append(stiffness_4[0])
    stiffness_SSE_4_err_list.append(stiffness_4[1])

    stiffness_array_3 = stiffness_arr_1 / rho_lsw
    stiffness_3 = stats.statistics_binning(stiffness_array_3, auto_corr_drop, eq_drop)
    stiffness_SSE_3_list.append(stiffness_3[0])
    stiffness_SSE_3_err_list.append(stiffness_3[1])

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    energy_2 = stats.statistics_binning(energy_array_2, auto_corr_drop, eq_drop)
    energy_3 = stats.statistics_binning(energy_array_3, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)
    s_k_pi = stats.statistics_binning(s_k_pi_array, auto_corr_drop, eq_drop)

    stiffness_SSE_1_list.append(stiffness[0])
    stiffness_SSE_1_err_list.append(stiffness[1])

    energy_SSE_1_list.append(energy[0])
    energy_SSE_1_err_list.append(energy[1])
    mag_z_SSE_1_list.append(mag_z_exp[0])
    mag_z_SSE_1_err_list.append(mag_z_exp[1])

    energy_SSE_2_list.append(energy_2[0])
    energy_SSE_2_err_list.append(energy_2[1])

    energy_SSE_3_list.append(energy_3[0])
    energy_SSE_3_err_list.append(energy_3[1])

    print(step_number_1[-1])

large_alpha_lsw_energy = lsw.E_0_LSW_NN(100 * N, 1.0)
print(large_alpha_lsw_energy / (100 * N))
delta_E = large_alpha_lsw_energy / (100 * N) + 1 / np.pi
ratio_E = (-1 / np.pi) / (large_alpha_lsw_energy / (100 * N))
print(delta_E)
print(100 * (delta_E * np.pi))

norm = 0.0
alpha_new = 20
for r in range(1, int((N + 1) / 2)):
    norm += 1.0 / (r ** (alpha - 2))

large_alpha_lsw_rho = lsw.rho_2(N, alpha_new, theta, 1.0) / norm
print(large_alpha_lsw_rho)
ratio_rho_s = (1 / np.pi) / (-large_alpha_lsw_energy / (100 * N))

classical_energy_nn = -0.25
ratio_E_classical = (-1 / np.pi) / (classical_energy_nn)
ratio_rho_classical = (1 / np.pi) / (-classical_energy_nn)

mft_E0 = -J * (N * (N - 1)) / 8.0
lmg_E0 = -J * (N**2 - 1) / 8.0

print(lmg_E0 / N)
print(mft_E0 / N)
print(lsw_energy[0])

ratio_lsw_mf = []
ratio_lsw_normalized = []
N_new = 101
theta = 0.01 * (1 / (N_new - 1))
for i in range(len(alpha_list)):
    alpha = alpha_list[i]
    norm = 0.0
    J_0 = 0.0
    J_theta = 0.0
    for r in range(1, int((N_new + 1) / 2)):
        norm += 1.0 / (r ** (alpha - 2))
        J_theta = np.cos(r * theta) / (r ** (alpha))
        J_0 = 1.0 / (r ** (alpha))

    lsw_rho_large_N = lsw.rho_2(N_new, alpha, theta, 1.0, 1.0)
    mf_rho_large_N = lsw.rho_mf(N_new, alpha, theta, 1.0)

    ratio_lsw_mf.append(lsw_rho_large_N / mf_rho_large_N)
    ratio_lsw_normalized.append(lsw_rho_large_N / norm)

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
plt.scatter(alpha_list, np.array(energy_SSE_1_list), label=r"$\beta={}$".format(int(beta_1)), color="C0")
plt.errorbar(
    alpha_list, np.array(energy_SSE_1_list), np.array(energy_SSE_1_err_list), fmt="None", capsize=5, color="C0"
)
plt.scatter(alpha_list, np.array(lsw_energy), label="LSW", color="C3")
# plt.plot(alpha_list, [large_alpha_lsw_energy/(100*N)] * len(alpha_list), color='green', linestyle='dashed', label=r'LSW($\alpha \rightarrow \infty$)')
# plt.plot(alpha_list, [-1.0/np.pi] * len(alpha_list), color='black', linestyle='dashed', label=r'FFS')
# plt.plot(alpha_list, [mft_E0/N] * len(alpha_list), color='black', linestyle='dashed', label=r'MFT')
# plt.plot(alpha_list, [lmg_E0/N] * len(alpha_list), color='green', linestyle='dashed', label=r'LMG')
plt.legend(prop={"size": 12})

plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$E$", fontsize=16)
plt.savefig("Energy_LSW_Comparison_2.png", bbox_inches="tight")

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

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
plt.scatter(alpha_list, np.array(stiffness_SSE_1_list), label=r"$\beta={}$".format(int(beta_1)), color="C0")
plt.errorbar(
    alpha_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt="None", capsize=5, color="C0"
)
plt.scatter(alpha_list, np.array(lsw_stiffness_2), label="LSW", color="C3")
plt.plot(
    alpha_list,
    [large_alpha_lsw_rho] * len(alpha_list),
    color="green",
    linestyle="dashed",
    label=r"LSW($\alpha \rightarrow \infty$)",
)
plt.plot(alpha_list, [1.0 / np.pi] * len(alpha_list), color="black", linestyle="dashed", label=r"FFS")
# plt.plot(alpha_list, [mft_E0/N] * len(alpha_list), color='black', linestyle='dashed', label=r'MFT')
# plt.plot(alpha_list, [lmg_E0/N] * len(alpha_list), color='green', linestyle='dashed', label=r'LMG')
# plt.scatter(alpha_shifted_list, stiffness_SSE_3_list, label=r'E($\alpha-2$)', color='C6')
# plt.errorbar(alpha_shifted_list, np.array(stiffness_SSE_3_list), np.array(stiffness_SSE_3_err_list), fmt='None', capsize=5, color='C6')
plt.legend(prop={"size": 12})

plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$\rho_s$", fontsize=16)
plt.savefig("Stiffness_LSW_Comparison_2.png", bbox_inches="tight")

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
plt.scatter(alpha_list, np.array(stiffness_SSE_3_list), label=r"$Classical$", color="C3")
plt.errorbar(
    alpha_list, np.array(stiffness_SSE_3_list), np.array(stiffness_SSE_3_err_list), fmt="None", capsize=5, color="C3"
)
plt.scatter(alpha_list, np.array(stiffness_SSE_4_list), label=r"$LSW$", color="C0")
plt.errorbar(
    alpha_list, np.array(stiffness_SSE_4_list), np.array(stiffness_SSE_4_err_list), fmt="None", capsize=5, color="C0"
)
plt.plot(alpha_list, [ratio_rho_s] * len(alpha_list), color="C0", linestyle="dashed")
plt.plot(alpha_list, [ratio_rho_classical] * len(alpha_list), color="C3", linestyle="dashed")
plt.scatter(alpha_list, ratio_lsw_mf, color="black")
plt.legend(prop={"size": 12})

plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$\frac{\rho_s^{Exact}}{\rho_s^{Approx.}}$", fontsize=16)
plt.savefig("Stiffness_LSW_Comparison.pdf", bbox_inches="tight", dpi=1200)

plt.figure()
plt.scatter(alpha_list, ratio_lsw_normalized, color="C0")
plt.plot(alpha_list, [1.0 / np.pi] * len(alpha_list), color="C3")
plt.plot(alpha_list, [0.25] * len(alpha_list), color="black")
plt.savefig("Normalized_LSW_Stiffness.png")
