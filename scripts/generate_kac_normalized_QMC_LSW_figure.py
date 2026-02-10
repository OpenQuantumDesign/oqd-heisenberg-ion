import matplotlib.pyplot as plt
import numpy as np

from heisenberg_ion.common.postprocess import swt as lsw
from heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = []
alpha_list_2 = [0.0 + i * 0.1 for i in range(81)]
alpha_list += alpha_list_2
print(alpha_list)

J = 1.0

T_1 = 0.005
T_2 = 0.05
T_3 = 5.0

beta_1 = J / T_1
beta_2 = J / T_2
beta_3 = J / T_3

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
mag_z_SSE_2_list = []
mag_z_SSE_2_err_list = []
mag_z_SSE_3_list = []
mag_z_SSE_3_err_list = []

stiffness_SSE_1_list = []
stiffness_SSE_1_err_list = []
stiffness_SSE_2_list = []
stiffness_SSE_2_err_list = []
stiffness_SSE_3_list = []
stiffness_SSE_3_err_list = []

init_start_config = 1
start_config = init_start_config

loop_type = "deterministic"
dist_dep_offset = 0
boundary = 1
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
num_bonds = int(N * (N - 1) / 2)

lsw_stiffness = []
lsw_energy = []
mf_energy = []
mf_stiffness = []

folder = "SSE"
auto_corr_drop = 2
theta = 0.01 * (1 / (N - 1))

for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    E0_lsw = lsw.E_0_LSW(N, alpha, 1.0)
    E0_mf = lsw.E_0_MF(N, alpha, 1.0, 0.0)

    rho_lsw = lsw.rho_2(N, alpha, theta, J)
    rho_mf = lsw.rho_mf(N, alpha, theta, J)

    norm = 0.0
    norm_energy = 0.0
    for r in range(1, int((N + 1) / 2)):
        norm += 1.0 / (r ** (alpha - 2))
        # norm += N/(r**(alpha))
        norm_energy += 1.0 / (r ** (alpha))

    norm = norm_energy

    lsw_energy.append(E0_lsw / norm_energy)
    mf_energy.append(E0_mf / norm_energy)

    mf_stiffness.append(rho_mf / norm)
    lsw_stiffness.append(rho_lsw / norm)

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

    drop_1 = int(step_number_1[-1] / 2)

    print(step_number_1[-1])

    energy_array = energy_arr_1 / norm_energy
    stiffness_array = np.array(stiffness_arr_1) / norm

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

    stiffness_SSE_1_list.append(stiffness[0])
    stiffness_SSE_1_err_list.append(stiffness[1])

    energy_SSE_1_list.append(energy[0])
    energy_SSE_1_err_list.append(energy[1])

    mag_z_SSE_1_list.append(mag_z_exp[0])
    mag_z_SSE_1_err_list.append(mag_z_exp[1])

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
        N,
        hamiltonian_type,
        Delta,
        h,
        alpha,
        gamma,
        dist_dep_offset,
        boundary,
        T_2,
        loop_type,
        start_config,
        init_start_config,
    )
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(
        file_1, delimiter=",", skiprows=2, unpack=True
    )

    drop_1 = int(step_number_1[-1] / 2)
    auto_corr_drop = 2

    energy_array = energy_arr_1 / norm_energy
    stiffness_array = np.array(stiffness_arr_1) / norm

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

    stiffness_SSE_2_list.append(stiffness[0])
    stiffness_SSE_2_err_list.append(stiffness[1])

    energy_SSE_2_list.append(energy[0])
    energy_SSE_2_err_list.append(energy[1])

    mag_z_SSE_2_list.append(mag_z_exp[0])
    mag_z_SSE_2_err_list.append(mag_z_exp[1])

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
        N,
        hamiltonian_type,
        Delta,
        h,
        alpha,
        gamma,
        dist_dep_offset,
        boundary,
        T_3,
        loop_type,
        start_config,
        init_start_config,
    )
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(
        file_1, delimiter=",", skiprows=2, unpack=True
    )

    drop_1 = int(step_number_1[-1] / 2)
    auto_corr_drop = 2

    energy_array = energy_arr_1 / norm_energy
    stiffness_array = np.array(stiffness_arr_1) / norm

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

    stiffness_SSE_3_list.append(stiffness[0])
    stiffness_SSE_3_err_list.append(stiffness[1])

    energy_SSE_3_list.append(energy[0])
    energy_SSE_3_err_list.append(energy[1])

    mag_z_SSE_3_list.append(mag_z_exp[0])
    mag_z_SSE_3_err_list.append(mag_z_exp[1])

norm = 0.0
alpha_new = 50
N_new = 1001
for r in range(1, int((N_new + 1) / 2)):
    norm += 1.0 / (r ** (alpha_new - 2))
    norm_energy += 1.0 / (r ** (alpha_new))

theta = 0.01 * (1 / (1001))

large_alpha_lsw_energy = lsw.E_0_LSW_NN(N_new, 1.0)

large_alpha_lsw_rho = lsw.rho_2(N_new, alpha_new, theta, 1.0) / norm

large_alpha_lsw_rho_2 = -large_alpha_lsw_energy / (N_new)
print(large_alpha_lsw_rho_2)
print(large_alpha_lsw_rho)
print(100 * (large_alpha_lsw_rho_2 - large_alpha_lsw_rho) / large_alpha_lsw_rho)

plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
# plt.scatter(alpha_list, np.array(energy_SSE_3_list)/N, label=r'$\beta={}$'.format(beta_3), color='black')
# plt.errorbar(alpha_list, np.array(energy_SSE_3_list)/N, np.array(energy_SSE_3_err_list)/N, fmt='None', capsize=5, color='black')
plt.scatter(alpha_list, np.array(energy_SSE_2_list) / N, label=r"QMC ($\beta={}$)".format(int(beta_2)), color="C3")
plt.errorbar(
    alpha_list, np.array(energy_SSE_2_list) / N, np.array(energy_SSE_2_err_list) / N, fmt="None", capsize=5, color="C3"
)
plt.scatter(alpha_list, np.array(energy_SSE_1_list) / N, label=r"QMC ($\beta={}$)".format(int(beta_1)), color="C0")
plt.errorbar(
    alpha_list, np.array(energy_SSE_1_list) / N, np.array(energy_SSE_1_err_list) / N, fmt="None", capsize=5, color="C0"
)
plt.plot(
    alpha_list,
    [-1.0 / np.pi] * len(alpha_list),
    color="black",
    linestyle="dashed",
    label=r"Exact ($\alpha \rightarrow \infty$)",
)
plt.plot(alpha_list, np.array(mf_energy) / N, color="C4", linestyle="dashed", label="MFT")
plt.plot(
    alpha_list,
    [large_alpha_lsw_energy / (N_new)] * len(alpha_list),
    color="C5",
    linestyle="dashed",
    label=r"LSW ($\alpha \rightarrow \infty$)",
)
plt.plot(alpha_list, np.array(lsw_energy) / N, label="LSW", color="C2")
plt.legend(prop={"size": 10}, bbox_to_anchor=(0.65, 0.8))
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$E/(N \mathcal{N}_E)$", fontsize=16)
plt.savefig("Kac_normalized_Energy_per_site.pdf", bbox_inches="tight", dpi=1200)

fig, ax1 = plt.subplots()
ax1.plot(
    alpha_list,
    [1 / np.pi] * len(alpha_list),
    color="black",
    linestyle="dashed",
    label=r"Exact ($\alpha \rightarrow \infty$)",
)
# ax1.plot(alpha_list, mf_stiffness, color='C4', linestyle='dashed', label='MFT')
ax1.plot(alpha_list, lsw_stiffness, color="C2", label="LSW")
# ax1.plot(alpha_list, [large_alpha_lsw_rho_2]*len(alpha_list), color='C5', linestyle='dashed', label=r'LSW ($\alpha \rightarrow \infty$)')
ax1.scatter(
    alpha_list, np.array(stiffness_SSE_3_list), color="black", label=r"QMC ($\beta={}$)".format(round(beta_3, 1))
)
ax1.errorbar(
    alpha_list, np.array(stiffness_SSE_3_list), np.array(stiffness_SSE_3_err_list), fmt="None", capsize=5, color="black"
)
ax1.scatter(alpha_list, np.array(stiffness_SSE_2_list), color="C3", label=r"QMC ($\beta={}$)".format(int(beta_2)))
ax1.errorbar(
    alpha_list, np.array(stiffness_SSE_2_list), np.array(stiffness_SSE_2_err_list), fmt="None", capsize=5, color="C3"
)
ax1.scatter(alpha_list, np.array(stiffness_SSE_1_list), color="C0", label=r"QMC ($\beta={}$)".format(int(beta_1)))
ax1.errorbar(alpha_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt="None", capsize=5)

# ax1.legend(prop={'size': 10}, bbox_to_anchor=(0.65, 0.25))

drop_from_top = 25
drop_from_bottom = -30
left, bottom, width, height = [0.385, 0.25, 0.50, 0.37]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.scatter(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_3_list)[drop_from_top:drop_from_bottom],
    color="black",
    label=r"SSE ($\beta={}$)".format(beta_3),
    s=25,
)
ax2.errorbar(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_3_list)[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_3_err_list)[drop_from_top:drop_from_bottom],
    fmt="None",
    capsize=5,
    color="black",
)
ax2.scatter(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_2_list)[drop_from_top:drop_from_bottom],
    color="C3",
    label=r"SSE ($\beta={}$)".format(int(beta_2)),
    s=25,
)
ax2.errorbar(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_2_list)[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_2_err_list)[drop_from_top:drop_from_bottom],
    fmt="None",
    capsize=5,
    color="C3",
)
ax2.scatter(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_1_list)[drop_from_top:drop_from_bottom],
    color="C0",
    label=r"SSE ($\beta={}$)".format(int(beta_1)),
    s=25,
)
ax2.errorbar(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_1_list)[drop_from_top:drop_from_bottom],
    np.array(stiffness_SSE_1_err_list)[drop_from_top:drop_from_bottom],
    fmt="None",
    capsize=5,
)
ax2.plot(
    alpha_list[drop_from_top:drop_from_bottom],
    [1 / np.pi] * len(alpha_list[drop_from_top:drop_from_bottom]),
    color="black",
    linestyle="dashed",
)
# ax2.plot(alpha_list[drop_from_top:drop_from_bottom], [large_alpha_lsw_rho_2]*len(alpha_list[drop_from_top:drop_from_bottom]), color='C5', linestyle='dashed')
ax2.plot(alpha_list[drop_from_top:drop_from_bottom], lsw_stiffness[drop_from_top:drop_from_bottom], color="C2")
ax2.locator_params(axis="x", nbins=8)

ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$\rho_s/\mathcal{N}_E$", fontsize=16)

ax1.legend(prop={"size": 10})

plt.savefig("Kac_normalized_Superfluid_Density.pdf", bbox_inches="tight", dpi=1200)
