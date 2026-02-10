import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

from heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = []
# alpha_list.append(0.0)
alpha_list_2 = [0.0 + i * 0.1 for i in range(81)]
alpha_list += alpha_list_2
alpha_list_3 = [0.0 + i * 0.1 for i in range(81)]
# alpha_list.append(10.0)
alpha_list = alpha_list_3
print(alpha_list)
J = 1.0
T_1 = 0.005
T_2 = 0.05
T_3 = 5.0
# T_4 = 0.1
beta_1 = J / T_1
beta_2 = J / T_2
beta_3 = J / T_3
# beta_4 = J/T_4
energy_ED_list = []
mag_z_ED_list = []
stiffness_ED_list = []
energy_SSE_1_list = []
energy_SSE_1_err_list = []
energy_SSE_2_list = []
energy_SSE_2_err_list = []
energy_SSE_3_list = []
energy_SSE_3_err_list = []
energy_SSE_4_list = []
energy_SSE_4_err_list = []
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
stiffness_SSE_4_list = []
stiffness_SSE_4_err_list = []
s_k_SSE_list = []
s_k_SSE_err_list = []
theta = 0.1
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
num_bonds = int(N * (N - 1) / 2)
lsw_stiffness = []
theta = 0.01
pct_diff_list = []
lsw_energy = []
gamma_0_list = []

for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    norm = 1.0

    folder = "SSE"

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

    drop_1 = int(step_number_1[-1] / 2)
    auto_corr_drop = 2

    if step_number_1[-1] < 1000000:
        print("alpha = ", alpha)
        print("T = ", T_1)
        print("steps = ", step_number_1[-1])

    if hamiltonian_type == 0:
        energy_array = energy_arr_1
        stiffness_array = np.array(stiffness_arr_1) / norm
        # stiffness_array = np.array(stiffness_arr_1)
        s_k_pi_array = s_k_pi_arr_1

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)
    s_k_pi = stats.statistics_binning(s_k_pi_array, auto_corr_drop, eq_drop)

    stiffness_SSE_1_list.append(stiffness[0])
    stiffness_SSE_1_err_list.append(stiffness[1])

    energy_SSE_1_list.append(energy[0])
    energy_SSE_1_err_list.append(energy[1])
    mag_z_SSE_1_list.append(mag_z_exp[0])
    mag_z_SSE_1_err_list.append(mag_z_exp[1])

    init_start_config = 1
    start_config = init_start_config
    boundary = 1
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

    if step_number_1[-1] < 1000000:
        print("alpha = ", alpha)
        print("T = ", T_2)
        print("steps = ", step_number_1[-1])

    energy_array = energy_arr_1
    stiffness_array = stiffness_arr_1

    if hamiltonian_type == 0:
        energy_array = energy_arr_1
        stiffness_array = np.array(stiffness_arr_1) / norm
        s_k_pi_array = s_k_pi_arr_1

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)
    s_k_pi = stats.statistics_binning(s_k_pi_array, auto_corr_drop, eq_drop)

    stiffness_SSE_2_list.append(stiffness[0])
    stiffness_SSE_2_err_list.append(stiffness[1])

    energy_SSE_2_list.append(energy[0])
    energy_SSE_2_err_list.append(energy[1])

    mag_z_SSE_2_list.append(mag_z_exp[0])
    mag_z_SSE_2_err_list.append(mag_z_exp[1])

    init_start_config = 1
    start_config = init_start_config
    boundary = 1
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

    if step_number_1[-1] < 1000000:
        print("alpha = ", alpha)
        print("T = ", T_3)
        print("steps = ", step_number_1[-1])

    energy_array = energy_arr_1
    stiffness_array = stiffness_arr_1

    if hamiltonian_type == 0:
        energy_array = energy_arr_1
        stiffness_array = np.array(stiffness_arr_1) / norm
        # stiffness_array = np.array(stiffness_arr_1)
        s_k_pi_array = s_k_pi_arr_1

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)
    s_k_pi = stats.statistics_binning(s_k_pi_array, auto_corr_drop, eq_drop)

    stiffness_SSE_3_list.append(stiffness[0])
    stiffness_SSE_3_err_list.append(stiffness[1])

    energy_SSE_3_list.append(energy[0])
    energy_SSE_3_err_list.append(energy[1])

    mag_z_SSE_3_list.append(mag_z_exp[0])
    mag_z_SSE_3_err_list.append(mag_z_exp[1])


print(gamma_0_list)
print(np.array(stiffness_SSE_1_list))
plt.rcParams["mathtext.fontset"] = "stix"
fig, ax1 = plt.subplots()

ax1.scatter(alpha_list_3, np.array(energy_SSE_3_list) / N, label=r"$\beta={}$".format(beta_3), color="black")
ax1.errorbar(
    alpha_list_3,
    np.array(energy_SSE_3_list) / N,
    np.array(energy_SSE_3_err_list) / N,
    fmt="None",
    capsize=5,
    color="black",
)
# ax1.scatter(alpha_list_3, np.array(energy_SSE_2_list)/N, label=r'$\beta={}$'.format(int(beta_2)), color='C3')
# ax1.errorbar(alpha_list_3, np.array(energy_SSE_2_list)/N, np.array(energy_SSE_2_err_list)/N, fmt='None', capsize=5, color='C3')
ax1.scatter(alpha_list, np.array(energy_SSE_1_list) / N, label=r"$\beta={}$".format(int(beta_1)), color="C0")
ax1.errorbar(
    alpha_list, np.array(energy_SSE_1_list) / N, np.array(energy_SSE_1_err_list) / N, fmt="None", capsize=5, color="C0"
)
ax1.plot(alpha_list, [-1.0 / np.pi] * len(alpha_list), color="green", linestyle="dashed")

drop_from_top = 10
drop_from_bottom = -30
left, bottom, width, height = [0.39, 0.38, 0.49, 0.37]
# left,bottom,width,height = [1.1, 0.0, 0.8, 0.8]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.scatter(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(energy_SSE_3_list)[drop_from_top:drop_from_bottom] / N,
    color="black",
    label=r"SSE ($\beta={}$)".format(beta_3),
    s=25,
)
ax2.errorbar(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(energy_SSE_3_list)[drop_from_top:drop_from_bottom] / N,
    np.array(energy_SSE_3_err_list)[drop_from_top:drop_from_bottom] / N,
    fmt="None",
    capsize=5,
    color="black",
)
# ax2.scatter(alpha_list[drop_from_top:drop_from_bottom], np.array(energy_SSE_2_list)[drop_from_top:drop_from_bottom]/N, color='C3', label=r'SSE ($\beta={}$)'.format(int(beta_2)), s=25)
# ax2.errorbar(alpha_list[drop_from_top:drop_from_bottom], np.array(energy_SSE_2_list)[drop_from_top:drop_from_bottom]/N, np.array(energy_SSE_2_err_list)[drop_from_top:drop_from_bottom]/N, fmt='None', capsize=5, color='C3')
ax2.scatter(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(energy_SSE_1_list)[drop_from_top:drop_from_bottom] / N,
    color="C0",
    label=r"SSE ($\beta={}$)".format(int(beta_1)),
    s=25,
)
ax2.errorbar(
    alpha_list[drop_from_top:drop_from_bottom],
    np.array(energy_SSE_1_list)[drop_from_top:drop_from_bottom] / N,
    np.array(energy_SSE_1_err_list)[drop_from_top:drop_from_bottom] / N,
    fmt="None",
    capsize=5,
)
ax2.plot(
    alpha_list[drop_from_top:drop_from_bottom],
    [-1 / np.pi] * len(alpha_list[drop_from_top:drop_from_bottom]),
    color="black",
    linestyle="dashed",
)
ax2.plot(
    alpha_list[drop_from_top:drop_from_bottom],
    [-1 / np.pi] * len(alpha_list[drop_from_top:drop_from_bottom]),
    color="green",
    linestyle="dashed",
)
ax2.locator_params(axis="x", nbins=8)
ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$E/N$", fontsize=16)

ax2.set_xlabel(r"$\alpha$", fontsize=14)
ax2.set_ylabel(r"$E/N$", fontsize=14)

ax1.legend(prop={"size": 12})

plt.savefig("Energy_per_site.pdf", bbox_inches="tight", dpi=1200)

plt.figure()
plt.scatter(alpha_list_3, np.array(mag_z_SSE_3_list), label=r"$\beta={}$".format(beta_3), color="black")
plt.errorbar(
    alpha_list_3, np.array(mag_z_SSE_3_list), np.array(mag_z_SSE_3_err_list), fmt="None", capsize=5, color="black"
)
plt.scatter(alpha_list_3, np.array(mag_z_SSE_2_list), label=r"$\beta={}$".format(int(beta_2)), color="C3")
plt.errorbar(
    alpha_list_3, np.array(mag_z_SSE_2_list), np.array(mag_z_SSE_2_err_list), fmt="None", capsize=5, color="C3"
)
plt.scatter(alpha_list, np.array(mag_z_SSE_1_list), label=r"$\beta={}$".format(int(beta_1)), color="C0")
plt.errorbar(alpha_list, np.array(mag_z_SSE_1_list), np.array(mag_z_SSE_1_err_list), fmt="None", capsize=5, color="C0")
plt.plot(alpha_list, [0.0] * len(alpha_list), color="green", linestyle="dashed")

formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-4, -4))
plt.gca().yaxis.set_major_formatter(formatter)
plt.legend(prop={"size": 12})
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$M_z$", fontsize=16)
plt.savefig("Magnetization.pdf", bbox_inches="tight", dpi=1200)

fig, ax1 = plt.subplots()
ax1.scatter(alpha_list_3, np.array(stiffness_SSE_3_list), color="black", label=r"$\beta={}$".format(beta_3))
ax1.errorbar(
    alpha_list_3,
    np.array(stiffness_SSE_3_list),
    np.array(stiffness_SSE_3_err_list),
    fmt="None",
    capsize=5,
    color="black",
)
ax1.scatter(alpha_list_3, np.array(stiffness_SSE_2_list), color="C3", label=r"$\beta={}$".format(int(beta_2)))
ax1.errorbar(
    alpha_list_3, np.array(stiffness_SSE_2_list), np.array(stiffness_SSE_2_err_list), fmt="None", capsize=5, color="C3"
)
ax1.scatter(alpha_list, np.array(stiffness_SSE_1_list), color="C0", label=r"$\beta={}$".format(int(beta_1)))
ax1.errorbar(alpha_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt="None", capsize=5)
ax1.legend(prop={"size": 12})
drop_from_top = 25
drop_from_bottom = -30
left, bottom, width, height = [0.39, 0.28, 0.49, 0.37]
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
ax2.plot(
    alpha_list[drop_from_top:drop_from_bottom],
    [1 / np.pi] * len(alpha_list[drop_from_top:drop_from_bottom]),
    color="green",
    linestyle="dashed",
)
ax2.locator_params(axis="x", nbins=8)

ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$\rho_s$", fontsize=16)

ax2.set_xlabel(r"$\alpha$", fontsize=14)
ax2.set_ylabel(r"$\rho_s$", fontsize=14)

plt.savefig("Superfluid_Density.pdf", bbox_inches="tight", dpi=1200)
