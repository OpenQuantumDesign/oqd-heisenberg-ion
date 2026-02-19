import matplotlib.pyplot as plt
import numpy as np

from oqd_heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = [0.0 + i * 0.1 for i in range(81)]
print(alpha_list)
J = 1.0
T_1 = 0.005
T_2 = 0.05
T_3 = 5.0
beta_1 = J / T_1
beta_2 = J / T_2
beta_3 = J / T_3

energy_SSE_1_list = []
energy_SSE_1_err_list = []
energy_SSE_2_list = []
energy_SSE_2_err_list = []
energy_SSE_3_list = []
energy_SSE_3_err_list = []

stiffness_SSE_1_list = []
stiffness_SSE_1_err_list = []
stiffness_SSE_2_list = []
stiffness_SSE_2_err_list = []
stiffness_SSE_3_list = []
stiffness_SSE_3_err_list = []

start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
boundary = 1
folder = "SSE"
init_start_config = 1
start_config = 1


def process_qmc_data(T, alpha, energy_list, energy_err_list, stiffness_list, stiffness_err_list, norm=1.0):

    # Updates input lists in place

    file = "Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
        folder,
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

    qmc_data = np.loadtxt(file, delimiter=",", skiprows=2)

    energy_array = qmc_data[:, 1] / norm
    stiffness_array = qmc_data[:, 3] / norm

    energy = stats.statistics_binning(energy_array, 2, 0)
    stiffness = stats.statistics_binning(stiffness_array, 2, 0)

    stiffness_list.append(stiffness[0])
    stiffness_err_list.append(stiffness[1])

    energy_list.append(energy[0])
    energy_err_list.append(energy[1])


for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    process_qmc_data(
        T_1, alpha, energy_SSE_1_list, energy_SSE_1_err_list, stiffness_SSE_1_list, stiffness_SSE_1_err_list
    )
    process_qmc_data(
        T_2, alpha, energy_SSE_2_list, energy_SSE_2_err_list, stiffness_SSE_2_list, stiffness_SSE_2_err_list
    )
    process_qmc_data(
        T_3, alpha, energy_SSE_3_list, energy_SSE_3_err_list, stiffness_SSE_3_list, stiffness_SSE_3_err_list
    )

print(np.array(stiffness_SSE_1_list))
plt.rcParams["mathtext.fontset"] = "stix"
fig, ax1 = plt.subplots()

ax1.scatter(alpha_list, np.array(energy_SSE_3_list) / N, label=r"$\beta={}$".format(beta_3), color="black")
ax1.errorbar(
    alpha_list,
    np.array(energy_SSE_3_list) / N,
    np.array(energy_SSE_3_err_list) / N,
    fmt="None",
    capsize=5,
    color="black",
)
ax1.scatter(alpha_list, np.array(energy_SSE_1_list) / N, label=r"$\beta={}$".format(int(beta_1)), color="C0")
ax1.errorbar(
    alpha_list, np.array(energy_SSE_1_list) / N, np.array(energy_SSE_1_err_list) / N, fmt="None", capsize=5, color="C0"
)
ax1.plot(alpha_list, [-1.0 / np.pi] * len(alpha_list), color="green", linestyle="dashed")

drop_from_top = 10
drop_from_bottom = -30
left, bottom, width, height = [0.39, 0.38, 0.49, 0.37]
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

fig, ax1 = plt.subplots()
ax1.scatter(alpha_list, np.array(stiffness_SSE_3_list), color="black", label=r"$\beta={}$".format(beta_3))
ax1.errorbar(
    alpha_list,
    np.array(stiffness_SSE_3_list),
    np.array(stiffness_SSE_3_err_list),
    fmt="None",
    capsize=5,
    color="black",
)
ax1.scatter(alpha_list, np.array(stiffness_SSE_2_list), color="C3", label=r"$\beta={}$".format(int(beta_2)))
ax1.errorbar(
    alpha_list, np.array(stiffness_SSE_2_list), np.array(stiffness_SSE_2_err_list), fmt="None", capsize=5, color="C3"
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
