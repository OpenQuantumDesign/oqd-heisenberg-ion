import matplotlib.pyplot as plt
import numpy as np

from oqd_heisenberg_ion.common.postprocess import utils as stats

N = 11
T_list = [
    "5.0",
    "2.0",
    "1.0",
    "0.5",
    "0.333",
    "0.25",
    "0.2",
    "0.167",
    "0.143",
    "0.122",
    "0.111",
    "0.1",
    "0.0909",
    "0.0833",
    "0.0769",
    "0.0714",
    "0.0666",
    "0.0625",
    "0.0588",
    "0.0555",
    "0.0526",
    "0.05",
]
beta_list = [1 / float(T) for T in T_list]

energy_ED_list = []
stiffness_ED_list = []

energy_SSE_list = []
energy_SSE_err_list = []
stiffness_SSE_list = []
stiffness_SSE_err_list = []

eq_drop = 0

Delta = 0.0
hamiltonian_type = int(Delta)
h = 0.0
alpha = 10.0
alpha_identifier_locs = {1.0: [8.75, 0.15], 10.0: [8.6, 0.01]}
gamma = 0.0
loop_type = "deterministic"
dist_dep_offset = 0
J = 1.0
theta = 0.02
auto_corr_drop = 2


def ed_properties(beta, evals, evals_theta, evals_minus_theta, energy_list, stiffness_list):

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

    energy_theta /= partition_theta
    energy_zero /= partition_zero
    energy_minus_theta /= partition_minus_theta
    free_energy_theta = (-1.0 / beta) * np.log(partition_theta)
    free_energy_zero = (-1.0 / beta) * np.log(partition_zero)
    free_energy_minus_theta = (-1.0 / beta) * np.log(partition_minus_theta)
    second_derivative = (free_energy_theta / N + free_energy_minus_theta / N - 2.0 * free_energy_zero / N) / (theta**2)
    second_derivative = (free_energy_theta / N + free_energy_minus_theta / N - 2.0 * free_energy_zero / N) / (theta**2)
    spin_stiffness = second_derivative

    energy_list.append(energy_zero / N)
    stiffness_list.append(spin_stiffness)


for i in range(len(T_list)):
    T = T_list[i]
    beta = beta_list[i]
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
        T,
        loop_type,
        start_config,
        init_start_config,
    )
    qmc_data = np.loadtxt(file_1, delimiter=",", skiprows=2)

    energy_array = qmc_data[:, 1]
    stiffness_array = qmc_data[:, 3]

    energy = stats.statistics_binning(energy_array, auto_corr_drop, eq_drop)
    stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)

    stiffness_SSE_list.append(stiffness[0])
    stiffness_SSE_err_list.append(stiffness[1])

    energy_SSE_list.append(energy[0] / N)
    energy_SSE_err_list.append(energy[1] / N)

comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_0.0_Heisenberg_PBC.csv".format(
    N, Delta, h, J, J, alpha
)
state_index, evals = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)

comp_file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(
    N, Delta, h, J, J, alpha, theta
)
state_index_theta, evals_theta, mag_z_theta, mag_x_theta = np.loadtxt(
    comp_file_1, skiprows=1, delimiter=",", unpack=True
)

comp_file_2 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(
    N, Delta, h, J, J, alpha, -theta
)
state_index_minus_theta, evals_minus_theta, mag_z_minus_theta, mag_x_minus_theta = np.loadtxt(
    comp_file_2, skiprows=1, delimiter=",", unpack=True
)

beta_ED_list = np.array(beta_list)
for i in range(len(beta_ED_list)):
    beta = beta_ED_list[i]

    ed_properties(beta, evals, evals_theta, evals_minus_theta, energy_ED_list, stiffness_ED_list)


beta_ED_list_2 = np.linspace(beta_list[0], beta_list[-1], 100, True)
stiffness_ED_list_2 = []
energy_ED_list_2 = []
for i in range(len(beta_ED_list_2)):
    beta = beta_ED_list_2[i]

    ed_properties(beta, evals, evals_theta, evals_minus_theta, energy_ED_list_2, stiffness_ED_list_2)

fig, ax1 = plt.subplots()
plt.rcParams["mathtext.fontset"] = "stix"
if alpha == 10.0:
    left, bottom, width, height = [0.5, 0.3, 0.33, 0.33]
else:
    left, bottom, width, height = [0.47, 0.35, 0.33, 0.33]
ax2 = fig.add_axes([left, bottom, width, height])

ax1.errorbar(beta_list, stiffness_SSE_list, stiffness_SSE_err_list, capsize=5, fmt="None")
ax1.scatter(beta_ED_list, stiffness_ED_list, color="C3", label="ED")
ax1.plot(beta_ED_list_2, stiffness_ED_list_2, color="C3")
ax1.scatter(beta_list, stiffness_SSE_list, color="C0", label="SSE")
ax1.legend(prop={"size": 12})
ax1.text(
    alpha_identifier_locs[alpha][0],
    alpha_identifier_locs[alpha][1],
    r"$\alpha={}$".format(alpha),
    fontsize=14,
    bbox=dict(facecolor="w"),
)
ax1.set_xlabel(r"$\beta$", fontsize=16)
ax1.set_ylabel(r"$\rho_s$", fontsize=16)

ax2.errorbar(beta_list, energy_SSE_list, energy_SSE_err_list, capsize=5, fmt="None")
ax2.scatter(beta_ED_list, energy_ED_list, color="C3", label="ED")
ax2.plot(beta_ED_list_2, energy_ED_list_2, color="C3")
ax2.scatter(beta_list, energy_SSE_list, color="C0", label="SSE")
ax2.set_xlabel(r"$\beta$")
ax2.set_ylabel(r"$E/N$")
plt.savefig("beta_convergence_alpha_{}.pdf".format(alpha), dpi=1200, bbox_inches="tight")
