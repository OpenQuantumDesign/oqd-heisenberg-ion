import matplotlib.pyplot as plt
import numpy as np

from heisenberg_ion.common.postprocess import swt as lsw
from heisenberg_ion.common.postprocess import utils as stats

N = 101
h = 0.0
alpha_list = [0.0 + i * 0.1 for i in range(81)]

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

init_start_config = 1
start_config = 1

loop_type = "deterministic"
dist_dep_offset = 0
boundary = 1
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 10000

lsw_stiffness = []
lsw_energy = []
mf_energy = []
mf_stiffness = []

folder = "SSE"
auto_corr_drop = 2
theta = 0.01 * (1 / (N - 1))


def process_qmc_data(T, alpha, stiffness_list, stiffness_err_list, norm_stiffness):

    # Updates input lists in place

    file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(
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

    stiffness_array = qmc_data[:, 3] / norm_stiffness

    stiffness = stats.statistics_binning(stiffness_array, 2, 0)

    stiffness_list.append(stiffness[0])
    stiffness_err_list.append(stiffness[1])


for i in range(len(alpha_list)):
    alpha = round(alpha_list[i], 1)

    E0_lsw = lsw.E_0_LSW(N, alpha, 1.0)
    E0_mf = lsw.E_0_MF(N, alpha, 1.0, 0.0)

    rho_lsw = lsw.rho_2(N, alpha, theta, J)
    rho_mf = lsw.rho_mf(N, alpha, theta, J)

    norm_stiffness = 0.0
    norm_energy = 0.0
    for r in range(1, int((N + 1) / 2)):
        norm_stiffness += 1.0 / (r ** (alpha - 2))
        norm_energy += 1.0 / (r ** (alpha))

    lsw_energy.append(E0_lsw / norm_energy)
    mf_energy.append(E0_mf / norm_energy)

    mf_stiffness.append(rho_mf / norm_stiffness)
    lsw_stiffness.append(rho_lsw / norm_stiffness)

    process_qmc_data(
        T_1,
        alpha,
        stiffness_SSE_1_list,
        stiffness_SSE_1_err_list,
        norm_stiffness,
    )

    process_qmc_data(
        T_2,
        alpha,
        stiffness_SSE_2_list,
        stiffness_SSE_2_err_list,
        norm_stiffness,
    )

    process_qmc_data(
        T_3,
        alpha,
        stiffness_SSE_3_list,
        stiffness_SSE_3_err_list,
        norm_stiffness,
    )

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


def get_free_fermion_energy(N, phi):

    s_min = -int(N / 2)
    E = 0.0
    for i in range(N):
        s = s_min + i
        eval = -np.cos(2 * np.pi * s / N - phi)
        if eval < 0.0:
            E += eval

    return E


theta = 0.0001
E_phi_ff = get_free_fermion_energy(101, theta)
E_zero_ff = get_free_fermion_energy(101, 0.0)
E_minus_phi_ff = get_free_fermion_energy(101, -theta)

rho_ff = (E_phi_ff / N + E_minus_phi_ff / N - 2.0 * E_zero_ff / N) / (theta**2)
print(rho_ff)
print(1.0 / np.pi)

fig, ax1 = plt.subplots()
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
ax1.plot(
    alpha_list,
    [1.0 / np.pi] * len(alpha_list),
    color="black",
    linestyle="dashed",
    label=r"Exact ($\alpha \rightarrow \infty$)",
)
ax1.plot(alpha_list, mf_stiffness, color="C4", linestyle="dashed", label="MFT")
ax1.plot(alpha_list, lsw_stiffness, color="C2", label="LSW")
ax1.plot(
    alpha_list,
    [large_alpha_lsw_rho_2] * len(alpha_list),
    color="C5",
    linestyle="dashed",
    label=r"LSW ($\alpha \rightarrow \infty$)",
)
ax1.legend(prop={"size": 10}, bbox_to_anchor=(0.65, 0.25))

ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$\rho_s/\mathcal{N}_\rho$", fontsize=16)

plt.savefig("Renormalized_Superfluid_Density.pdf", bbox_inches="tight", dpi=1200)
