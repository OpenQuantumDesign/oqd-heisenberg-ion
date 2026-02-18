import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

alpha_list = [50.0]
N_list = [5, 7, 9, 11, 13]
h = 0.0
J = 1.0
T = 0.001
theta = 0.02
beta = J / T
energy_ED_list = []
mag_z_ED_list = []
stiffness_ED_list = []
Delta = 0.0
eq_drop = 0
exact_results = []
energy_exact_results = []
one_over_N_squared = []


def get_free_fermion_energy(N, phi):

    s_min = -int(N / 2)
    E = 0.0
    for i in range(N):
        s = s_min + i
        eval = -np.cos(2 * np.pi * s / N - phi)
        if eval < 0.0:
            E += eval

    return E


for k in range(len(alpha_list)):
    for j in range(len(N_list)):
        alpha = alpha_list[k]
        N = N_list[j]
        one_over_N_squared.append(1.0 / N**2)

        theta_str = round(theta, 3)

        comparison_file = "Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_0.0_Heisenberg_PBC.csv".format(
            N, Delta, h, J, J, alpha
        )
        state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)
        comp_file_1 = "Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(
            N, Delta, h, J, J, alpha, theta_str
        )
        state_index_theta, evals_theta, mag_z_theta, mag_x_theta = np.loadtxt(
            comp_file_1, skiprows=1, delimiter=",", unpack=True
        )
        comp_file_2 = "Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_theta_{}_Heisenberg_PBC.csv".format(
            N, Delta, h, J, J, alpha, -theta_str
        )
        state_index_minus_theta, evals_minus_theta, mag_z_minus_theta, mag_x_minus_theta = np.loadtxt(
            comp_file_2, skiprows=1, delimiter=",", unpack=True
        )

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
        free_energy_theta = (-1.0 / beta) * np.log(partition_theta)
        free_energy_zero = (-1.0 / beta) * np.log(partition_zero)
        free_energy_minus_theta = (-1.0 / beta) * np.log(partition_minus_theta)
        second_derivative = (free_energy_theta / N + free_energy_minus_theta / N - 2.0 * free_energy_zero / N) / (
            theta**2
        )
        second_derivative = (evals_theta[0] / N + evals_minus_theta[0] / N - 2.0 * evals[0] / N) / (theta**2)
        spin_stiffness = second_derivative
        stiffness_ED_list.append(spin_stiffness)

        energy_ED_list.append(evals[0] / N)

        E_phi = get_free_fermion_energy(N, theta)
        E_zero = get_free_fermion_energy(N, 0.0)
        E_minus_phi = get_free_fermion_energy(N, -theta)

        rho = (E_phi / N + E_minus_phi / N - 2.0 * E_zero / N) / (theta**2)

        exact_results.append(rho)
        energy_exact_results.append(E_zero / N)
        print(alpha)


def linear_func(x, a, b):

    return a * x + b


popt, pcov = sp.optimize.curve_fit(linear_func, one_over_N_squared, stiffness_ED_list)
popt_2, pcov_2 = sp.optimize.curve_fit(linear_func, one_over_N_squared, energy_ED_list)
print(popt)
one_over_N = [1.0 / n for n in N_list]

extrapolated_stiffness = []
one_over_N_extrapolated = []
extrapolated_energy = []
for i in range(N_list[0], 1000):
    one_over_N_extrapolated.append(1.0 / i)
    rho_s_extrapolated = linear_func(1.0 / i**2, popt[0], popt[1])
    E_over_N_extrapolated = linear_func(1.0 / i**2, popt_2[0], popt_2[1])
    extrapolated_stiffness.append(rho_s_extrapolated)
    extrapolated_energy.append(E_over_N_extrapolated)

print(stiffness_ED_list)
plt.figure()
plt.rcParams["mathtext.fontset"] = "stix"
plt.scatter(one_over_N, exact_results, color="black", marker="s", s=55, label="FFS")
plt.scatter(one_over_N, np.array(stiffness_ED_list), color="C0", label=r"ED ($\alpha \rightarrow \infty$)")
plt.plot(
    one_over_N_extrapolated,
    [J / np.pi] * len(one_over_N_extrapolated),
    label=r"FFS ($N \rightarrow \infty$)",
    color="black",
)
plt.plot(one_over_N_extrapolated, extrapolated_stiffness, color="C3", label="ED Fit")
plt.scatter([0.0], [popt[1]], marker="o", color="C2", label=r"ED Extrapolation")
plt.errorbar([0.0], [popt[1]], [np.sqrt(np.diag(pcov))[1]], capsize=5, color="C0", fmt="None")
plt.legend(prop={"size": 12})
plt.xlabel(r"$1/N$", fontsize=16)
plt.ylabel(r"$\rho_s$", fontsize=16)
plt.savefig("Superfluid Density ED.pdf", dpi=1200, bbox_inches="tight")

plt.figure()
plt.scatter(one_over_N, energy_exact_results, color="black", marker="s", s=55, label="FFS")
plt.scatter(one_over_N, np.array(energy_ED_list), color="C0", label=r"ED ($\alpha \rightarrow \infty$)")
plt.plot(
    one_over_N_extrapolated,
    [-J / np.pi] * len(one_over_N_extrapolated),
    label=r"FFS ($N \rightarrow \infty$)",
    color="black",
)
plt.plot(one_over_N_extrapolated, extrapolated_energy, color="C3", label="ED Fit")
plt.scatter([0.0], [popt_2[1]], marker="o", color="C2", label=r"ED Extrapolation")
plt.errorbar([0.0], [popt_2[1]], [np.sqrt(np.diag(pcov_2))[1]], capsize=5, color="C0", fmt="None")
plt.legend(prop={"size": 12})
plt.xlabel(r"$1/N$", fontsize=16)
plt.ylabel(r"$E_0/N$", fontsize=16)
plt.savefig("Energy ED.pdf", dpi=1200, bbox_inches="tight")
