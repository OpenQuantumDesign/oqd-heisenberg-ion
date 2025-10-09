import numpy as np
import matplotlib.pyplot as plt
import statistical_analysis as stats

N_list = [100]
h = 0.0
#alpha_list=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
alpha = 1.0
J = 1.0
T_list = [0.01, 0.008, 0.005, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0]
energy_SSE_list = np.zeros((len(N_list), len(T_list)))
energy_SSE_err_list = np.zeros((len(N_list), len(T_list)))
mag_z_SSE_list = np.zeros((len(N_list), len(T_list)))
mag_z_SSE_err_list = np.zeros((len(N_list), len(T_list)))
stiffness_SSE_list = np.zeros((len(N_list), len(T_list)))
stiffness_SSE_err_list = np.zeros((len(N_list), len(T_list)))
s_k_SSE_list = np.zeros((len(N_list), len(T_list)))
s_k_SSE_err_list = np.zeros((len(N_list), len(T_list)))
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 1
one_over_N_squared = []
exact_stiffness = []
exact_energy = []
init_start_config = 1
start_config = init_start_config
boundary=1
auto_corr_drop = 2

for i in range(len(N_list)):

    N = N_list[i]
    one_over_N_squared.append(1.0/N**2)
    exact_stiffness.append(1.0/(N*np.sin(np.pi/N)))
    exact_energy.append(-1.0/(N*np.sin(np.pi/N)))

    for j in range(len(T_list)):

        T = T_list[j]

        beta = J/T
    
        file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T, loop_type, start_config, init_start_config)
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
        
        energy_array = energy_arr_1
        stiffness_array = stiffness_arr_1/1.5
        s_k_pi_array = s_k_pi_arr_1

        energy = stats.statistics_binning(energy_array/N, auto_corr_drop, eq_drop)
        mag_z_exp = stats.statistics_binning(magnetization_arr_1, auto_corr_drop, eq_drop)
        stiffness = stats.statistics_binning(stiffness_array, auto_corr_drop, eq_drop)
        s_k_pi = stats.statistics_binning(s_k_pi_array, auto_corr_drop, eq_drop)

        energy_SSE_list[i,j] = energy[0]
        energy_SSE_err_list[i,j] = energy[1]
        mag_z_SSE_list[i,j] = mag_z_exp[0]
        mag_z_SSE_err_list[i,j] = mag_z_exp[1]
        stiffness_SSE_list[i,j] = stiffness[0]
        stiffness_SSE_err_list[i,j] = stiffness[1]
        s_k_SSE_list[i,j] = s_k_pi[0]
        s_k_SSE_err_list[i,j] = s_k_pi[1]

        print(step_number_1[-1])

    plt.figure()
    plt.scatter(T_list, energy_SSE_list[i,:], label="SSE", color='C0')
    plt.errorbar(T_list, energy_SSE_list[i,:], energy_SSE_err_list[i,:], fmt='None', capsize=5)
    plt.legend()
    plt.xlabel(r"$T$")
    plt.ylabel("E/N")
    plt.title("E/N vs T (N={})".format(N))
    plt.tight_layout()
    plt.savefig("Energy_per_site_N_{}.png".format(N))

    plt.figure()
    plt.scatter(T_list, mag_z_SSE_list[i,:], color='C0', label='SSE')
    plt.errorbar(T_list, mag_z_SSE_list[i,:], mag_z_SSE_err_list[i,:], fmt='None', capsize=5)
    plt.legend()
    plt.xlabel(r"$T$")
    plt.ylabel(r"$M_z$")
    plt.title("Magnetization vs T (N = {})".format(N))
    plt.tight_layout()
    plt.savefig("Magnetization_N_{}.png".format(N))

    plt.figure()
    plt.scatter(T_list, stiffness_SSE_list[i,:], color='C0', label='SSE')
    plt.errorbar(T_list, stiffness_SSE_list[i,:], stiffness_SSE_err_list[i,:], fmt='None', capsize=5)
    plt.legend()
    plt.xlabel(r"$T$")
    plt.ylabel(r"$\rho_s$")
    plt.title(r"Superfluid Density ($\Delta$={}, $\alpha$={})".format(Delta, alpha))
    plt.tight_layout()
    plt.savefig("Superfluid_Density_N_{}.png".format(N))
    
print(stiffness_SSE_list)
print(stiffness_SSE_err_list)
print(s_k_SSE_list)

plt.figure()
plt.scatter(N_list, energy_SSE_list[:,-1], label="SSE", color='C0')
plt.errorbar(N_list, energy_SSE_list[:,-1], energy_SSE_err_list[:,-1], fmt='None', capsize=5)
plt.plot(N_list, exact_energy, color='C3', label="Exact")
plt.legend()
plt.xlabel(r"$N$")
plt.ylabel("E/N")
plt.title(r"E/N vs N ($\Delta$={}, $\alpha$={})".format(Delta, alpha))
plt.tight_layout()
plt.savefig("Energy_per_site_vs_N_NN.png")

plt.figure()
plt.scatter(N_list, mag_z_SSE_list[:,-1], color='C0', label='SSE')
plt.plot(N_list, [0.0]*len(N_list), color='C3', label="Exact")
plt.errorbar(N_list, mag_z_SSE_list[:,-1], mag_z_SSE_err_list[:,-1], fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$N$")
plt.ylabel(r"$M_z$")
plt.title(r"E/N vs N ($\Delta$={}, $\alpha$={})".format(Delta, alpha))
plt.tight_layout()
plt.savefig("Magnetization_NN.png")

plt.figure()
plt.scatter(N_list, stiffness_SSE_list[:,-1], color='C0', label='SSE')
plt.plot(N_list, exact_stiffness, color='C3', label="Exact")
plt.errorbar(N_list, stiffness_SSE_list[:,-1], stiffness_SSE_err_list[:,-1], fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$N$")
plt.ylabel(r"$\rho_s$")
plt.title(r"Spin Stiffness ($\Delta$={}, $\alpha$={})".format(Delta, alpha))
plt.tight_layout()
plt.savefig("Spin_Stiffness_NN.png")