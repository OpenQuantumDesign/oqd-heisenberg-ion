import os
import numpy as np
import math
import matplotlib.pyplot as plt
import statistical_analysis as stats
import spin_wave_analysis as lsw

N=101
h = 0.0
alpha_list = []
#alpha_list.append(0.0)
alpha_list_2 = [0.0 + i*0.1 for i in range(51)]
alpha_list += alpha_list_2
alpha_list_3 = [1.0 + i*0.1 for i in range(41)]
#alpha_list.append(10.0)
print(alpha_list)
J = 1.0
T_1 = 0.005
T_2 = 0.05
T_3 = 5.0
#T_4 = 0.1
beta_1 = J/T_1
beta_2 = J/T_2
beta_3 = J/T_3
#beta_4 = J/T_4
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
num_bonds = int(N * (N-1)/2)
lsw_stiffness = []
theta = 0.01
pct_diff_list = []
lsw_energy = []
gamma_0_list = []

for i in range(len(alpha_list)):

    alpha = round(alpha_list[i], 1)

    '''
    gamma_phi = 0.0
    gamma_0 = 0.0
    theta = 0.000001
    for r in range(1, int((N+1)/2)):
        gamma_phi += np.cos(r*theta)/(r**(alpha))
        gamma_0 += 1.0/(r**(alpha))
    gamma_0_0 = (-J/2) * (gamma_phi - gamma_0)/(theta**2)
    '''
    
    norm = 0.0
    for r in range(1, int((N+1)/2)):
        norm += 1.0/(r**(alpha-2))

    if alpha >= 0.0 and alpha <= 5.0:
        folder = "SSE"
        #folder = "SSE_18_08_25"
    else: 
        folder = "SSE_18_08_25"

    init_start_config = 1
    start_config = init_start_config
    boundary=1
    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/{}/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(folder, N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T_1, loop_type, start_config, init_start_config)
    step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
    
    drop_1 = int(step_number_1[-1]/2)
    auto_corr_drop = 2

    print(step_number_1[-1])

    if hamiltonian_type == 0:
        energy_array = energy_arr_1
        stiffness_array = np.array(stiffness_arr_1)/norm
        #stiffness_array = np.array(stiffness_arr_1)
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

    if alpha >= 1.0:
        init_start_config = 1
        start_config = init_start_config
        boundary=1
        file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T_2, loop_type, start_config, init_start_config)
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
        
        drop_1 = int(step_number_1[-1]/2)
        auto_corr_drop = 2

        energy_array = energy_arr_1
        stiffness_array = stiffness_arr_1

        if hamiltonian_type == 0:
            energy_array = energy_arr_1
            stiffness_array = np.array(stiffness_arr_1)/norm
            #stiffness_array = np.array(stiffness_arr_1)
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
        boundary=1
        file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Step Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, boundary, T_3, loop_type, start_config, init_start_config)
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1, s_k_pi_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
        
        drop_1 = int(step_number_1[-1]/2)
        auto_corr_drop = 2

        energy_array = energy_arr_1
        stiffness_array = stiffness_arr_1

        if hamiltonian_type == 0:
            energy_array = energy_arr_1
            stiffness_array = np.array(stiffness_arr_1)/norm
            #stiffness_array = np.array(stiffness_arr_1)
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

    '''
    lsw_rho, lsw_E = lsw.rho(N, theta, alpha)

    pct_diff = 100*(lsw_rho - stiffness[0])/stiffness[0]
    pct_diff_list.append(pct_diff)

    lsw_stiffness.append(lsw_rho)
    lsw_energy.append(lsw_E)

    gamma_0_list.append(gamma_0_0)

    #print(gamma_0)
    '''

#print(np.array(energy_SSE_1_list)/N)
print(gamma_0_list)
print(np.array(stiffness_SSE_1_list))
#print(np.array(stiffness_SSE_1_err_list))
#print(s_k_SSE_list)
#print(lsw_stiffness)
#print(pct_diff_list)

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.scatter(alpha_list_3, np.array(energy_SSE_3_list)/N, label=r'$\beta={}$'.format(beta_3), color='black')
plt.errorbar(alpha_list_3, np.array(energy_SSE_3_list)/N, np.array(energy_SSE_3_err_list)/N, fmt='None', capsize=5, color='black')
plt.scatter(alpha_list_3, np.array(energy_SSE_2_list)/N, label=r'$\beta={}$'.format(int(beta_2)), color='C3')
plt.errorbar(alpha_list_3, np.array(energy_SSE_2_list)/N, np.array(energy_SSE_2_err_list)/N, fmt='None', capsize=5, color='C3')
plt.scatter(alpha_list, np.array(energy_SSE_1_list)/N, label=r'$\beta={}$'.format(int(beta_1)), color='C0')
plt.errorbar(alpha_list, np.array(energy_SSE_1_list)/N, np.array(energy_SSE_1_err_list)/N, fmt='None', capsize=5, color='C0')
plt.plot(alpha_list, [-1.0/np.pi] * len(alpha_list), color='green', linestyle='dashed')
plt.legend(prop={'size': 12})

#plt.plot(alpha_list, np.array(lsw_energy)/N, label='SSE', color='C3')
plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$E/N$", fontsize=16)
#plt.title("E per site Comparison vs alpha")
#plt.tight_layout()
plt.savefig("Energy_per_site.png", bbox_inches='tight', dpi=1200)

plt.figure()
plt.scatter(alpha_list_3, np.array(mag_z_SSE_3_list), label=r'$\beta={}$'.format(beta_3), color='black')
plt.errorbar(alpha_list_3, np.array(mag_z_SSE_3_list), np.array(mag_z_SSE_3_err_list), fmt='None', capsize=5, color='black')
plt.scatter(alpha_list_3, np.array(mag_z_SSE_2_list), label=r'$\beta={}$'.format(int(beta_2)), color='C3')
plt.errorbar(alpha_list_3, np.array(mag_z_SSE_2_list), np.array(mag_z_SSE_2_err_list), fmt='None', capsize=5, color='C3')
plt.scatter(alpha_list, np.array(mag_z_SSE_1_list), label=r'$\beta={}$'.format(int(beta_1)), color='C0')
plt.errorbar(alpha_list, np.array(mag_z_SSE_1_list), np.array(mag_z_SSE_1_err_list), fmt='None', capsize=5, color='C0')
plt.plot(alpha_list, [0.0] * len(alpha_list), color='green', linestyle='dashed')
plt.legend(prop={'size': 12})

#plt.plot(alpha_list, np.array(lsw_energy)/N, label='SSE', color='C3')
plt.legend()
plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$M_z$", fontsize=16)
#plt.title("E per site Comparison vs alpha")
#plt.tight_layout()
plt.savefig("Magnetization.png", bbox_inches='tight', dpi=1200)

fig, ax1 = plt.subplots()
#plt.plot(alpha_list, gamma_0_list, color='black', label='MF')
ax1.scatter(alpha_list_3, np.array(stiffness_SSE_3_list), color='black', label=r'$\beta={}$'.format(beta_3))
ax1.errorbar(alpha_list_3, np.array(stiffness_SSE_3_list), np.array(stiffness_SSE_3_err_list), fmt='None', capsize=5, color='black')
ax1.scatter(alpha_list_3, np.array(stiffness_SSE_2_list), color='C3', label=r'$\beta={}$'.format(int(beta_2)))
ax1.errorbar(alpha_list_3, np.array(stiffness_SSE_2_list), np.array(stiffness_SSE_2_err_list), fmt='None', capsize=5, color='C3')
ax1.scatter(alpha_list, np.array(stiffness_SSE_1_list), color='C0', label=r'$\beta={}$'.format(int(beta_1)))
ax1.errorbar(alpha_list, np.array(stiffness_SSE_1_list), np.array(stiffness_SSE_1_err_list), fmt='None', capsize=5)
ax1.plot(alpha_list, [1/np.pi] * len(alpha_list), color='green', linestyle='dashed')
#plt.plot(alpha_list, lsw_stiffness, color='C3', label='LSW')
ax1.legend(prop={'size': 12})
drop_from_top = 14
#left,bottom,width,height = [0.39, 0.28, 0.49, 0.37]
left,bottom,width,height = [1.1, 0.0, 0.8, 0.8]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.scatter(alpha_list[drop_from_top+10:], np.array(stiffness_SSE_3_list)[drop_from_top:], color='black', label=r'SSE ($\beta={}$)'.format(beta_3), s=25)
ax2.errorbar(alpha_list[drop_from_top+10:], np.array(stiffness_SSE_3_list)[drop_from_top:], np.array(stiffness_SSE_3_err_list)[drop_from_top:], fmt='None', capsize=5, color='black')
ax2.scatter(alpha_list[drop_from_top+10:], np.array(stiffness_SSE_2_list)[drop_from_top:], color='C3', label=r'SSE ($\beta={}$)'.format(int(beta_2)), s=25)
ax2.errorbar(alpha_list[drop_from_top+10:], np.array(stiffness_SSE_2_list)[drop_from_top:], np.array(stiffness_SSE_2_err_list)[drop_from_top:], fmt='None', capsize=5, color='C3')
ax2.scatter(alpha_list[drop_from_top:], np.array(stiffness_SSE_1_list)[drop_from_top:], color='C0', label=r'SSE ($\beta={}$)'.format(int(beta_1)), s=25)
ax2.errorbar(alpha_list[drop_from_top:], np.array(stiffness_SSE_1_list)[drop_from_top:], np.array(stiffness_SSE_1_err_list)[drop_from_top:], fmt='None', capsize=5)
ax2.plot(alpha_list[drop_from_top:], [1/np.pi] * len(alpha_list[drop_from_top:]), color='green', linestyle='dashed')
ax2.locator_params(axis='x', nbins=8)

ax1.set_xlabel(r"$\alpha$", fontsize=16)
ax1.set_ylabel(r"$\rho_s$", fontsize=16)

ax2.set_xlabel(r"$\alpha$", fontsize=14)
ax2.set_ylabel(r"$\rho_s$", fontsize=14)

plt.savefig("Superfluid_Density.png", bbox_inches='tight', dpi=1200)

'''
plt.figure()
plt.scatter(alpha_list, np.array(s_k_SSE_list), color='C0', label='SSE')
plt.errorbar(alpha_list, np.array(s_k_SSE_list), np.array(s_k_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\rho_s$")
plt.title(r"S($\pi$) ($\Delta$={})".format(Delta))
plt.tight_layout()
plt.savefig("S_k_pi.png")
'''
