import numpy as np
import statistical_analysis as stats
import matplotlib.pyplot as plt

def compute_histogram_from_shot_data(N, num_shots, num_bits, shot_data):

    freqs = np.zeros(num_bits)
    for i in range(num_shots):
        bit_str = 0
        for j in range(N):
            signed_spin_config = shot_data[i,j]
            bit_i_j = int(0.5*(signed_spin_config + np.abs(signed_spin_config)))
            bit_str += (2**(j)) * (bit_i_j)
        freqs[bit_str] += 1.0/num_shots

    return freqs

def extract_shot_data(alpha, N, T):

    Delta = 0.0
    hamiltonian_type = int(Delta)
    h = 0.0
    gamma = 0.0
    loop_type = "deterministic"
    dist_dep_offset = 0
    boundary=0
    init_start_config = 1
    start_config = init_start_config

    if alpha[0:3] == "exp":
        J = "exp"
    else:
        J = 1.0

    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_{}_dist_dep_offset_{}_boundary_{}_T_{}_{}_input_config_{}_initial_input_config_{}/MC Spin Configurations.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, J, dist_dep_offset, boundary, T, loop_type, start_config, init_start_config)
    shot_data = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=False)

    return shot_data

N = 11
T = 0.005
beta = 1.0/T
num_bits = 2**N
'''
alpha_list = ["exp_mu_1.1", "exp_mu_1.2", "exp_mu_1.3", "exp_mu_1.4", "exp_mu_1.5", "exp_mu_1.6", "exp_mu_1.7", "exp_mu_1.8", "exp_mu_1.9", "exp_mu_2.0"]

shot_data = extract_shot_data("exp_mu_1.41", N, T)
num_steps = np.shape(shot_data)[0]

hist_freqs_exp_J_ij = compute_histogram_from_shot_data(N, num_steps, num_bits, shot_data)

plt.figure()
plt.stem(np.linspace(0,num_bits-1, num_bits), hist_freqs_exp_J_ij)
plt.savefig("Exp_J_ij_hist.png")

kl_divergence_array = np.zeros((len(alpha_list), 2))

for i in range(len(alpha_list)):

    alpha = alpha_list[i]

    shot_data = extract_shot_data(alpha, N, T)
    num_steps = np.shape(shot_data)[0]

    hist_freqs = compute_histogram_from_shot_data(N, num_steps, num_bits, shot_data)

    kl_divergence = 0.0
    for j in range(num_bits):
        if hist_freqs[j] != 0 and hist_freqs_exp_J_ij[j] != 0:
            kl_divergence += hist_freqs[j] * (np.log(hist_freqs[j]) - np.log(hist_freqs_exp_J_ij[j]))

    kl_divergence_array[i,0] = float(alpha[-3:])
    kl_divergence_array[i,1] = kl_divergence

np.savetxt("kl_divergence.csv", kl_divergence_array, header='alpha,KL\n')
'''
kl_data = np.loadtxt("kl_divergence.csv", skiprows=2)
plt.figure()
plt.scatter(kl_data[1:,0], kl_data[1:,1])
plt.xlabel(r"$\lambda$")
plt.ylabel(r"$D_{KL}$")
plt.savefig("kl_divergence.png")

    



    