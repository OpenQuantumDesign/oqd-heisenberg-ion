import numpy as np
import matplotlib.pyplot as plt
import statistical_analysis as stats

def structure_factor(k, pair_correlations, N, auto_corr_drop, eq_drop):

    array_shape = np.shape(pair_correlations[-1000:,:])
    num_steps = array_shape[0]
    print(num_steps)
    S_k_array = np.zeros(num_steps, dtype=np.complex64)
    b=0
    for i in range(num_steps):
        b = 0
        for i_b in range(N):
            for j_b in range(i_b+1,N):
                S_k_array[i] = np.exp(1j*k*(i_b - j_b)) * pair_correlations[i,b]
                b += 1
    

    S_k = stats.statistics_binning(S_k_array, auto_corr_drop, eq_drop,  np.complex64)

    return S_k

N=10
h=0.0
alpha_list = [0.0]
J = 1.0
T = 0.1
beta = J/T
energy_ED_list = []
mag_z_ED_list = []
stiffness_ED_list = []
energy_SSE_list = []
energy_SSE_err_list = []
mag_z_SSE_list = []
mag_z_SSE_err_list = []
stiffness_SSE_list = []
stiffness_SSE_err_list = []
theta = 0.1
start_config = 1
loop_type = "deterministic"
dist_dep_offset = 0
gamma = 0.0
Delta = 0.0
hamiltonian_type = int(Delta)
eq_drop = 0
auto_corr_drop = 1

S_k = np.zeros((len(alpha_list), 3))

for i in range(len(alpha_list)):

    alpha = alpha_list[i]

    init_start_config = 0
    start_config = init_start_config
    file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_dist_dep_offset_{}_{}_input_config_{}_initial_input_config_{}/Pair Correlation Outputs.csv".format(N, hamiltonian_type, Delta, h, alpha, gamma, dist_dep_offset, loop_type, start_config, init_start_config)
    data = np.loadtxt(file_1, delimiter=",", skiprows=49000, unpack=False)

    pair_correlations = data[:,1:]
    num_bonds = int(N*(N-1)/2)
    print(num_bonds)
    pair_correlation_stats = np.zeros((2,num_bonds))
    #S_k_SSE = structure_factor(np.pi, pair_correlations, N, auto_corr_drop, eq_drop)
    #k_list = np.linspace(-np.pi,np.pi)
    k_list = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
    S_k_list = np.zeros((2,len(k_list)), dtype=np.complex64)

    '''
    S_k[i,0] = alpha
    S_k[i,1] = S_k_SSE[0]
    S_k[i,2] = S_k_SSE[1]
    '''
    for i in range(len(k_list)):
        k = k_list[i]
        S_k_list[:,i] = structure_factor(k, pair_correlations, N, auto_corr_drop, eq_drop)
    
plt.plot(k_list,S_k_list)
plt.show()


np.savetxt("Structure_Factor_pi.csv", S_k, delimiter=",", header='alpha,S_k_mean,S_k_err')

S_k_data = np.loadtxt("Structure_Factor_pi.csv", delimiter=",", skiprows=1, unpack=False)

plt.plot(alpha_list, S_k_data[:,1], color='C0')
plt.errorbar(alpha_list, S_k_data[:,1], S_k_data[:,2], fmt='None', capsize=5)
plt.ylabel(r"$S(\pi)$")
plt.xlabel(r"$\alpha$")
plt.title(r"S($\pi$) vs $\alpha$")
plt.savefig("structure_factor.png")



