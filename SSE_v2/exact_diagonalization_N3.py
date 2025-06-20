import numpy as np
import math
from matplotlib import pyplot as plt
import statistical_analysis as stats

# Simple construction for small scale testing
def construct_Hamiltonian_3(Delta, alpha, J=1.0, h=0.0):

    sigma_z = np.diag([1,-1])
    identity_matrix = np.eye(2)
    J_z = Delta * J * 0.25
    H_z = J_z * (np.kron(np.kron(sigma_z, sigma_z), identity_matrix) + np.kron(identity_matrix, np.kron(sigma_z, sigma_z)) + 
                ((0.5**alpha) * np.kron(np.kron(sigma_z, identity_matrix), sigma_z)))
    
    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0
    
    H_int = J*0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) + 
                   np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + 
                   np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) + 
                   np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix) + 
                   ((0.5**alpha) * (np.kron(np.kron(sigma_plus, identity_matrix), sigma_minus) + 
                                    np.kron(np.kron(sigma_minus, identity_matrix), sigma_plus))))
    
    j_x = 0.5j * (np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus)) - 
                   np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus)) + 
                   np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix) - 
                   np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix) + 
                   ((0.5**alpha) * (np.kron(np.kron(sigma_plus, identity_matrix), sigma_minus) - 
                                    np.kron(np.kron(sigma_minus, identity_matrix), sigma_plus))))
    
    return -H_z - H_int, H_int, j_x

# Simple construction for small scale testing
def construct_twisted_Hamiltonian_3(theta, Delta, alpha, J=1.0, h=0.0):

    sigma_z = np.diag([1,-1])
    identity_matrix = np.eye(2)
    J_z = Delta * J * 0.25
    H_z = J_z * (np.kron(np.kron(sigma_z, sigma_z), identity_matrix) + np.kron(identity_matrix, np.kron(sigma_z,sigma_z)) + 
                ((0.5**alpha) * np.kron(np.kron(sigma_z, identity_matrix), sigma_z)))
    
    sigma_plus = np.zeros((2,2))
    sigma_minus = np.zeros((2,2))
    sigma_plus[0,1] = 1.0
    sigma_minus[1,0] = 1.0
    
    H_int = J*0.5*(np.kron(identity_matrix, np.kron(sigma_plus, sigma_minus))*np.exp(theta*1j) + 
                   np.kron(identity_matrix, np.kron(sigma_minus, sigma_plus))*np.exp(-theta*1j) + 
                   np.kron(np.kron(sigma_plus, sigma_minus), identity_matrix)*np.exp(theta*1j) + 
                   np.kron(np.kron(sigma_minus, sigma_plus), identity_matrix)*np.exp(-theta*1j) + 
                   ((0.5**alpha) * (np.kron(np.kron(sigma_plus, identity_matrix), sigma_minus)*np.exp(theta*1j) + 
                                    np.kron(np.kron(sigma_minus, identity_matrix), sigma_plus)*np.exp(-theta*1j))))
    
    return -H_z - H_int

if __name__=="__main__":

    #Delta_list = [-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
    Delta_list = [0.0]
    h = 0.0
    alpha=1.0
    J = 1.0
    N = 3
    theta = 1.0
    T = 0.1
    beta = J/T
    spin_stiffness_list = []
    stiffness_SSE_list = []
    stiffness_SSE_err_list = []

    for Delta in Delta_list:

        if Delta < 0.0:
            gamma = 0.1
        else:
            gamma = 0.0

        '''
        Hamiltonian, T_x, j_x = construct_Hamiltonian_3(Delta, alpha, J, h)
        evals, evecs = np.linalg.eigh(Hamiltonian)
        exp_negative_T_x = 0.5*np.matmul(np.matmul(np.transpose(evecs[0,:]), (1.0/N)*T_x), evecs[0,:])
        perturbation_term = 0.0
        for i in range(len(evals)):
            if (np.abs(evals[i] - evals[0]) > 1e-15):
                j_x_overlap = np.matmul(np.matmul(np.transpose(evecs[0,:]), (1.0/N)*j_x), evecs[i,:])
                perturbation_term += np.absolute(j_x_overlap)**2/(evals[i] - evals[0])
        spin_stiffness = (3.0/2.0*N)*(exp_negative_T_x + perturbation_term)
        spin_stiffness_list.append(spin_stiffness)
        '''
        H_theta = construct_twisted_Hamiltonian_3(theta, Delta, alpha, J, h)
        H_zero, tx, jx = construct_Hamiltonian_3(Delta, alpha, J, h)
        H_minus_theta = construct_twisted_Hamiltonian_3(-theta, Delta, alpha, J, h)
        '''
        for row_index in range(8):
            print(H_theta[row_index,:])
        for row_index in range(8):
            print(H_zero[row_index,:])
        for row_index in range(8):
            print(H_minus_theta[row_index,:])
        '''
        evals_theta, evecs = np.linalg.eigh(H_theta)
        print(evals_theta)
        evals, evecs = np.linalg.eigh(H_zero)
        print(evals)
        evals_minus_theta, evecs = np.linalg.eigh(H_minus_theta)
        print(evals_minus_theta)
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
        second_derivative = (energy_theta + energy_minus_theta - 2.0*energy_zero)/(theta**2)
        #second_derivative = (evals_theta[0] - evals[0])/(theta**2)
        spin_stiffness = (3.0/(2.0 * N))*second_derivative
        spin_stiffness_list.append(spin_stiffness)

        file_1 = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_0.0_J_1.0_directed_loops_input_config_-1/MC Step Outputs.csv".format(N, Delta, h, alpha, gamma)
        step_number_1, energy_arr_1, magnetization_arr_1, stiffness_arr_1 = np.loadtxt(file_1, delimiter=",", skiprows=2, unpack=True)
        
        stiffness = stats.statistics_binning(stiffness_arr_1)

        stiffness_SSE_list.append(stiffness[0])
        stiffness_SSE_err_list.append(stiffness[1])


    print(spin_stiffness_list)
    plt.figure()
    plt.scatter(Delta_list, stiffness_SSE_list, label='SSE', color='C0')
    plt.errorbar(Delta_list, stiffness_SSE_list, stiffness_SSE_err_list, fmt="None", capsize=5, color='C0')
    plt.plot(Delta_list, spin_stiffness_list, label='ED', color='C3')
    plt.legend()
    plt.xlabel(r"$\Delta$")
    plt.ylabel(r"$\rho_s$")
    plt.title("Superfluid Density (N={})".format(N))
    plt.tight_layout()
    plt.savefig("Superfluid_Density_N_{}.png".format(N))
    plt.show()