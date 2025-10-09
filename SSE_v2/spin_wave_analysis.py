import numpy as np
import matplotlib.pyplot as plt

def gamma_k(N, k, theta, alpha):

    gamma = 0.0
    num_terms = int((N-1)/2)
    for r in range(1,num_terms+1):
        gamma += (1.0/(r**alpha)) * np.cos(2.0 * np.pi * r * k/N) * np.cos(r * theta)
    
    return gamma

def E_0_LSW(N, alpha, J, theta=0.0):

    gamma_0 = gamma_k(N, 0, theta, alpha)
    energy = (-3.0/2.0) * N

    for j in range(N):
        gamma_j = gamma_k(N, j, theta, alpha)
        energy += np.sqrt(1.0 - gamma_j/gamma_0) + gamma_j/(2.0 * gamma_0)

    energy *= J*gamma_0/2.0

    return energy

def E_0_MF(N, alpha, J, theta=0.0):

    gamma_0 = gamma_k(N, 0, theta, alpha)

    energy = -J * N * gamma_0/4.0

    return energy

def E_0_LSW_NN(N, J):

    energy = (-3.0/2.0) * N

    for j in range(N):
        gamma_j = np.cos(2.0 * np.pi * j/N)
        energy += np.sqrt(1.0 - gamma_j) + gamma_j/(2.0)

    energy *= J/2.0

    return energy

def rho(N, alpha, J):

    rho_N_alpha = -E_0_LSW(N, alpha-2.0, J)/N

    return rho_N_alpha

def rho_2(N, alpha, theta, J_1, J_2=1.0):

    energy_theta = E_0_LSW(N, alpha, J_1, theta)
    energy_0 = E_0_LSW(N, alpha, J_2, 0.0)
    rho_N_alpha = 2.0 * (energy_theta - energy_0)/(N*(theta**2))

    return rho_N_alpha

def rho_mf(N, alpha, theta, J):

    energy_theta = E_0_MF(N, alpha, J, theta)
    energy_0 = E_0_MF(N, alpha, J, 0.0)
    rho_N_alpha = 2.0 * (energy_theta - energy_0)/(N*(theta**2))

    return rho_N_alpha



if __name__ == "__main__":

    N = 101
    phi = 0.01
    alpha_list = [1.0,2.0,3.0,4.0,5.0,10.0,20.0,50.0]
    rho_list = []
    for i in range(len(alpha_list)):
        alpha = alpha_list[i]
        rho_alpha = rho(N, phi, alpha)
        rho_list.append(rho_alpha[0])

    print(rho_list)

    plt.plot(alpha_list, rho_list)
    plt.show()