import numpy as np


def gamma_k(N, k, theta, alpha):
    """
    computes the coefficient gamma corresponding to the kth Fourier mode in linear spin wave theory for the power law decay spin-1/2 XY model
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        k (int): Fourier mode
        theta (float): twist angle
        alpha (float): power law decay exponent in interaction couplings

    Returns:
        (float): gamma_k
    """

    gamma = 0.0
    num_terms = int((N - 1) / 2)
    for r in range(1, num_terms + 1):
        gamma += (1.0 / (r**alpha)) * np.cos(2.0 * np.pi * r * k / N) * np.cos(r * theta)

    return gamma


def E_0_LSW(N, alpha, J, theta=0.0):
    """
    computes the ground state energy in linear spin wave theory for the power law decay spin-1/2 XY model.
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        alpha (float): power law decay exponent
        J (float): Energy scale
        theta (float, optional): Boundary twist angle. Defaults to 0.0.

    Returns:
        (float): ground state energy in linear spin wave theory
    """

    gamma_0 = gamma_k(N, 0, theta, alpha)
    energy = (-3.0 / 2.0) * N

    for j in range(N):
        gamma_j = gamma_k(N, j, theta, alpha)
        energy += np.sqrt(1.0 - gamma_j / gamma_0) + gamma_j / (2.0 * gamma_0)

    energy *= J * gamma_0 / 2.0

    return energy


def E_0_MF(N, alpha, J, theta=0.0):
    """
    computes the ground state energy in mean field theory for the power law decay spin-1/2 XY model.
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        alpha (float): power law decay exponent
        J (float): Energy scale
        theta (float, optional): Boundary twist angle. Defaults to 0.0.

    Returns:
        (float): ground state energy in mean field theory
    """

    gamma_0 = gamma_k(N, 0, theta, alpha)

    energy = -J * N * gamma_0 / 4.0

    return energy


def E_0_LSW_NN(N, J):
    """
    computes the ground state energy in linear spin wave theory for the nearest neighbor spin-1/2 XY model.
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        J (float): Energy scale

    Returns:
        (float): ground state energy in linear spin wave theory
    """

    energy = (-3.0 / 2.0) * N

    for j in range(N):
        gamma_j = np.cos(2.0 * np.pi * j / N)
        energy += np.sqrt(1.0 - gamma_j) + gamma_j / (2.0)

    energy *= J / 2.0

    return energy


def rho(N, alpha, J):

    rho_N_alpha = -E_0_LSW(N, alpha - 2.0, J) / N

    return rho_N_alpha


def rho_2(N, alpha, theta, J_1, J_2=1.0):
    """
    computes the spin stiffness for the power law decay spin-1/2 XY model in linear spin wave theory
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        alpha (float): power law decay exponent
        theta (float): twist angle
        J_1 (float): energy scale
        J_2 (float, optional): energy scale. Defaults to 1.0. Should typically be set equal to J_1

    Returns:
        (float): spin wave theory spin stiffness
    """

    energy_theta = E_0_LSW(N, alpha, J_1, theta)
    energy_0 = E_0_LSW(N, alpha, J_2, 0.0)
    rho_N_alpha = 2.0 * (energy_theta - energy_0) / (N * (theta**2))

    return rho_N_alpha


def rho_mf(N, alpha, theta, J):
    """
    computes the spin stiffness for the power law decay spin-1/2 XY model in mean field theory
    For details, see: https://arxiv.org/pdf/2601.20058

    Args:
        N (int): number of sites
        alpha (float): power law decay exponent
        theta (float): twist angle
        J (float): energy scale

    Returns:
        (float): mean field theory spin stiffness
    """

    energy_theta = E_0_MF(N, alpha, J, theta)
    energy_0 = E_0_MF(N, alpha, J, 0.0)
    rho_N_alpha = 2.0 * (energy_theta - energy_0) / (N * (theta**2))

    return rho_N_alpha
