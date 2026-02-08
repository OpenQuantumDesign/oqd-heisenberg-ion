import numpy as np


def free_fermion_energy(N, phi):

    # N should not be divisible by 4

    s_min = -int(N / 2)
    E = 0.0
    for i in range(N):
        s = s_min + i
        eval = -np.cos(2 * np.pi * s / N - phi / N)
        if eval < 0.0:
            E += eval

    """
    if (s != s_max):
        raise Exception("s not equal to s_max")
    """
    return E / N


def free_fermion_stiffness(N, phi):

    E_phi = free_fermion_energy(N, phi)
    E_zero = free_fermion_energy(N, 0.0)
    E_minus_phi = free_fermion_energy(N, -phi)

    rho = (N**2) * ((E_phi + E_minus_phi - 2.0 * E_zero) / (phi**2))

    return rho
