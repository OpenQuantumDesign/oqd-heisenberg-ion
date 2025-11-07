import numpy as np
import matplotlib.pyplot as plt

def get_energy(N, phi):

    s_min = -int(N/2)
    E = 0.0
    for i in range(N):
        s = s_min + i
        eval = -np.cos(2 * np.pi * s/N - phi/N)
        if eval < 0.0:
            E += eval
    
    '''
    if (s != s_max):
        raise Exception("s not equal to s_max")
    '''
    return E/N

N_list = [5,7,9,11,13,102] # Should be even but not divisible by 4
phi = 0.01
rho_list = []
for j in range(len(N_list)):

    N = N_list[j]

    E_phi = get_energy(N, phi)
    E_zero = get_energy(N, 0.0)
    E_minus_phi = get_energy(N, -phi)

    #E_phi = -np.cos(phi/N)/(N*np.sin(np.pi/N))
    #E_zero = -np.cos(0.0/N)/(N*np.sin(np.pi/N))
    #E_minus_phi = -np.cos(-phi/N)/(N*np.sin(np.pi/N))

    rho =(N**2) * ((E_phi + E_minus_phi - 2.0*E_zero)/((phi**2)))
    rho_list.append(rho)


print(rho)
print(1.0/np.pi)
plt.figure()
plt.scatter(N_list, rho_list, color='C0')
plt.plot(N_list, [1.0/np.pi] * len(N_list), color='C3')
plt.show()    