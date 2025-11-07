import numpy as np

def update_diagonal_estimators(spin_array, magnetization, step, N):

    M_z = 0.0
    for i in range(N):
        M_z += spin_array[i]
        
    M_z *= 0.5/N
    magnetization[step] = M_z

    return magnetization