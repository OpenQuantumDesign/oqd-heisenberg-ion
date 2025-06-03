import os
import numpy as np
import math
import matplotlib.pyplot as plt

def statistics_binning(arr):
    # Average and standard error using the binning method
    workingNdim  = int(math.log(len(arr))/math.log(2))
    trunc = int(len(arr)-2**workingNdim)
    mean = np.mean(arr[trunc:])
    standardError = max_error_binning(arr[trunc:], workingNdim-6)
    return mean, standardError

def error_propagation(data):
    ndim   = len(data)
    error = np.std(data,ddof=0)/np.sqrt(ndim)
    return error

def max_error_binning(data, workingNdim):
    if(workingNdim<=1):
        raise Exception('Not enough points MC steps were used for the binning method, please increase the number of MC steps')
    error = np.zeros(workingNdim)
    i = 0
    error[0] = error_propagation(data)

    for i in range(1,workingNdim):
        ndim = int(len(data)/2)
        data1 = np.zeros(ndim)

        for j in range(ndim):
            data1[j] = 0.5*(data[2*j]+data[2*j+1])
        data = data1
        error[i] = error_propagation(data)
    return np.max(error)

N=8
Delta = 1.1
h_list = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
alpha=1.0
J = 1.0
T = 0.5
beta = J/T
energy_ED_list = []
mag_z_ED_list = []
energy_SSE_list = []
energy_SSE_err_list = []
mag_z_SSE_list = []
mag_z_SSE_err_list = []

for h in h_list:

    file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/SSE/N_{}_Delta_{}_h_{}_alpha_{}_gamma_0.0_ksi_0.0_J_1.0_directed_loops/MC Step Outputs.csv".format(N, Delta, h, alpha)
    step_number, energy_arr, magnetization_arr = np.loadtxt(file, delimiter=",", skiprows=2, unpack=True)

    comparison_file = "/Users/shaeermoeed/Github/Heisenberg_Ion/Results/Exact_Diagonalization/ED_N_{}_Delta_{}_h_{}_Jx_{}_Jy_{}_alpha_{}_B_0.0_Heisenberg_OBC.csv".format(N, Delta, h, J, J, alpha)
    state_index, evals, mag_z, mag_x = np.loadtxt(comparison_file, skiprows=1, delimiter=",", unpack=True)

    ED_energy = 0.0
    ED_magnetization = 0.0
    partition_func = 0.0
    for i in range(len(evals)):
        ED_energy += evals[i]*np.exp(-beta * evals[i])
        ED_magnetization +=  mag_z[i] * np.exp(-beta * evals[i])
        partition_func += np.exp(-beta * evals[i])
    ED_energy /= partition_func
    ED_magnetization /= partition_func

    energy = statistics_binning(energy_arr)
    mag_z_exp = statistics_binning(magnetization_arr)
    print(energy)
    print(mag_z_exp)
    print(ED_energy)
    print(ED_magnetization)

    energy_ED_list.append(ED_energy)
    mag_z_ED_list.append(ED_magnetization)
    energy_SSE_list.append(energy[0])
    energy_SSE_err_list.append(energy[1])
    mag_z_SSE_list.append(mag_z_exp[0])
    mag_z_SSE_err_list.append(mag_z_exp[1])

plt.figure()
plt.plot(h_list, [0.0]*len(h_list), color="C3")
plt.scatter(h_list, np.array(energy_ED_list)/N-np.array(energy_SSE_list)/N, label="ED-SSE", color='C0')
plt.errorbar(h_list, np.array(energy_ED_list)/N - np.array(energy_SSE_list)/N, np.array(energy_SSE_err_list)/N, fmt='None', capsize=5)
plt.legend()
plt.xlabel("h")
plt.ylabel("E/N Difference")
plt.title("E/N Comparison")
plt.tight_layout()
plt.savefig("Energy per site comparison.png")

plt.figure()
plt.plot(h_list, [0.0]*len(h_list), color="C3")
plt.scatter(h_list, np.array(mag_z_ED_list)-np.array(mag_z_SSE_list), color='C0', label='ED-SSE')
plt.errorbar(h_list, np.array(mag_z_ED_list)-np.array(mag_z_SSE_list), np.array(mag_z_SSE_err_list), fmt='None', capsize=5)
plt.legend()
plt.xlabel("h")
plt.ylabel("M_z Difference")
plt.title("Magnetization Comparison")
plt.tight_layout()
plt.savefig("Magnetization_comparison.png")

