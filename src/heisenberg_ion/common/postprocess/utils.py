import math

import numpy as np


def combine_different_runs(data_1, data_2, drop_samples_1, drop_samples_2):

    data_3 = []
    data_1 = data_1[drop_samples_1:]
    data_2 = data_2[drop_samples_2:]

    for i in range(len(data_1)):
        data_3.append(data_1[i])

    for i in range(len(data_2)):
        data_3.append(data_2[i])

    return data_3


def statistics_binning(arr, auto_corr_drop, eq_drop):

    arr2 = arr[eq_drop:]
    arr3 = arr2[0::auto_corr_drop]
    workingNdim = int(math.log(len(arr3)) / math.log(2))
    mean = np.mean(arr3)
    standardError = max_error_binning(arr3, workingNdim - 6)
    return mean, standardError


def error_propagation(data):

    ndim = len(data)
    error = np.std(data, ddof=0) / np.sqrt(ndim)

    return error


def max_error_binning(data, dim):

    if dim <= 1:
        raise Exception("Not enough points MC steps were used for the binning method\n")

    error = np.zeros(dim)
    error[0] = error_propagation(data)

    for i in range(1, dim):
        bin_dim = int(len(data) / 2)
        data_bin = np.zeros(bin_dim)

        for j in range(bin_dim):
            data_bin[j] = 0.5 * (data[2 * j] + data[2 * j + 1])
        data = data_bin
        error[i] = error_propagation(data)

    return np.max(error)


def combine_dataset_stats(data_1, data_2):

    num_data_size_1 = len(data_1)
    stats_1 = statistics_binning(data_1)

    num_data_size_2 = len(data_2)
    stats_2 = statistics_binning(data_2)

    data_mean = np.mean([stats_1[0], stats_2[0]])
    numerator = np.sqrt((stats_1[1] * np.sqrt(num_data_size_1)) ** 2 + (stats_2[1] * np.sqrt(num_data_size_2)) ** 2)
    denominator = np.sqrt(num_data_size_1 + num_data_size_2)

    data_error = numerator / denominator

    return data_mean, data_error


def ed_energy(evals, T):

    beta = 1.0 / T

    energy = 0.0
    partition = 0.0

    for i in range(len(evals)):
        energy += np.exp(-beta * evals[i]) * evals[i]
        partition += np.exp(-beta * evals[i])

    energy /= partition

    return energy


def compute_histogram_from_shot_data(N, num_shots, shot_data):

    num_bits = 2**N
    freqs = np.zeros(num_bits)
    for i in range(num_shots):
        bit_str = 0
        for j in range(N):
            signed_spin_config = shot_data[i, j]
            bit_i_j = int(0.5 * (signed_spin_config + np.abs(signed_spin_config)))
            bit_str += (2 ** (j)) * (bit_i_j)
        freqs[bit_str] += 1.0 / num_shots

    return freqs


def kl_divergence(hist1, hist2, num_bits):

    kl_divergence = 0.0
    for j in range(num_bits):
        if hist1[j] != 0 and hist2[j] != 0:
            kl_divergence += hist1[j] * (np.log(hist1[j]) - np.log(hist2[j]))

    return kl_divergence
