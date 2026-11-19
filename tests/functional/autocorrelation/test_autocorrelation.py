import math
import os

import numpy as np
from matplotlib import pyplot as plt

from oqd_heisenberg_ion.common.driver.factory import DriverFactory
from oqd_heisenberg_ion.common.inputs.input_reader import InputReader
from oqd_heisenberg_ion.common.preprocessor.factory import PreprocessorFactory


def compare_runs(folder_old, folder_new):

    exclude_keys = {"root_folder", "output_folder_name", "simulation_folder", "run_folder", "uuid"}

    kvp_settings1 = read_input_file(folder_old)
    kvp_settings2 = read_input_file(folder_new)

    for i in range(len(kvp_settings1.keys())):
        key = list(kvp_settings1.keys())[i]
        if key not in exclude_keys and kvp_settings1[key] != kvp_settings2[key]:
            raise Exception("Values for key {} not equal".format(key))

    num_parameter_sets = len((kvp_settings1['simulator'].strip().split(",")))

    alpha_list = kvp_settings1['alpha'].strip().split(",")

    energy1_list = []
    energy1_err_list = []

    stiffness1_list = []
    stiffness1_err_list = []

    energy2_list = []
    energy2_err_list = []

    stiffness2_list = []
    stiffness2_err_list = []

    for i in range(num_parameter_sets):

        alpha = alpha_list[i]

        uuid1 = kvp_settings1['uuid'].strip().split(",")[i]
        estimators1 = extract_estimators(folder_old, uuid1)

        uuid2 = kvp_settings2['uuid'].strip().split(",")[i]
        estimators2 = extract_estimators(folder_new, uuid2)

        energy_stats1 = statistics_binning(estimators1[:,1], 1, 0)
        energy_stats2 = statistics_binning(estimators2[:,1], 1, 0)
        stiffness_stats1 = statistics_binning(estimators1[:,3], 1, 0)
        stiffness_stats2 = statistics_binning(estimators2[:,3], 1, 0)

        plt.figure()
        plt.scatter(np.arange(len(energy_stats1[1])), (energy_stats1[1]/energy_stats1[1][0])**2, label='old')
        plt.scatter(np.arange(len(energy_stats2[1])), (energy_stats2[1]/energy_stats2[1][0])**2, label='new')
        plt.xlabel("bin level")
        plt.ylabel("Error")
        plt.legend()
        plt.title("Energy Error (alpha = {})".format(alpha))

        plt.figure()
        plt.scatter(np.arange(len(stiffness_stats1[1])), (stiffness_stats1[1]/stiffness_stats1[1][0])**2, label='old')
        plt.scatter(np.arange(len(stiffness_stats2[1])), (stiffness_stats2[1]/stiffness_stats2[1][0])**2, label='new')
        plt.xlabel("bin level")
        plt.ylabel("Error")
        plt.legend()
        plt.title("Stiffness Error (alpha = {})".format(alpha))

        energy1_list.append(energy_stats1[0])
        energy2_list.append(energy_stats2[0])
        stiffness1_list.append(stiffness_stats1[0])
        stiffness2_list.append(stiffness_stats2[0])

        energy1_err_list.append(np.max(energy_stats1[1]))
        energy2_err_list.append(np.max(energy_stats2[1]))
        stiffness1_err_list.append(np.max(stiffness_stats1[1]))
        stiffness2_err_list.append(np.max(stiffness_stats2[1]))

    plt.figure()
    plt.scatter(alpha_list, energy1_list, label='old')
    plt.errorbar(alpha_list, energy1_list, energy1_err_list, fmt='None', capsize=5)
    plt.scatter(alpha_list, energy2_list, label='new')
    plt.errorbar(alpha_list, energy2_list, energy2_err_list, fmt='None', capsize=5)
    plt.xlabel('alpha')
    plt.ylabel('Energy')
    plt.legend()

    plt.figure()
    plt.scatter(alpha_list, stiffness1_list, label='old')
    plt.errorbar(alpha_list, stiffness1_list, stiffness1_err_list, fmt='None', capsize=5)
    plt.scatter(alpha_list, stiffness2_list, label='new')
    plt.errorbar(alpha_list, stiffness2_list, stiffness2_err_list, fmt='None', capsize=5)
    plt.xlabel('alpha')
    plt.ylabel('Stiffness')
    plt.legend()



def read_input_file(folder):

    input_file = os.path.join(folder, "inputs.txt")
    data_dict = {}
    with open(input_file, 'r') as f:
            for line in f.readlines():
                words = line.strip().split("\t")
                key = words[0]
                vals = words[1]
                data_dict[key] = vals

    return data_dict


def extract_estimators(folder, uuid):

    uuid_folder = os.path.join(folder, uuid)
    out_folder = os.path.join(uuid_folder, 'qmc_output')
    estimators = np.loadtxt(os.path.join(out_folder, 'estimators.csv'), skiprows=2, delimiter=',')

    return estimators


def statistics_binning(arr, auto_corr_drop, eq_drop):

    arr2 = arr[eq_drop:]
    arr3 = arr2[0::auto_corr_drop]
    workingNdim = int(math.log(len(arr3)) / math.log(2))
    mean = np.mean(arr3)
    standardErrors = max_error_binning(arr3, workingNdim-6)
    return mean, standardErrors


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

    return error

input_file = "tests/functional/autocorrelation/input_files/long_range_XY.txt"
inputs = InputReader(input_file_path=input_file)
inputs.read_inputs_from_file()

parameter_set_list = inputs.parameter_set_list
preprocessor_qmc = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
driver_inputs = preprocessor_qmc.preprocess()
driver = DriverFactory.create("long_range_qmc", preprocessor_qmc.simulation_folder, driver_inputs)
driver.simulate()

folder_old = "tests/functional/autocorrelation/results/benchmark"
folder_new = preprocessor_qmc.simulation_folder
compare_runs(folder_old, folder_new)
plt.show()