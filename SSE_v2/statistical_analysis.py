import numpy as np
import math

def statistics_binning(arr, auto_corr_drop, eq_drop):
    # Average and standard error using the binning method
    arr2 = arr[eq_drop:]
    arr3 = arr2[0::auto_corr_drop]
    workingNdim  = int(math.log(len(arr3))/math.log(2))
    trunc = int(len(arr3)-2**workingNdim)
    mean = np.mean(arr3)
    standardError = max_error_binning(arr3, workingNdim-6)
    return mean, standardError

def error_propagation(data):
    ndim = len(data)
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

def combine_dataset_stats(data_1, data_2):

    num_data_size_1 = len(data_1)
    stats_1 = statistics_binning(data_1)

    num_data_size_2 = len(data_2)
    stats_2 = statistics_binning(data_2)

    data_mean = np.mean([stats_1[0], stats_2[0]])
    data_error = np.sqrt((stats_1[1]*np.sqrt(num_data_size_1))**2 + (stats_2[1]*np.sqrt(num_data_size_2))**2)/(np.sqrt(num_data_size_1 + num_data_size_2))

    return data_mean, data_error
    