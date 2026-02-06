import numpy as np


def set_probability(val):

    if val < 0.0:
        if np.abs(val) < 1e-15:
            print("Correcting small negative")
            return 0.0
        else:
            print("invalid probability encountered: ", val)
    else:
        return val


def enforce_positive(table):

    size = np.shape(table)[0]
    for i in range(size):
        table[i] = set_probability(table[i])

    return 0
