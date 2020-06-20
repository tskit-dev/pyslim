import numpy as np

def unique_labels_by_group(group, label, minlength=0):
    '''
    Given an array of integers ("group") from -1 and up and an array of integer
    "labels" of the same length, returns a logical value for each distinct
    value of ``group``, except for -1, indicating if all the entries of
    ``label`` that correspond to that value of ``group`` are identical. The
    return value will be of length at least ``minlength``.

    In other words, if the result is ``x``, then
    ``x[j]`` is ``len(set(label[group == j])) == 1``.
    '''
    n = np.bincount(1 + group, minlength=minlength + 1)[1:]
    x = np.bincount(1 + group, weights=label, minlength=minlength + 1)[1:]
    x2 = np.bincount(1 + group, weights=label.astype("int64") ** 2, minlength=minlength + 1)[1:]
    # (a * x)**2 = a * (a * x**2)
    return np.logical_and(n > 0, x**2 == n * x2)


