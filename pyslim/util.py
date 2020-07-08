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
    w = label.astype("float64")
    n = np.bincount(1 + group, minlength=minlength + 1)
    x = np.bincount(1 + group, weights=w, minlength=minlength + 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        xm = x/n
    xm[n == 0] = 0
    w -= xm[1 + group]
    gw = np.bincount(1 + group, weights=np.abs(w), minlength=minlength + 1)[1:]
    # after subtracting groupwise means, should be all zero
    return np.logical_and(n[1:] > 0, np.abs(gw) < 1e-7)
