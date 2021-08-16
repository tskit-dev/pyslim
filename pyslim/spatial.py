import msprime
import tskit
import warnings
import numpy as np

from .slim_tree_sequence import *


def _in_location_bin(locations, x0, x1, y0, y1):
    '''
    Takes the locations of individuals and tests if each is within [x0, x1) and
    [y0, y1). Returns a list the same length as the number of individuals with
    True if the corresponding individual is within the bounds and False
    otherwise.  This is a helper function for `population_size`, and may change
    or be removed in the future.

    :param array locations: An n x 3 array with rows individuals, the first
        column the x coordinates, and the second column the y coordinates.
    :param float x0: The lower x coordinate boundary.
    :param float x1: The upper x coordinate boundary.
    :param float y0: The lower y coordinate boundary.
    :param float y1: The upper y coordinate boundary.
    '''
    return np.logical_and(
            np.logical_and(
                locations[:,0] < x1,
                locations[:,0] >= x0),
            np.logical_and(
                locations[:,1] < y1,
                locations[:,1] >= y0
            )
        )


def _average_time_alive(birth_times, death_times, t0, t1):
    '''
    Calculates the average population size in the time interval [`t0`, `t1`)
    from birth times and death times of individuals. This is a helper function for
    `population_size`, and may change or be removed in the future.

    'Average population size' is the sum over all individuals of the length of
    time each individual was alive in [`t0`, `t1`), divided by `t1 - t0`. For
    integer t0 and t1, this is equivalent to the average over all integer time
    steps in [t0, t1) of the number of individuals alive at that time step;
    however, this is not restricted to integer times.

    Time goes bacwards, so t1 is further in the past than t0 and birth times
    are greater than death times. Birth times and death times should already
    include an offset for the stage that individuals are remembered at. See
    `population_size` code.

    For example, suppose t0 = 1, t1 = 4, and there are four individuals with
    birth times [3, 1, 5, 4] and death times [0, 1, 2, 4]. The first individual
    is alive for three time steps in [t0, t1): 1, 2, and 3; the second
    individual for one time step: 1; the third individual for two time steps: 2
    and 3, and the fourth individual for zero (since the time interval is open
    on the right and it was born at t1). So the average population size is (3 +
    1 + 2 + 0)/(4-1) = 2.

    :param array birth_times: Birth times of individuals.
    :param array death_times: Death times of individuals. Should be the same
        length as `birth_times`. Death times are less than or equal to birth times.
    :param float t0: Lower time endpoint.
    :param float t1: Upper time endpoint. t1 is greater than t0.

    '''
    # length of the segment of time in [t0, t1) that the individual
    # was alive for, i.e., of the intersection with [birth, death].
    a = np.maximum(0, np.minimum(t1, 1 + birth_times) -  np.maximum(t0, death_times))
    return sum(a)/(t1-t0)


def population_size(ts, x_bins, y_bins, time_bins, stage='late', remembered_stage=None):
    '''
    Calculates the population size in each of the spatial bins defined by grid lines
    at ``x_bins`` and ``y_bins``, averaged over each of the time intervals separated by
    ``time_bins``. To obtain actual (census) sizes, the tree sequence must contain
    all individuals alive, e.g., from a SLiM simulation with all individuals
    permanently remembered.

    With ``nx``, ``ny`` and ``nt`` the number of bins in the ``x``, ``y`` and time directions
    (so, ``nx`` is one less than the length of ``x_bins``),
    returns a 3-d array with dimensions ``(nx, ny, nt)``.
    The ``(i,j,k)``th element of the array is the average number of
    individuals with ``x`` coordinate in the half-open interval ``[x_bins[i], x_bins[i + 1])``
    and ``y`` coordinate in the half-open interval ``[y_bins[i], y_bins[i + 1])``,
    averaged across all times in the half-open interval ``[time_bins[k], time_bins[k + 1])``.

    For integer endpoints of the time bins, this average is equivalent to
    recording the number of indivduals that are alive at each time in the time
    interval and have location in the relevant location bin, then taking the
    mean of these recorded population sizes.

    :param TreeSequence ts: The tree sequence to calculate population size from.
    :param array x_bins: The x-coordinates of the boundaries of the location bins.
    :param array y_bins: The y-coordinates of the boundaries of the location bins.
    :param array time_bins: The endpoints of the time bins.
    :param str stage: The stage in the SLiM life cycle that the endpoints of
        the time bins refer to (either "early" or "late"; defaults to "late").
    :param str remembered_stage: The stage in the SLiM life cycle during which
        individuals were Remembered (defaults to the stage the tree sequence was
        recorded at, stored in metadata).
    '''

    # TODO: make a helper function to do this
    if stage not in ("late", "early"):
        raise ValueError(f"Unknown stage '{stage}': "
                          "should be either 'early' or 'late'.")

    if remembered_stage is None:
        remembered_stage = ts.metadata['SLiM']['stage']

    if remembered_stage not in ("late", "early"):
        raise ValueError(f"Unknown remembered_stage '{remembered_stage}': "
                          "should be either 'early' or 'late'.")
    if remembered_stage != ts.metadata['SLiM']['stage']:
        warnings.warn(f"Provided remembered_stage '{remembered_stage}' does not"
                      " match the stage at which the tree sequence was saved"
                      f" ('{ts.metadata['SLiM']['stage']}'). This is not necessarily"
                      " an error, but mismatched stages will lead to inconsistencies:"
                      " make sure you know what you're doing.")

    # see individuals_alive_at for explanation
    if stage == "early" and ts.metadata['SLiM']['model_type'] == "WF":
        birth_offset = 1
    else:
        birth_offset = 0
    birth_times = ts.individual_times - birth_offset
    if (ts.metadata['SLiM']['model_type'] == "WF"
            or stage == remembered_stage):
        age_offset = 0
    else:
        if (remembered_stage == "early"
                and stage == "late"):
            age_offset = -1
        else:
            age_offset = 1
    ages = ts.individual_ages + age_offset
    death_times = birth_times - ages

    time_breaks = time_bins
    x_breaks = x_bins
    y_breaks = y_bins

    nxbins = len(x_breaks) - 1
    nybins = len(y_breaks) - 1
    ntbins = len(time_breaks) - 1
    popsize = np.empty((nxbins, nybins, ntbins))

    locations = ts.individual_locations

    for i in np.arange(nxbins):
        for j in np.arange(nybins):
            x0, x1 = x_breaks[i], x_breaks[i + 1]
            y0, y1 = y_breaks[j], y_breaks[j + 1]
            in_bin = _in_location_bin(locations, x0, x1, y0, y1)
            for k in np.arange(ntbins):
                t0, t1 = time_breaks[k], time_breaks[k + 1]
                popsize[i, j, k] = _average_time_alive(
                        birth_times[in_bin],
                        death_times[in_bin],
                        t0,
                        t1,
                )

    return popsize
