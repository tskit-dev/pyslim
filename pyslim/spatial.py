import msprime
import tskit
import warnings
import numpy as np

# probably don't need all of these
from .slim_tree_sequence import *
from .slim_metadata import *
from .provenance import *
from .util import *


def in_location_bin(locations, x0, x1, y0, y1):
    '''
    Takes the locations of individuals and tests if each is within [x0, x1) and [y0, y1). Returns a list the same length as the number of individuals with True if the corresponding individual is within the bounds and False otherwise.
    
    :param array locations: An n x 3 array with rows individuals, the first column the x coordinates, and the second column the y coordinates.
    :param float x0: The lower x coordinate boundary.
    :param float x1: The upper x coordinate boundary.
    :param float y0: The lower y coordinate boundary.
    :param float y1: The upper y coordinate boundary.
    '''
    return(np.logical_and(np.logical_and(locations[:,0] < x1, locations[:,0] >= x0),
               np.logical_and(locations[:,1] < y1, locations[:,1] >= y0)))

def average_time_alive(birth_times, death_times, t0, t1):
    '''
    Calculates the average population size in [`t0`, `t1`) from birth times and death times of individuals. Helper function for `population_size`. 

    'Average population size' is the sum over all individuals of the length of time each individual was alive in [`t0`, `t1`), divided by `t1 - t0`. For integer t0 and t1, this is equivalent to the average over all integer time steps in [t0, t1) of the number of individuals alive at that time step.

    Time goes bacwards, so t1 is further in the past than t0 and birth times are greater than death times. Birth times and death times should already include an offset for the stage that individuals are remembered at. See `population_size` code.

    For example, suppose t0 = 1, t1 = 4, and there are four individuals with birth times [3, 1, 5, 4] and death times [0, 1, 2, 4]. The first individual is alive for three time steps in [t0, t1): 1, 2, and 3; the second individual for one tome step: 1; the third individual for two time steps: 2 and 3, and the fourth individual for zero (since the time interval is open on the right and it was born at t1). So the average population size is (3 + 1 + 2 + 0)/(4-1) = 2.

    :param array birth_times: Birth times of individuals.
    :param array death_times: Death times of individuals. Should be the same length as `birth_times`. Death times are less than or equal to birth times.
    :param float t0: Lower time endpoint.
    :param float t1: Upper time endpoint. t1 is greater than t0. 

    '''
    # # time steps each individual is alive in [t0, t1)
    a = np.maximum(0, np.minimum(t1, 1 + birth_times) -  np.maximum(t0, death_times))
    # Average population size
    return(sum(a)/(t1-t0))

def population_size(ts, x_bins, y_bins, time_bins, stage='late', remembered_stage=None):
    '''
    Calculates the population size in each location bin averaged over the times in each time bin from a tree sequence. The tree sequence must be from a SLiM simulation with all individuals remembered.
    
    Returns a 3-d array with dimensions (`len(x_bins) - 1`,`len(y_bins) - 1`,`len(time_bins)-1`). The i,j,k element of the array is the average population size in [`time_bins[k]`, `time_bins[k + 1]`) for individuals located in [`x_bins[i]`, `x_bins[i + 1]`) and [`y_bins[j]`, `y_bins[j + 1]`). 

    'Average population size' is the sum over all individuals located in the location bin of the length of time each individual was alive in the time bin divided by the length of the time bin. For integer endpoints of the time bins this is equivalent to: for every integer time step in [`time_bins[k]`, `time_bins[k + 1]`), record the number of individuals that are alive at stage `stage` in that time step and are located in [`x_bins[i]`, `x_bins[i + 1]`) and [`y_bins[j]`, `y_bins[j + 1]`), then take the mean of the recorded population sizes.

    :param TreeSequence ts: The tree sequence to calculate population size from.
    :param array x_bins: The x-coordinates of the boundaries of the location bins.
    :param array y_bins: The y-coordinates of the boundaries of the location bins.
    :param array time_bins: The endpoints of the time bins.
    :param str stage: The stage in the SLiM life cycle that the endpoints of the time bins refer to (either "early" or "late"; defaults to "late"). 
    :param str remembered_stage: The stage in the SLiM life cycle during which individuals were Remembered (defaults to the stage the tree sequence was recorded at, stored in metadata).
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
  
    # Initialize array to store popoluationsize
    nxbins = len(x_breaks) - 1
    nybins = len(y_breaks) - 1
    ntbins = len(time_breaks) - 1
    popsize = np.empty((nxbins, nybins, ntbins))
    popsize_validate = np.empty((nxbins, nybins, ntbins))
    #print(np.shape(popsize))

    # Location, times, and ages of individuals
    locations = ts.individual_locations
    times = ts.individual_times

    # Iterate through location bins and time bins
    for i in np.arange(nxbins):
        for j in np.arange(nybins):
            # Endpoints of bins
            x0, x1 = x_breaks[i], x_breaks[i + 1]
            y0, y1 = y_breaks[j], y_breaks[j + 1]
            # Individuals in the location bin
            in_bin = in_location_bin(locations, x0, x1, y0, y1)
            for k in np.arange(ntbins):
                #print(i, j, k)
                # Endpoints of bins
                t0, t1 = time_breaks[k], time_breaks[k + 1]
                #print(x0, x1, y0, y1, t0, t1)
                # Individuals in the location bin
                popsize[i, j, k] = average_time_alive(birth_times[in_bin], death_times[in_bin], t0, t1)
    return(popsize)
