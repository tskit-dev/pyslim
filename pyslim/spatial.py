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
    Takes the locations of individuals and  tests if each is within [x0, x1), [y0, y1). 
    Returns a vector the same length as locations with True if they are within the bounds and False otherwise
    locations is a n x 3 array with rows individuals and [:, 0] the x coordinate and [:, 1] the y coordinate
    Returns an n x 1 array, TRUE if individual is in the bin
    Each row corresponds to an individual
    '''
    return(np.logical_and(np.logical_and(locations[:,0] < x1, locations[:,0] >= x0),
               np.logical_and(locations[:,1] < y1, locations[:,1] >= y0)))

def average_time_alive(birth_times, death_times, t0, t1):
    '''
    Finds average population size in [t0, t1) using time alive
    '''
    # # time steps each individual is alive in [t0, t1)
    a = np.maximum(0, np.minimum(t1, 1 + birth_times) -  np.maximum(t0, death_times))
    # Average population size
    return(sum(a)/(t1-t0))

def population_size(ts, x_bins, y_bins, time_bins, stage='late', remembered_stage=None):
    '''
    Calculates population size in each location bin averaged over each time_bin.   
    '''
    
    # Want to return
    # [[[px0y0t0, px0y0t1, ...], [px0y1t0, px0y1t1, ...], ...],
    #  [[px1y0t0, px0y0t1, ...], [px1y1t0, px1y1t1, ...], ...],
    #  [[px2y0t0, px0y0t1, ...], [px2y1t0, px2y1t1, ...], ...]]
    # where pxiyjtk is the population size in [x_breaks[i], x_breaks[i + 1]) and [y_breaks[j], x_breaks[j + 1]),
    # averaged over each time point in [t[k], t[k + 1])
    # different offsets for different stages, from pyslim.individuals_alive_at code

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
