"""
Test cases for tree sequences.
"""
import pickle
import random
import os

import numpy as np
import pytest
import tskit
import msprime
import pyslim

import tests

from .recipe_specs import recipe_eq

class TestPopulationSize(tests.PyslimTestCase):

    def population_size_simple(self, ts, x_bins, y_bins, time_bins, **kwargs):
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
               
        time_breaks = time_bins
        x_breaks = x_bins
        y_breaks = y_bins

        # Initialize array to store popoluationsize
        nxbins = len(x_breaks) - 1
        nybins = len(y_breaks) - 1
        ntbins = len(time_breaks) - 1
        popsize = np.empty((nxbins, nybins, ntbins))

        # Location, times, and ages of individuals
        locations = ts.individual_locations
        times = ts.individual_times

        # Iterate through location bins and time bins
        for i in np.arange(nxbins):
            for j in np.arange(nybins):
                # Endpoints of bins
                x0, x1 = x_breaks[i], x_breaks[i + 1]
                y0, y1 = y_breaks[j], y_breaks[j + 1]
                for k in np.arange(ntbins):
                    # Endpoints of bins
                    t0, t1 = time_breaks[k], time_breaks[k + 1]
                    alive = 0
                    for t in np.arange(np.ceil(t0), t1):
                        for ind_id in ts.individuals_alive_at(t, **kwargs):
                            ind = ts.individual(ind_id)
                            if (ind.location[0] < x1 and ind.location[0] >= x0 and ind.location[1] < y1 and ind.location[1] >= y0):
                                alive += 1/(t1-t0)
                    popsize[i, j, k] = alive     
        return(popsize)

    def make_bins(self, ts):
        for nx in np.arange(1, 20, 10):
            for ny in np.arange(1, 20, 10):
                for nt in np.arange(1, 60, 30):
                    yield [np.linspace(0, round(max(ts.individual_locations[:,0])), nx + 1), 
                           np.linspace(0, round(max(ts.individual_locations[:,1])), ny + 1),
                           np.round(np.linspace(0, max(nt, max(ts.individual_times)), nt + 1))]  

    def verify(self, ts, remembered_stage):
        for bins in self.make_bins(ts):
            for stage in ('early', 'late'):
                x_bins, y_bins, time_bins = bins
                # as computed by pyslim
                popsize0 = pyslim.population_size(ts, x_bins, y_bins, time_bins, stage=stage, remembered_stage=remembered_stage)
                # as computed in a simple way
                popsize1 = self.population_size_simple(ts, x_bins, y_bins, time_bins, stage=stage, remembered_stage=remembered_stage)
                assert(np.allclose(popsize1, popsize0))

    @pytest.mark.parametrize('recipe', [next(recipe_eq("everyone"))], indirect=True)
    def test_errors(self, recipe):
        ts = recipe["ts"]
        x_bins = [0, 1.0]
        y_bins = [0, 1.0]
        time_bins = [0, 10.0]
        for stage in ['abcd', 10, []]:
            with pytest.raises(ValueError):
                pyslim.population_size(ts, x_bins, y_bins, time_bins, stage=stage)
            with pytest.raises(ValueError):
                pyslim.population_size(ts, x_bins, y_bins, time_bins, remembered_stage=stage)

    @pytest.mark.parametrize('recipe', [next(recipe_eq("pedigree", "WF"))], indirect=True)
    def test_mismatched_remembered_stage(self, recipe):
        ts = recipe["ts"]
        x_bins = [0, 1.0]
        y_bins = [0, 1.0]
        time_bins = [0, 10.0]
        info = recipe["info"]
        if "remembered_early" in recipe:
            with pytest.warns(UserWarning):
                pyslim.population_size(ts, x_bins, y_bins, time_bins, remembered_stage="late")
        else:
            with pytest.warns(UserWarning):
                pyslim.population_size(ts, x_bins, y_bins, time_bins, remembered_stage="early")

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_population_size(self, recipe):
        # compare output to the right answer
        ts = recipe["ts"]
        remembered_stage = 'early' if 'remembered_early' in recipe else 'late'
        self.verify(ts, remembered_stage=remembered_stage)

    def test_known_answer(self):
        # a simple example to make sure we've got the edge cases right
        tables = tskit.TableCollection(sequence_length=1)
        pyslim.set_tree_sequence_metadata(tables,  model_type='nonWF', generation=0)
        pyslim.set_metadata_schemas(tables)
        locs = [[0, 0], # alive at 0, 1
                [0, 1], # alive at 0, 1, 2
                [2, 0], # alive at 0, 1, 2
                [1, 1], # alive at 0, 1, 2
                [0, 0], # alive at 1
                [0.5, 1], # alive at 1
                [2, 2], # alive at 1
                [3, 2] # alive at 0, 1, 2, 3
        ]
        births = [1, 2, 2, 2, 1, 1, 1, 3]
        ages = [1, 2, 2, 2, 0, 0, 0, 3]
        x_bins = [0, 1, 3]
        for xy, a in zip(locs, ages):
            md = pyslim.default_slim_metadata('individual')
            md['age'] = a
            tables.individuals.add_row(
                    location=xy + [np.nan],
                    metadata=md,
            )

        for j, b in enumerate(births):
            tables.nodes.add_row(time=b, individual=j)

        ts = pyslim.SlimTreeSequence(
                tables.tree_sequence()
        )

        # check we've got this right
        for k, n in enumerate([
                [0, 1, 2, 3, 7],
                [0, 1, 2, 3, 4, 5, 6, 7],
                [1, 2, 3, 7],
                [7]]):
            np.testing.assert_array_equal(
                    n,
                    ts.individuals_alive_at(k)
            )

        # no-one
        counts = pyslim.population_size(
                ts,
                x_bins=np.arange(10),
                y_bins=np.arange(10),
                time_bins=[100, 200, 300],
        )
        np.testing.assert_array_equal(
                counts,
                np.zeros((9, 9, 2)),
        )

        # everyone at the start
        counts = pyslim.population_size(
                ts,
                x_bins=[0, 10],
                y_bins=[0, 10],
                time_bins=[0, 1],
        )
        np.testing.assert_array_equal(
                counts,
                [[[5]]]
        )

        # should omit the last one
        counts = pyslim.population_size(
                ts,
                x_bins=[0, 3],
                y_bins=[0, 3],
                time_bins=[0, 1]
        )
        np.testing.assert_array_equal(
                counts,
                [[[4]]]
        )

        # now should omit the ones at the boundaries
        counts = pyslim.population_size(
                ts,
                x_bins=[0, 1, 2],
                y_bins=[0, 1, 2],
                time_bins=[0, 1, 2, 5]
        )
        np.testing.assert_array_equal(
                counts,
                [[[1, 2, 0],
                  [1, 2, 1/3]],
                 [[0, 0, 0],
                  [1, 1, 1/3]]]
        )
