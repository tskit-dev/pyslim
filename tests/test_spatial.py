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

class TestPopulationSize(tests.PyslimTestCase):

    def verify(self, ts):
        # as computed by pyslim
        popsize0 = pyslim.population_size(ts)
        # as computed by testing code
        popsize1 = # compute in a simple way here
        assert np.array_equal(popsize0, popsize1)

    def test_population_size_errors(self):
        # test that it fails appropriately on bad input

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_population_size(self):
        # compare output to the right answer
        ts = recipe["ts"]
        self.verify(ts)
