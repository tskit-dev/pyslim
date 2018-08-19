"""
Test cases for recapitation.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import tests
import unittest
import random
import _msprime

def get_slim_example_files():
    for filename in ["tests/examples/recipe_WF.trees",
                     "tests/examples/recipe_nonWF.trees"]:
        yield filename

def get_slim_examples():
    for filename in get_slim_example_files():
        yield pyslim.load(filename, slim_format=True)


class TestRecapitation(tests.PyslimTestCase):
    '''
    Tests for recapitation.
    '''

    def check_recap_consistency(self, ts, recap):
        self.assertEqual(ts.num_samples, recap.num_samples)
        self.assertLessEqual(ts.num_nodes, recap.num_nodes)
        for k in range(ts.num_nodes):
            n1 = ts.node(k)
            n2 = recap.node(k)
            self.assertEqual(n1.time, n2.time)
            self.assertEqual(n1.individual, n2.individual)
            self.assertEqual(n1.flags, n2.flags)
            self.assertEqual(n1.metadata, n2.metadata)
            self.assertEqual(n1.population, n2.population)
        self.assertLessEqual(ts.num_populations, recap.num_populations)
        for k in range(ts.num_populations):
            p1 = ts.population(k)
            p2 = recap.population(k)
            self.assertEqual(p1.metadata, p2.metadata)

    def test_recapitation(self):
        for ts in get_slim_examples():
            assert ts.num_populations == 2
            # if not we need migration rates
            recap = ts.recapitate(recombination_rate = 1.0)
            # there should be no new mutations
            self.assertEqual(ts.tables.mutations, recap.tables.mutations)
            self.assertEqual(ts.tables.sites, recap.tables.sites)
            self.check_recap_consistency(ts, recap)
            for t in recap.trees():
                self.assertEqual(t.num_roots, 1)

            recap = ts.recapitate(recombination_rate = 1.0,
                                  Ne = 1e-6)
            self.check_recap_consistency(ts, recap)
            for t in recap.trees():
                self.assertEqual(t.num_roots, 1)
                self.assertAlmostEqual(recap.node(t.root).time, 
                                       recap.slim_generation, 
                                       delta = 1e-4)

    def test_simplify(self):
        for ts in get_slim_examples():
            sts = ts.simplify([0, 1])
            self.assertEqual(ts.slim_generation, sts.slim_generation)
            mtabs = sts.dump_tables()
            mtabs.simplify([0, 1])
            mts = mtabs.tree_sequence()
            self.assertEqual(mts.tables, sts.tables)

