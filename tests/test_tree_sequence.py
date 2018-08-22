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
import os


class TestMutate(tests.PyslimTestCase):
    '''
    Tests for pyslim.mutate().
    '''

    def test_mutate(self):
        for ts in self.get_slim_examples():
            mts = pyslim.mutate(ts, rate=1/ts.sequence_length, keep=True)
            self.assertGreaterEqual(mts.num_sites, ts.num_sites)
            self.assertGreaterEqual(mts.num_mutations, ts.num_mutations)


class TestRecapitation(tests.PyslimTestCase):
    '''
    Tests for recapitation.
    '''

    def check_recap_consistency(self, ts, recap,
                                keep_first_generation):
        self.assertEqual(ts.slim_generation, recap.slim_generation)
        sample_map = {}
        k = 0
        for u in ts.samples():
            if (keep_first_generation or
                 ts.node(u).time < ts.slim_generation):
                sample_map[u] = k
                k += 1
        recap_samples = list(recap.samples())
        for u in sample_map:
            n1 = ts.node(u)
            self.assertTrue(sample_map[u] in recap_samples)
            n2 = recap.node(sample_map[u])
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
        for ts in self.get_slim_examples():
            assert ts.num_populations == 2
            # if not we need migration rates
            for keep_first in [True, False]:
                print("keep?", keep_first)
                recap = ts.recapitate(recombination_rate = 1.0,
                                      keep_first_generation=keep_first)
                # there should be no new mutations
                self.assertEqual(ts.num_mutations, recap.num_mutations)
                self.assertEqual(ts.num_sites, recap.num_sites)
                self.assertListEqual(list(ts.tables.sites.position),
                                     list(recap.tables.sites.position))
                self.check_recap_consistency(ts, recap, keep_first)
                for t in recap.trees():
                    self.assertEqual(t.num_roots, 1)

                recap = ts.recapitate(recombination_rate = 1.0,
                                      Ne = 1e-6, keep_first_generation=keep_first)
                self.check_recap_consistency(ts, recap, keep_first)
                for t in recap.trees():
                    self.assertEqual(t.num_roots, 1)
                    self.assertAlmostEqual(recap.node(t.root).time, 
                                           recap.slim_generation, 
                                           delta = 1e-4)

    def test_simplify(self):
        for ts in self.get_slim_examples():
            sts = ts.simplify([0, 1])
            self.assertEqual(ts.slim_generation, sts.slim_generation)
            mtabs = sts.dump_tables()
            mtabs.simplify([0, 1])
            mts = mtabs.tree_sequence()
            self.assertEqual(mts.tables, sts.tables)

