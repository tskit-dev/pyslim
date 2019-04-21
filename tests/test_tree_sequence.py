"""
Test cases for recapitation.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import _msprime
import tests
import unittest
import random
import os


class TestRecapitation(tests.PyslimTestCase):
    '''
    Tests for recapitation.
    '''

    def check_recap_consistency(self, ts, recap,
                                keep_first_generation):
        self.assertEqual(ts.slim_generation, recap.slim_generation)
        ts_samples = list(ts.samples())
        for u in recap.samples():
            n1 = recap.node(u)
            self.assertGreaterEqual(n1.individual, 0)
            i1 = recap.individual(n1.individual)
            firstgen = ((pyslim.INDIVIDUAL_FIRST_GEN & i1.flags) > 0)
            remembered = ((pyslim.INDIVIDUAL_REMEMBERED & i1.flags) > 0)
            alive = ((pyslim.INDIVIDUAL_ALIVE & i1.flags) > 0)
            self.assertTrue(alive or remembered or (firstgen and keep_first_generation))
            self.assertTrue((keep_first_generation and firstgen) or (u in ts_samples))
            n2 = ts.node(u)
            self.assertEqual(n1.time, n2.time)
            self.assertEqual(n1.individual, n2.individual)
            recap_flags = n1.flags
            if keep_first_generation and firstgen and not remembered and not alive:
                recap_flags = (recap_flags & ~msprime.NODE_IS_SAMPLE)
            self.assertEqual(recap_flags, n2.flags)
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


class TestSimplify(tests.PyslimTestCase):
    '''
    Our simplify() is just a wrapper around the tskit simplify.
    '''

    def test_simplify(self):
        for ts in self.get_slim_examples():
            sts = ts.simplify()
            self.assertEqual(ts.sequence_length, sts.sequence_length)
            self.assertEqual(type(ts), type(sts))
            self.assertEqual(sts.samples()[0], 0)    

class TestReferenceSequence(tests.PyslimTestCase):
    '''
    Test for operations involving the reference sequence
    '''
    def test_reference_sequence(self):
        for ts in self.get_slim_examples():
            mut_md = pyslim.decode_mutation(ts.mutation(0).metadata)
            has_nucleotides = (mut_md[0].nucleotide >= 0)
            if not has_nucleotides:
                self.assertEqual(ts.reference_sequence, None)
            else:
                self.assertEqual(type(ts.reference_sequence), type(''))
                self.assertEqual(len(ts.reference_sequence), ts.sequence_length)
            sts = ts.simplify(ts.samples()[:2])
            self.assertEqual(sts.reference_sequence, ts.reference_sequence)


