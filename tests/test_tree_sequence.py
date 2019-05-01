"""
Test cases for tree sequences.
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
import numpy as np


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


class TestIndividualMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Individuals.
    '''

    def test_individual_derived_info(self):
        for ts in self.get_slim_examples():
            for j, ind in enumerate(ts.individuals()):
                a = ts.tables.individuals.metadata_offset[j]
                b = ts.tables.individuals.metadata_offset[j+1]
                raw_md = ts.tables.individuals.metadata[a:b]
                md = pyslim.decode_individual(raw_md)
                self.assertEqual(ind.metadata, md)
                for n in ind.nodes:
                    self.assertEqual(ts.node(n).population, ind.population)
                    self.assertEqual(ts.node(n).time, ind.time)

    def test_individual_embellishments(self):
        '''
        Test the individual additional information.
        '''
        for ts in self.get_slim_examples():
            for j, ind in enumerate(ts.individuals()):
                self.assertEqual(ts.individual_times[j], ind.time)
                self.assertEqual(ts.individual_ages[j], ind.metadata.age)
                self.assertEqual(ts.individual_populations[j], ind.population)
                self.assertArrayEqual(ts.individual_locations[j], ind.location)


class TestMutationMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Mutations.
    '''

    def test_mutation_derived_info(self):
        for ts in self.get_slim_examples():
            for j, mut in enumerate(ts.mutations()):
                a = ts.tables.mutations.metadata_offset[j]
                b = ts.tables.mutations.metadata_offset[j+1]
                raw_md = ts.tables.mutations.metadata[a:b]
                md = pyslim.decode_mutation(raw_md)
                self.assertEqual(mut.metadata, md)


class TestEveryone(tests.PyslimTestCase):
    '''
    Test things we calculate from tree sequences where 
    everyone is remembered at all times.
    '''

    def test_alive_ages(self):
        for ts in self.get_slim_everyone_examples():
            ages_mat = np.zeros((ts.num_individuals, ts.slim_generation))
            alive_mat = np.zeros((ts.num_individuals, ts.slim_generation))
            alive_now = ts.individuals_alive_at(0)
            for time in range(ts.slim_generation):
                ages_mat[:, time] = ts.individual_ages_at(time)
                for j in ts.individuals_alive_at(time):
                    alive_mat[j, time] = 1

            for j, ind in enumerate(ts.individuals()):
                self.assertEqual(j in alive_now,
                                 ind.flags & pyslim.INDIVIDUAL_ALIVE > 0)
                self.assertEqual(sum(alive_mat[j, :]), 1 + ind.metadata.age)
                self.assertEqual(alive_mat[j, int(ind.time)], 1)
                self.assertEqual(ages_mat[j, int(ind.time)], 0)
                for t in range(ts.slim_generation):
                    if t > ind.time:
                        self.assertTrue(np.isnan(ages_mat[j, t]))
                    if t == ind.time:
                        self.assertEqual(ages_mat[j, t], 0)
                if ind.time + 1 < ts.slim_generation:
                    self.assertEqual(alive_mat[j, int(ind.time + 1)], 0)


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
            mut_md = ts.mutation(0).metadata
            has_nucleotides = (mut_md[0].nucleotide >= 0)
            if not has_nucleotides:
                self.assertEqual(ts.reference_sequence, None)
            else:
                self.assertEqual(type(ts.reference_sequence), type(''))
                self.assertEqual(len(ts.reference_sequence), ts.sequence_length)
                for u in ts.reference_sequence:
                    self.assertTrue(u in pyslim.NUCLEOTIDES)
            sts = ts.simplify(ts.samples()[:2])
            self.assertEqual(sts.reference_sequence, ts.reference_sequence)
