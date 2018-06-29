"""
Test cases for the metadata reading/writing of pyslim.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import tests
import unittest
import random


def get_msprime_examples():
    for n in [2, 10, 100]:
        for mutrate in [0.0]:
            for recrate in [0.0, 1.0]:
                yield msprime.simulate(n, mutation_rate=mutrate,
                                       recombination_rate=recrate)

class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests for tools to annotate existing msprime-derived tree sequences.
    '''

    def verify_annotated_tables(self, ts1, ts2, time_offset):
        '''
        Verify that the tables returned after annotation are equal, up to the
        expected forgetting of metadata.
        '''
        tables1 = ts1.tables
        tables2 = ts2.tables
        # compare nodes
        self.assertArrayEqual(tables1.nodes.flags, tables2.nodes.flags)
        self.assertArrayAlmostEqual(tables1.nodes.time, tables2.nodes.time + time_offset)
        self.assertArrayEqual(tables1.nodes.population, tables2.nodes.population)
        # compare edges
        self.assertEqual(tables1.edges, tables2.edges)
        # compare sites
        self.assertArrayEqual(tables1.sites.position, tables2.sites.position)
        self.assertArrayEqual(tables1.sites.ancestral_state, tables2.sites.ancestral_state)
        self.assertArrayEqual(tables1.sites.ancestral_state_offset,
                              tables2.sites.ancestral_state_offset)
        # compare mutations
        self.assertArrayEqual(tables1.mutations.site, tables2.mutations.site)
        self.assertArrayEqual(tables1.mutations.node, tables2.mutations.node)
        self.assertArrayEqual(tables1.mutations.derived_state, tables2.mutations.derived_state)
        self.assertArrayEqual(tables1.mutations.derived_state_offset,
                              tables2.mutations.derived_state_offset)

    def verify_annotated_trees(self, ts1, ts2):
        '''
        Verify the *trees* returned before and after annotation are equal.
        '''
        self.assertEqual(ts1.num_trees, ts2.num_trees)
        for t1, t2 in zip(ts1.trees(), ts2.trees()):
            self.assertEqual(t1.length, t2.length)
            self.assertEqual(t1.get_parent_dict(), t2.get_parent_dict())
            self.assertAlmostEqual(t1.total_branch_length, t2.total_branch_length)

    def verify_alleles(self, ts, slim_ts):
        '''
        Verify that haplotypes agree between tree sequences, after translation
        through slim_ts.alleles.
        '''

    def test_basic_annotation(self):
        for ts in get_msprime_examples():
            slim_gen = 1
            slim_ts = pyslim.annotate_defaults(ts, model_type="WF",
                                               slim_generation=slim_gen)
            self.verify_annotated_tables(ts, slim_ts, time_offset=slim_gen)
            self.verify_annotated_trees(ts, slim_ts)
            self.verify_alleles(ts, slim_ts)

    def test_annotate_individuals(self):
        for ts in get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
            metadata = list(pyslim.extract_individual_metadata(tables))
            self.assertEqual(len(metadata), slim_ts.num_individuals)
            sexes = [random.choice([pyslim.INDIVIDUAL_TYPE_FEMALE, pyslim.INDIVIDUAL_TYPE_MALE])
                     for _ in metadata]
            for j in range(len(metadata)):
                metadata[j].sex = sexes[j]
            pyslim.annotate_individual_metadata(tables, metadata)
            new_ts = pyslim.load_tables(tables, slim_format=True)
            for j, ind in enumerate(new_ts.individuals()):
                md = pyslim.decode_individual(ind.metadata)
                self.assertEqual(md.sex, sexes[j])

    def test_annotate_mutations(self):
        for ts in get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
            metadata = list(pyslim.extract_mutation_metadata(tables))
            self.assertEqual(len(metadata), slim_ts.num_mutations)
            selcoefs = [random.uniform(0, 1) for _ in metadata]
            for j in range(len(metadata)):
                metadata[j].selection_coeff = selcoefs[j]
            pyslim.annotate_mutation_metadata(tables, metadata)
            new_ts = pyslim.load_tables(tables, slim_format=True)
            for j, x in enumerate(new_ts.mutations()):
                md = pyslim.decode_mutation(x.metadata)
                self.assertEqual(md.selection_coeff, selcoefs[j])

    def test_annotate_nodes(self):
        for ts in get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
            metadata = list(pyslim.extract_node_metadata(tables))
            self.assertEqual(len(metadata), slim_ts.num_nodes)
            gtypes = [random.choice([pyslim.GENOME_TYPE_X, pyslim.GENOME_TYPE_Y])
                      for _ in metadata]
            for j in range(len(metadata)):
                if metadata[j] is not None:
                    metadata[j].genome_type = gtypes[j]
            pyslim.annotate_node_metadata(tables, metadata)
            new_ts = pyslim.load_tables(tables, slim_format=True)
            for j, x in enumerate(new_ts.nodes()):
                md = pyslim.decode_node(x.metadata)
                if md is not None:
                    self.assertEqual(md.genome_type, gtypes[j])

