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
import json


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

    def verify_annotated_tables(self, ts1, ts2):
        '''
        Verify that the tables returned after annotation are equal, up to the
        expected forgetting of metadata.
        '''
        tables1 = ts1.tables
        tables2 = ts2.tables
        # compare nodes
        self.assertArrayEqual(tables1.nodes.flags, tables2.nodes.flags)
        self.assertArrayAlmostEqual(tables1.nodes.time, tables2.nodes.time)
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

    def verify_consistency(self, ts):
        '''
        Check that individuals exist, and populations agree between nodes and individuals.
        '''

    def verify_defaults(self, ts):
        '''
        Verify the default values have been entered into metadata.
        '''
        mut_md = pyslim.extract_mutation_metadata(ts.tables)
        for md in mut_md:
            self.assertEqual(md.mutation_type, 1)
            self.assertEqual(md.selection_coeff, 0.0)
            self.assertEqual(md.population, msprime.NULL_POPULATION)
            self.assertEqual(md.slim_time, 0)
        node_md = pyslim.extract_node_metadata(ts.tables)
        for md, node in zip(node_md, ts.nodes()):
            if not node.is_sample():
                self.assertEqual(md, None)
            else:
                self.assertEqual(md.is_null, False)
                self.assertEqual(md.genome_type, pyslim.GENOME_TYPE_AUTOSOME)
        for ind in ts.individuals(): 
            self.assertArrayEqual(ind.location, [0, 0, 0])
            self.assertEqual(ind.flags, pyslim.INDIVIDUAL_ALIVE)
        ind_md = pyslim.extract_individual_metadata(ts.tables)
        for md in ind_md:
            self.assertEqual(md.sex, pyslim.INDIVIDUAL_TYPE_HERMAPHRODITE)
            self.assertEqual(md.flags, 0)
        pop_md = pyslim.extract_population_metadata(ts.tables)
        for md in pop_md:
            self.assertEqual(md.selfing_fraction, 0.0)
            self.assertEqual(md.female_cloning_fraction, 0.0)
            self.assertEqual(md.male_cloning_fraction, 0.0)
            self.assertEqual(md.sex_ratio, 0.5)
            self.assertEqual(md.bounds_x0, 0.0)
            self.assertEqual(md.bounds_x1, 0.0)
            self.assertEqual(md.bounds_y0, 0.0)
            self.assertEqual(md.bounds_y1, 0.0)
            self.assertEqual(md.bounds_z0, 0.0)
            self.assertEqual(md.bounds_z1, 0.0)
            self.assertEqual(len(md.migration_records), 0)

    def verify_provenance(self, ts):
        for u in ts.provenances():
            msprime.validate_provenance(json.loads(u.record))

    def test_basic_annotation(self):
        for ts in get_msprime_examples():
            slim_gen = 4
            slim_ts = pyslim.annotate_defaults(ts, model_type="WF",
                                               slim_generation=slim_gen)
            self.verify_annotated_tables(ts, slim_ts)
            self.verify_annotated_trees(ts, slim_ts)
            self.verify_haplotype_equality(ts, slim_ts)
            self.verify_defaults(slim_ts)
            self.verify_provenance(slim_ts)

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
            new_ts = pyslim.load_tables(tables)
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
            new_ts = pyslim.load_tables(tables)
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
            new_ts = pyslim.load_tables(tables)
            for j, x in enumerate(new_ts.nodes()):
                md = pyslim.decode_node(x.metadata)
                if md is not None:
                    self.assertEqual(md.genome_type, gtypes[j])

