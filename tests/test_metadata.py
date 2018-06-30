"""
Test cases for the metadata reading/writing of pyslim.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import tests
import unittest

def get_slim_example_files():
    for filename in ["tests/examples/recipe_WF.trees",
                     "tests/examples/recipe_nonWF.trees"]:
        yield filename

def get_slim_examples():
    for filename in get_slim_example_files():
        yield pyslim.load(filename, slim_format=True)


class TestEncodeDecode(tests.PyslimTestCase):
    '''
    Tests for conversion to/from binary representations of metadata.
    '''

    def test_mutation_metadata(self):
        for md_length in [0, 1, 5]:
            md = [pyslim.MutationMetadata(
                     mutation_type=j, selection_coeff=0.5, population=j,
                     slim_time=10 + j) for j in range(md_length)]
            md_bytes = pyslim.encode_mutation(md)
            new_md = pyslim.decode_mutation(md_bytes)
            self.assertEqual(len(md), len(new_md))
            for x, y in zip(md, new_md):
                self.assertEqual(x, y)

    def test_node_metadata(self):
        md = pyslim.NodeMetadata(slim_id=2, is_null=False,
                                 genome_type=pyslim.GENOME_TYPE_X)
        md_bytes = pyslim.encode_node(md)
        new_md = pyslim.decode_node(md_bytes)
        self.assertEqual(md, new_md)

    def test_individual_metadata(self):
        md = pyslim.IndividualMetadata(
                age=2, pedigree_id=23, population=0,
                sex=pyslim.INDIVIDUAL_TYPE_MALE,
                flags=pyslim.INDIVIDUAL_FLAG_MIGRATED)
        md_bytes = pyslim.encode_individual(md)
        new_md = pyslim.decode_individual(md_bytes)
        self.assertEqual(md, new_md)

    def test_population_metadata(self):
        mrs = [pyslim.PopulationMigrationMetadata(source_subpop=j, migration_rate=0.2)
               for j in range(3)]
        for mr_list in [[], mrs]:
            md = pyslim.PopulationMetadata(
                    slim_id=1, selfing_fraction=0.75, female_cloning_fraction=0.2,
                    male_cloning_fraction=0.8, sex_ratio=0.6, bounds_x0=-1.0,
                    bounds_x1=2.0, bounds_y0=0.0, bounds_y1=0.0, bounds_z0=0.0,
                    bounds_z1=1.0, migration_records=mr_list)
            md_bytes = pyslim.encode_population(md)
            new_md = pyslim.decode_population(md_bytes)
            self.assertEqual(md, new_md)


class TestAnnotate(unittest.TestCase):
    '''
    Tests the table annotation methods.
    '''

    def test_annotate_mutations(self):
        for ts in get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in msprime.unpack_bytes(tables.mutations.metadata, 
                                           tables.mutations.metadata_offset):
                dm = pyslim.decode_mutation(md)
                edm = pyslim.encode_mutation(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_mutation_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_nodes(self):
        for ts in get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in msprime.unpack_bytes(tables.nodes.metadata, 
                                           tables.nodes.metadata_offset):
                dm = pyslim.decode_node(md)
                edm = pyslim.encode_node(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_node_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_individuals(self):
        for ts in get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in msprime.unpack_bytes(tables.individuals.metadata, 
                                           tables.individuals.metadata_offset):
                dm = pyslim.decode_individual(md)
                edm = pyslim.encode_individual(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_individual_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_populations(self):
        for ts in get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in msprime.unpack_bytes(tables.populations.metadata, 
                                           tables.populations.metadata_offset):
                dm = pyslim.decode_population(md)
                edm = pyslim.encode_population(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_population_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)


class TestDumpLoad(tests.PyslimTestCase):
    '''
    Test reading and writing.
    '''

    def verify_time_offset(self, ts, slim_ts):
        gen = slim_ts.slim_generation
        self.assertEqual(ts.num_nodes, slim_ts.num_nodes)
        for n1, n2 in zip(ts.nodes(), slim_ts.nodes()):
            self.assertEqual(n1.time, n2.time - gen)

    def test_load_tables(self):
        for ts in get_slim_examples():
            self.assertTrue(type(ts) is pyslim.SlimTreeSequence)
            tables = ts.tables
            new_ts = pyslim.load_tables(tables, slim_format=True)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            new_tables = new_ts.tables
            self.assertTablesAlmostEqual(tables, new_tables)

    def test_load(self):
        for fn in get_slim_example_files():
            msp_ts = pyslim.load(fn, slim_format=False)
            self.assertTrue(type(msp_ts) is msprime.TreeSequence)
            tables = msp_ts.tables
            new_ts = pyslim.load_tables(tables, slim_format=True)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            self.verify_time_offset(msp_ts, new_ts)
            new_tables = new_ts.tables
            self.assertTablesAlmostEqual(tables, new_tables)
            slim_ts = pyslim.load(fn, slim_format=True)
            self.assertTrue(type(slim_ts) is pyslim.SlimTreeSequence)
            slim_tables = slim_ts.tables
            self.assertTablesAlmostEqual(slim_tables, new_tables)
            new_msp_ts = pyslim.load_tables(tables, slim_format=False)
            self.assertTrue(type(new_msp_ts) is msprime.TreeSequence)
            self.assertTablesAlmostEqual(tables, new_msp_ts.tables)


class TestAlleles(tests.PyslimTestCase):
    '''
    Test allele translation.
    '''

    def test_haplotypes(self):
        for slim_ts in get_slim_examples():
            tables = slim_ts.tables
            ts = tables.tree_sequence()
            self.verify_haplotype_equality(ts, slim_ts)

