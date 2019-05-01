"""
Test cases for the metadata reading/writing of pyslim.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import tskit
import tests
import unittest
import os
import tempfile


class TestEncodeDecode(tests.PyslimTestCase):
    '''
    Tests for conversion to/from binary representations of metadata.
    '''

    def test_decode_errors(self):
        with self.assertRaises(ValueError):
            pyslim.decode_mutation(2.0)
        with self.assertRaises(ValueError):
            pyslim.decode_mutation([2.0, 3.0])
        with self.assertRaises(ValueError):
            pyslim.decode_node(2.0)
        with self.assertRaises(ValueError):
            pyslim.decode_node([1, 2])
        with self.assertRaises(ValueError):
            pyslim.decode_individual(3.0)
        with self.assertRaises(ValueError):
            pyslim.decode_individual([1, 2])
        with self.assertRaises(ValueError):
            pyslim.decode_population(1.0)
        with self.assertRaises(ValueError):
            pyslim.decode_population([2, 3])

    def test_decode_already_mutation(self):
        m = [pyslim.MutationMetadata(mutation_type = 0,
                                     selection_coeff = 0.2,
                                     population = k,
                                     slim_time = 130,
                                     nucleotide = 2) for k in range(4)]
        dm = pyslim.decode_mutation(m)
        self.assertEqual(type(dm), type([]))
        for a, b in zip(m, dm):
            self.assertEqual(a, b)

    def test_decode_already_node(self):
        m = pyslim.NodeMetadata(slim_id=123, is_null=True, genome_type=0)
        dm = pyslim.decode_node(m)
        self.assertEqual(m, dm)

    def test_decode_already_population(self):
        m = pyslim.PopulationMetadata(slim_id=1, selfing_fraction=0.2,
                                      female_cloning_fraction=0.3,
                                      male_cloning_fraction=0.4,
                                      sex_ratio=0.5, bounds_x0=0, bounds_x1=10,
                                      bounds_y0=2, bounds_y1=20, bounds_z0=3,
                                      bounds_z1=30, migration_records=[])
        dm = pyslim.decode_population(m)
        self.assertEqual(m, dm)

    def test_decode_already_individual(self):
        m = pyslim.IndividualMetadata(pedigree_id=24, age=8, population=1,
                                      sex=1, flags=0)
        dm = pyslim.decode_individual(m)
        self.assertEqual(m, dm)

    def test_mutation_metadata(self):
        for md_length in [0, 1, 5]:
            md = [pyslim.MutationMetadata(
                     mutation_type=j, selection_coeff=0.5, population=j,
                     slim_time=10 + j, nucleotide=(j % 5) - 1) 
                     for j in range(md_length)]
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


class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests the table annotation methods.
    '''

    def test_annotate_mutations(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.mutations.metadata,
                                           tables.mutations.metadata_offset):
                dm = pyslim.decode_mutation(md)
                edm = pyslim.encode_mutation(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_mutation_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_nodes(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.nodes.metadata,
                                           tables.nodes.metadata_offset):
                dm = pyslim.decode_node(md)
                edm = pyslim.encode_node(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_node_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_individuals(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.individuals.metadata,
                                           tables.individuals.metadata_offset):
                dm = pyslim.decode_individual(md)
                edm = pyslim.encode_individual(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            pyslim.annotate_individual_metadata(new_tables, metadata)
            self.assertEqual(tables, new_tables)

    def test_annotate_populations(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.populations.metadata,
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

    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_ll_ts_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)

    def verify_times(self, ts, slim_ts):
        gen = slim_ts.slim_generation
        self.assertEqual(ts.num_nodes, slim_ts.num_nodes)
        # verify internal consistency
        for j in range(slim_ts.num_nodes):
            self.assertEqual(slim_ts.node(j).time,
                             slim_ts.tables.nodes.time[j])
        # verify consistency between tree sequences
        for n1, n2 in zip(ts.nodes(), slim_ts.nodes()):
            self.assertEqual(n1.time, n2.time)

    def verify_dump_equality(self, ts):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(self.temp_file)
        ts2 = pyslim.load(self.temp_file)
        self.assertEqual(ts.num_samples, ts2.num_samples)
        self.assertEqual(ts.sequence_length, ts2.sequence_length)
        self.assertEqual(ts.tables, ts2.tables)

    def test_load_tables(self):
        for ts in self.get_slim_examples():
            self.assertTrue(type(ts) is pyslim.SlimTreeSequence)
            tables = ts.tables
            new_ts = pyslim.load_tables(tables)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            new_tables = new_ts.tables
            self.assertEqual(tables, new_tables)

    def test_load(self):
        for fn in self.get_slim_example_files():
            # load in msprime then switch
            msp_ts = tskit.load(fn)
            self.assertTrue(type(msp_ts) is msprime.TreeSequence)
            # transfer tables
            msp_tables = msp_ts.tables
            new_ts = pyslim.load_tables(msp_tables)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            self.verify_times(msp_ts, new_ts)
            self.assertEqual(msp_tables, new_ts.tables)
            # convert directly
            new_ts = pyslim.SlimTreeSequence(msp_ts)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            self.verify_times(msp_ts, new_ts)
            self.assertEqual(msp_tables, new_ts.tables)
            # load to pyslim from file
            slim_ts = pyslim.load(fn)
            self.assertTrue(type(slim_ts) is pyslim.SlimTreeSequence)
            self.assertEqual(msp_tables, slim_ts.tables)
            self.assertEqual(slim_ts.slim_generation, new_ts.slim_generation)

    def test_dump_equality(self):
        for ts in self.get_slim_examples():
            self.verify_dump_equality(ts)


class TestAlleles(tests.PyslimTestCase):
    '''
    Test nothing got messed up with haplotypes.
    '''

    def test_haplotypes(self):
        for slim_ts in self.get_slim_examples():
            tables = slim_ts.tables
            ts = tables.tree_sequence()
            self.verify_haplotype_equality(ts, slim_ts)


class TestNucleotides(tests.PyslimTestCase):
    '''
    Test nucleotide support
    '''

    def check_nucleotides(self, ts):
        '''
        Check that nucleotides are all valid, i.e.,
        -1, 0, 1, 2, or 3.
        '''
        for mut in ts.mutations():
            for u in mut.metadata:
                self.assertGreaterEqual(u.nucleotide, -1)
                self.assertLessEqual(u.nucleotide, 3)

    def test_nucleotides(self):
        for ts in self.get_slim_examples():
            self.check_nucleotides(ts)

