"""
Test cases for the deprecated, legacy metadata representation of pyslim.
"""
from __future__ import print_function
from __future__ import division

import unittest
import os
import tempfile
import numpy as np
import warnings
import msprime
import tskit

import pyslim

import tests


class LegacyPyslimTestCase(tests.PyslimTestCase):

    # copied from __init__ but with legacy_metadata=True
    def get_slim_examples(self, return_info=False, **kwargs):
        for ex in tests.example_files.values():
            basename = ex['basename']
            use = True
            for a in kwargs:
                if a not in ex or ex[a] != kwargs[a]:
                    use = False
            if use:
                treefile = basename + ".trees"
                print("---->", treefile)
                self.assertTrue(os.path.isfile(treefile))
                ts = pyslim.load(treefile, legacy_metadata=True)
                if return_info:
                    infofile = treefile + ".pedigree"
                    if os.path.isfile(infofile):
                        ex['info'] = self.get_slim_info(infofile)
                    else:
                        ex['info'] = None
                    yield (ts, ex)
                else:
                    yield ts


class TestLegacyTypes(LegacyPyslimTestCase):

    def test_test_warnings(self):
        # check that our checking for warnings works
        # (it didn't with DeprecationWarnings)
        with self.assertWarns(FutureWarning):
            warnings.warn('hi', FutureWarning)

    def test_legacy_types(self):
        for ts in self.get_slim_examples():
            self.assertTrue(ts.legacy_metadata)
            for pop in ts.populations():
                if pop.metadata is not None:
                    self.assertTrue(isinstance(pop.metadata, pyslim.PopulationMetadata))
                    break
            for ind in ts.individuals():
                self.assertTrue(isinstance(ind.metadata, pyslim.IndividualMetadata))
                break
            for n in ts.nodes():
                if n.metadata is not None:
                    self.assertTrue(isinstance(n.metadata, pyslim.NodeMetadata))
                    break
            for mut in ts.mutations():
                self.assertTrue(isinstance(mut.metadata, list))
                if len(mut.metadata) > 0:
                    for u in mut.metadata:
                        self.assertTrue(isinstance(u, pyslim.MutationMetadata))
                    break


class TestEncodeDecode(LegacyPyslimTestCase):
    '''
    Tests for conversion to/from binary representations of metadata.
    '''

    def test_legacy_errors(self):
        defaults = pyslim.default_slim_metadata
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.decode_mutation(defaults['mutation'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.decode_population(defaults['population'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.decode_individual(defaults['individual'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.decode_node(defaults['node'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.encode_mutation(defaults['mutation'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.encode_population(defaults['population'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.encode_individual(defaults['individual'])
        with self.assertRaisesRegex(ValueError, "legacy"):
            pyslim.encode_node(defaults['node'])

    def test_decode_errors(self):
        with self.assertWarns(FutureWarning):
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
        with self.assertWarns(FutureWarning):
            dm = pyslim.decode_mutation(m)
        self.assertEqual(type(dm), type([]))
        for a, b in zip(m, dm):
            self.assertEqual(a, b)

    def test_decode_already_node(self):
        m = pyslim.NodeMetadata(slim_id=123, is_null=True, genome_type=0)
        with self.assertWarns(FutureWarning):
            dm = pyslim.decode_node(m)
        self.assertEqual(m, dm)

    def test_decode_already_population(self):
        m = pyslim.PopulationMetadata(slim_id=1, selfing_fraction=0.2,
                                      female_cloning_fraction=0.3,
                                      male_cloning_fraction=0.4,
                                      sex_ratio=0.5, bounds_x0=0, bounds_x1=10,
                                      bounds_y0=2, bounds_y1=20, bounds_z0=3,
                                      bounds_z1=30, migration_records=[])
        with self.assertWarns(FutureWarning):
            dm = pyslim.decode_population(m)
        self.assertEqual(m, dm)

    def test_decode_already_individual(self):
        m = pyslim.IndividualMetadata(pedigree_id=24, age=8, population=1,
                                      sex=1, flags=0)
        with self.assertWarns(FutureWarning):
            dm = pyslim.decode_individual(m)
        self.assertEqual(m, dm)

    def test_mutation_metadata(self):
        for md_length in [0, 1, 5]:
            md = [pyslim.MutationMetadata(
                     mutation_type=j, selection_coeff=0.5, population=j,
                     slim_time=10 + j, nucleotide=(j % 5) - 1) 
                     for j in range(md_length)]
            with self.assertWarns(FutureWarning):
                md_bytes = pyslim.encode_mutation(md)
            with self.assertWarns(FutureWarning):
                new_md = pyslim.decode_mutation(md_bytes)
            self.assertEqual(len(md), len(new_md))
            for x, y in zip(md, new_md):
                self.assertEqual(x, y)

    def test_node_metadata(self):
        md = pyslim.NodeMetadata(slim_id=2, is_null=False,
                                 genome_type=pyslim.GENOME_TYPE_X)
        with self.assertWarns(FutureWarning):
            md_bytes = pyslim.encode_node(md)
        with self.assertWarns(FutureWarning):
            new_md = pyslim.decode_node(md_bytes)
        self.assertEqual(md, new_md)

    def test_individual_metadata(self):
        md = pyslim.IndividualMetadata(
                age=2, pedigree_id=23, population=0,
                sex=pyslim.INDIVIDUAL_TYPE_MALE,
                flags=pyslim.INDIVIDUAL_FLAG_MIGRATED)
        with self.assertWarns(FutureWarning):
            md_bytes = pyslim.encode_individual(md)
        with self.assertWarns(FutureWarning):
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
            with self.assertWarns(FutureWarning):
                md_bytes = pyslim.encode_population(md)
            with self.assertWarns(FutureWarning):
                new_md = pyslim.decode_population(md_bytes)
            self.assertEqual(md, new_md)


class TestAnnotate(LegacyPyslimTestCase):
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
                with self.assertWarns(FutureWarning):
                    dm = pyslim.decode_mutation(md)
                with self.assertWarns(FutureWarning):
                    edm = pyslim.encode_mutation(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            with self.assertWarns(FutureWarning):
                pyslim.annotate_mutation_metadata(new_tables, metadata)
            self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_nodes(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.nodes.metadata,
                                           tables.nodes.metadata_offset):
                with self.assertWarns(FutureWarning):
                    dm = pyslim.decode_node(md)
                with self.assertWarns(FutureWarning):
                    edm = pyslim.encode_node(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            with self.assertWarns(FutureWarning):
                pyslim.annotate_node_metadata(new_tables, metadata)
            self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_individuals(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.individuals.metadata,
                                           tables.individuals.metadata_offset):
                with self.assertWarns(FutureWarning):
                    dm = pyslim.decode_individual(md)
                with self.assertWarns(FutureWarning):
                    edm = pyslim.encode_individual(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            with self.assertWarns(FutureWarning):
                pyslim.annotate_individual_metadata(new_tables, metadata)
            self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_populations(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            new_tables = ts.tables
            metadata = []
            for md in tskit.unpack_bytes(tables.populations.metadata,
                                         tables.populations.metadata_offset):
                with self.assertWarns(FutureWarning):
                    dm = pyslim.decode_population(md)
                with self.assertWarns(FutureWarning):
                    edm = pyslim.encode_population(dm)
                self.assertEqual(md, edm)
                metadata.append(dm)

            with self.assertWarns(FutureWarning):
                pyslim.annotate_population_metadata(new_tables, metadata)
            self.assertTableCollectionsEqual(tables, new_tables)


class TestDumpLoad(LegacyPyslimTestCase):
    '''
    Test reading and writing.
    '''

    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="pyslim_ts_")
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
        ts2 = pyslim.load(self.temp_file, legacy_metadata=True)
        self.assertEqual(ts.num_samples, ts2.num_samples)
        self.assertEqual(ts.sequence_length, ts2.sequence_length)
        self.assertEqual(ts.tables, ts2.tables)
        self.assertEqual(ts.reference_sequence, ts2.reference_sequence)

    def test_load_tables(self):
        for ts in self.get_slim_examples():
            self.assertTrue(type(ts) is pyslim.SlimTreeSequence)
            tables = ts.tables
            new_ts = pyslim.load_tables(tables, legacy_metadata=True)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            new_tables = new_ts.tables
            self.assertEqual(tables, new_tables)

    def test_load(self):
        for _, ex in self.get_slim_examples(return_info=True):
            fn = ex['basename'] + ".trees"
            # load in msprime then switch
            msp_ts = tskit.load(fn)
            self.assertTrue(type(msp_ts) is msprime.TreeSequence)
            # transfer tables
            msp_tables = msp_ts.tables
            new_ts = pyslim.load_tables(msp_tables, legacy_metadata=True)
            self.assertTrue(isinstance(new_ts, pyslim.SlimTreeSequence))
            self.verify_times(msp_ts, new_ts)
            new_tables = new_ts.tables
            self.assertTableCollectionsEqual(msp_tables, new_tables)
            # convert directly
            new_ts = pyslim.SlimTreeSequence(msp_ts)
            self.assertTrue(type(new_ts) is pyslim.SlimTreeSequence)
            self.verify_times(msp_ts, new_ts)
            new_tables = new_ts.tables
            self.assertTableCollectionsEqual(msp_tables, new_tables)
            # load to pyslim from file
            slim_ts = pyslim.load(fn, legacy_metadata=True)
            self.assertTrue(type(slim_ts) is pyslim.SlimTreeSequence)
            slim_tables = slim_ts.tables
            self.assertTableCollectionsEqual(msp_tables, slim_tables)
            self.assertEqual(slim_ts.slim_generation, new_ts.slim_generation)

    def test_dump_equality(self):
        for ts in self.get_slim_examples():
            self.verify_dump_equality(ts)


class TestIndividualMetadata(LegacyPyslimTestCase):
    # Tests for extra stuff related to Individuals.

    def test_individual_derived_info(self):
        for ts in self.get_slim_examples():
            for j, ind in enumerate(ts.individuals()):
                a = ts.tables.individuals.metadata_offset[j]
                b = ts.tables.individuals.metadata_offset[j+1]
                raw_md = ts.tables.individuals.metadata[a:b]
                with self.assertWarns(FutureWarning):
                    md = pyslim.decode_individual(raw_md)
                self.assertEqual(ind.metadata, md)
                self.assertEqual(ts.individual(j).metadata, md)
                for n in ind.nodes:
                    self.assertEqual(ts.node(n).population, ind.population)
                    self.assertEqual(ts.node(n).time, ind.time)


class TestNodeMetadata(LegacyPyslimTestCase):
    '''
    Tests for extra stuff related to Nodes.
    '''

    def test_node_derived_info(self):
        for ts in self.get_slim_examples():
            for j, node in enumerate(ts.nodes()):
                a = ts.tables.nodes.metadata_offset[j]
                b = ts.tables.nodes.metadata_offset[j+1]
                raw_md = ts.tables.nodes.metadata[a:b]
                with self.assertWarns(FutureWarning):
                    md = pyslim.decode_node(raw_md)
                self.assertEqual(node.metadata, md)
                self.assertEqual(ts.node(j).metadata, md)


class TestMutationMetadata(LegacyPyslimTestCase):
    '''
    Tests for extra stuff related to Mutations.
    '''

    def test_mutation_derived_info(self):
        for ts in self.get_slim_examples():
            for j, mut in enumerate(ts.mutations()):
                a = ts.tables.mutations.metadata_offset[j]
                b = ts.tables.mutations.metadata_offset[j+1]
                raw_md = ts.tables.mutations.metadata[a:b]
                with self.assertWarns(FutureWarning):
                    md = pyslim.decode_mutation(raw_md)
                self.assertEqual(mut.metadata, md)
                self.assertEqual(ts.mutation(j).metadata, md)


class TestPopulationMetadata(LegacyPyslimTestCase):
    '''
    Tests for extra stuff related to Populations.
    '''

    def test_population_derived_info(self):
        for ts in self.get_slim_examples():
            for j, pop in enumerate(ts.populations()):
                a = ts.tables.populations.metadata_offset[j]
                b = ts.tables.populations.metadata_offset[j+1]
                raw_md = ts.tables.populations.metadata[a:b]
                with self.assertWarns(FutureWarning):
                    md = pyslim.decode_population(raw_md)
                self.assertEqual(pop.metadata, md)
                self.assertEqual(ts.population(j).metadata, md)
