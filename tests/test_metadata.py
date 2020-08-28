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
import numpy as np

class TestMetadataSchemas(tests.PyslimTestCase):

    def test_default_metadata(self):
        for k in pyslim.slim_metadata_schemas:
            self.assertTrue(k in pyslim.default_slim_metadata)
            schema = pyslim.slim_metadata_schemas[k]
            entry = pyslim.default_slim_metadata[k]
            encoded = schema.validate_and_encode_row(entry)
            decoded = schema.decode_row(encoded)
            if entry is None:
                self.assertTrue(decoded is None)
            else:
                self.assertDictEqual(entry, decoded)

    def test_slim_metadata_schema_equality(self):
        for ts in self.get_slim_examples():
            t = ts.tables
            self.assertEqual(str(t.metadata_schema),
                             str(pyslim.slim_metadata_schemas['tree_sequence']))
            self.assertEqual(t.edges.metadata_schema,
                             pyslim.slim_metadata_schemas['edge'])
            self.assertEqual(t.sites.metadata_schema,
                             pyslim.slim_metadata_schemas['site'])
            self.assertEqual(
                    t.mutations.metadata_schema,
                    pyslim.slim_metadata_schemas['mutation'])
            self.assertEqual(
                    t.nodes.metadata_schema,
                    pyslim.slim_metadata_schemas['node'])
            self.assertEqual(
                    t.individuals.metadata_schema,
                    pyslim.slim_metadata_schemas['individual'])
            self.assertEqual(
                    t.populations.metadata_schema,
                    pyslim.slim_metadata_schemas['population'])
            break

class TestTreeSequenceMetadata(tests.PyslimTestCase):

    def validate_slim_metadata(self, t):
        # t could be tables or a tree sequence
        schema = t.metadata_schema.schema
        self.assertTrue('SLiM' in schema['properties'])
        self.assertTrue('SLiM' in t.metadata)
        for k in pyslim.default_slim_metadata['tree_sequence']['SLiM']:
            self.assertTrue(k in schema['properties']['SLiM']['properties'])
            self.assertTrue(k in t.metadata['SLiM'])

    def test_set_tree_sequence_metadata_errors(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            tables.metadata_schema = tskit.MetadataSchema(None)
            self.assertGreater(len(tables.metadata), 0)
            with self.assertRaises(ValueError):
                pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)
            break

    def test_set_tree_sequence_metadata_keeps(self):
        # make sure doesn't overwrite other stuff
        dummy_schema = tskit.MetadataSchema({
                'codec': 'json',
                'type': 'object',
                'properties': { 'abc': { 'type': 'string' } }
                })
        dummy_metadata = { 'abc': 'foo' }
        for ts in self.get_slim_examples():
            tables = ts.tables
            tables.metadata_schema = dummy_schema
            tables.metadata = dummy_metadata
            pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)
            schema = tables.metadata_schema.schema
            for k in dummy_metadata:
                self.assertTrue(k in schema['properties'])
                self.assertTrue(k in tables.metadata)
                self.assertEqual(tables.metadata[k], dummy_metadata[k])
            self.validate_slim_metadata(tables)
            self.assertEqual(tables.metadata['SLiM']['model_type'], "nonWF")
            self.assertEqual(tables.metadata['SLiM']['generation'], 0)
            break

    def test_set_tree_sequence_metadata(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            pyslim.set_tree_sequence_metadata(
                    tables, "WF", 99, 
                    spatial_dimensionality='xy',
                    spatial_periodicity='y',
                    separate_sexes=False,
                    nucleotide_based=True)
            self.validate_slim_metadata(tables)
            self.assertEqual(tables.metadata['SLiM']['model_type'], "WF")
            self.assertEqual(tables.metadata['SLiM']['generation'], 99)
            self.assertEqual(tables.metadata['SLiM']['spatial_dimensionality'], 'xy')
            self.assertEqual(tables.metadata['SLiM']['spatial_periodicity'], 'y')
            self.assertEqual(tables.metadata['SLiM']['separate_sexes'], False)
            self.assertEqual(tables.metadata['SLiM']['nucleotide_based'], True)
            break

    def test_model_type(self):
        for model_type in ["WF", "nonWF"]:
            args = {model_type: True}
            for ts in self.get_slim_examples(**args):
                self.assertEqual(ts.metadata['SLiM']['file_version'], pyslim.slim_file_version)
                self.assertEqual(ts.metadata['SLiM']['model_type'], model_type)
                self.assertGreater(ts.metadata['SLiM']['generation'], 0)
                self.assertGreaterEqual(ts.metadata['SLiM']['generation'],
                        np.max(ts.tables.nodes.time))

    def test_recover_metadata(self):
        # msprime <=0.7.5 discards metadata, but we can recover it from provenance
        for ts in self.get_slim_examples():
            t = ts.tables
            t.metadata_schema = tskit.MetadataSchema(None)
            t.metadata = b''
            new_ts = pyslim.load_tables(t)
            self.assertEqual(new_ts.metadata, ts.metadata)


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
            self.assertTableCollectionsEqual(tables, new_tables)


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
        self.assertEqual(ts.reference_sequence, ts2.reference_sequence)

    def test_load_tables(self):
        for ts in self.get_slim_examples():
            self.assertTrue(type(ts) is pyslim.SlimTreeSequence)
            tables = ts.tables
            new_ts = pyslim.load_tables(tables)
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
            new_ts = pyslim.load_tables(msp_tables)
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
            slim_ts = pyslim.load(fn)
            self.assertTrue(type(slim_ts) is pyslim.SlimTreeSequence)
            slim_tables = slim_ts.tables
            self.assertTableCollectionsEqual(msp_tables, slim_tables)
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
            print(mut)
            for u in mut.metadata:
                self.assertGreaterEqual(u.nucleotide, -1)
                self.assertLessEqual(u.nucleotide, 3)

    def test_nucleotides(self):
        for ts in self.get_slim_examples():
            self.check_nucleotides(ts)


class TestLegacyDecoding(tests.PyslimTestCase):
    '''
    Test by comparing decoding to our previous direct implementation of struct decoding.
    '''

    def verify_decoding(self, t, decoder):
        ms = tskit.MetadataSchema(None)
        nt = t.copy()
        nt.metadata_schema = ms
        for a, b in zip(t, nt):
            md = a.metadata
            omd = decoder(b.metadata)
            if md is None:
                self.assertTrue(omd is None)
            else:
                self.assertEqual(md, omd.asdict())

    def verify_mutation_decoding(self, t):
        ms = tskit.MetadataSchema(None)
        nt = t.copy()
        nt.metadata_schema = ms
        for a, b in zip(t, nt):
            md = a.metadata
            omd = pyslim.decode_mutation(b.metadata)
            self.assertEqual(md,
                    {"mutation_list": [u.asdict() for u in omd]})

    def test_decoding(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            self.verify_decoding(tables.populations, pyslim.decode_population)
            self.verify_decoding(tables.individuals, pyslim.decode_individual)
            self.verify_decoding(tables.nodes, pyslim.decode_node)
            self.verify_mutation_decoding(tables.mutations)
