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

    def validate_table_metadata(self, table):
        ms = table.metadata_schema
        for j, row in enumerate(table):
            a = table.metadata_offset[j]
            b = table.metadata_offset[j+1]
            raw_md = table.metadata[a:b]
            # this checks to make sure metadata follows the schema
            enc_md = ms.validate_and_encode_row(row.metadata)
            self.assertEqual(bytes(raw_md), enc_md)

    def validate_metadata(self, ts):
        tables = ts.tables
        for t in (tables.populations, tables.individuals, tables.nodes, tables.edges,
                  tables.sites, tables.mutations, tables.migrations):
            self.validate_table_metadata(t)

    def test_slim_metadata(self):
        for ts in self.get_slim_examples():
            self.validate_metadata(ts)

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


class TestDumpLoad(tests.PyslimTestCase):
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
            for u in mut.metadata['mutation_list']:
                self.assertGreaterEqual(u["nucleotide"], -1)
                self.assertLessEqual(u["nucleotide"], 3)

    def test_nucleotides(self):
        for ts in self.get_slim_examples():
            self.check_nucleotides(ts)


class TestDecoding(tests.PyslimTestCase):
    '''
    Test by comparing decoding to our previous direct implementation of struct decoding.
    '''

    def verify_decoding(self, t, decoder):
        ms = tskit.MetadataSchema(None)
        nt = t.copy()
        nt.metadata_schema = ms
        for a, b in zip(t, nt):
            md = a.metadata
            with self.assertWarns(FutureWarning):
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
            with self.assertWarns(FutureWarning):
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


class TestMetadataAttributeError(tests.PyslimTestCase):

    def test_population_error(self):
        for ts in self.get_slim_examples():
            for x in ts.populations():
                if x.metadata is not None:
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.slim_id
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.selfing_fraction
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.sex_ratio
                with self.assertRaisesRegex(AttributeError, "has no attribute 'ping'$"):
                    _ = x.metadata.ping
                break
            break

    def test_individual_error(self):
        for ts in self.get_slim_examples():
            for x in ts.individuals():
                with self.assertRaisesRegex(AttributeError, 'legacy'):
                    _ = x.metadata.pedigree_id
                with self.assertRaisesRegex(AttributeError, 'legacy'):
                    _ = x.metadata.age
                with self.assertRaisesRegex(AttributeError, 'legacy'):
                    _ = x.metadata.subpopulation
                with self.assertRaisesRegex(AttributeError, 'legacy'):
                    _ = x.metadata.sex
                with self.assertRaisesRegex(AttributeError, 'legacy'):
                    _ = x.metadata.flags
                with self.assertRaisesRegex(AttributeError, "has no attribute 'pong'$"):
                    _ = x.metadata.pong
                break
            break

    def test_node_error(self):
        for ts in self.get_slim_examples():
            for x in ts.nodes():
                if x.metadata is not None:
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.slim_id
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.is_null
                    with self.assertRaisesRegex(AttributeError, 'legacy'):
                        _ = x.metadata.genome_type
                with self.assertRaisesRegex(AttributeError, "has no attribute 'pang'$"):
                    _ = x.metadata.pang
                break
            break

    def test_mutation_error(self):
        for ts in self.get_slim_examples():
            for x in ts.mutations():
                with self.assertRaisesRegex(KeyError, 'legacy'):
                    _ = x.metadata[0]
                with self.assertRaisesRegex(KeyError, 'legacy'):
                    _ = x.metadata[999]
                with self.assertRaises(KeyError):
                    _ = x.metadata['ping']
                break
            break
