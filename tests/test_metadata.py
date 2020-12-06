"""
Test cases for the metadata reading/writing of pyslim.
"""
from __future__ import print_function
from __future__ import division

import tests
import unittest
import os
import tempfile
import numpy as np
import pytest
import msprime
import tskit
import pyslim


class TestMetadataSchemas(tests.PyslimTestCase):

    def validate_table_metadata(self, table):
        ms = table.metadata_schema
        for j, row in enumerate(table):
            a = table.metadata_offset[j]
            b = table.metadata_offset[j+1]
            raw_md = table.metadata[a:b]
            # this checks to make sure metadata follows the schema
            enc_md = ms.validate_and_encode_row(row.metadata)
            assert bytes(raw_md) == enc_md

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
            schema = pyslim.slim_metadata_schemas[k]
            entry = pyslim.default_slim_metadata(k)
            encoded = schema.validate_and_encode_row(entry)
            decoded = schema.decode_row(encoded)
            if entry is None:
                assert decoded is None
            else:
                assert entry == decoded

    def test_slim_metadata_schema_equality(self):
        for ts in self.get_slim_examples():
            t = ts.tables
            assert str(t.metadata_schema) == str(pyslim.slim_metadata_schemas['tree_sequence'])
            assert t.edges.metadata_schema == pyslim.slim_metadata_schemas['edge']
            assert t.sites.metadata_schema == pyslim.slim_metadata_schemas['site']
            assert t.mutations.metadata_schema == pyslim.slim_metadata_schemas['mutation']
            assert t.nodes.metadata_schema == pyslim.slim_metadata_schemas['node']
            assert t.individuals.metadata_schema == pyslim.slim_metadata_schemas['individual']
            assert t.populations.metadata_schema == pyslim.slim_metadata_schemas['population']


class TestTreeSequenceMetadata(tests.PyslimTestCase):

    def validate_slim_metadata(self, t):
        # t could be tables or a tree sequence
        schema = t.metadata_schema.schema
        assert 'SLiM' in schema['properties']
        assert 'SLiM' in t.metadata
        for k in pyslim.default_slim_metadata('tree_sequence')['SLiM']:
            assert k in schema['properties']['SLiM']['properties']
            assert k in t.metadata['SLiM']

    def test_set_tree_sequence_metadata_errors(self):
        for ts in self.get_slim_examples():
            tables = ts.tables
            tables.metadata_schema = tskit.MetadataSchema(None)
            assert len(tables.metadata) > 0
            with pytest.raises(ValueError):
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
                assert k in schema['properties']
                assert k in tables.metadata
                assert tables.metadata[k] == dummy_metadata[k]
            self.validate_slim_metadata(tables)
            assert tables.metadata['SLiM']['model_type'] == "nonWF"
            assert tables.metadata['SLiM']['generation'] == 0
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
            assert tables.metadata['SLiM']['model_type'] == "WF"
            assert tables.metadata['SLiM']['generation'] == 99
            assert tables.metadata['SLiM']['spatial_dimensionality'] == 'xy'
            assert tables.metadata['SLiM']['spatial_periodicity'] == 'y'
            assert tables.metadata['SLiM']['separate_sexes'] == False
            assert tables.metadata['SLiM']['nucleotide_based'] == True
            break

    def test_model_type(self):
        for model_type in ["WF", "nonWF"]:
            args = {model_type: True}
            for ts in self.get_slim_examples(**args):
                assert ts.metadata['SLiM']['file_version'] == pyslim.slim_file_version
                assert ts.metadata['SLiM']['model_type'] == model_type
                assert ts.metadata['SLiM']['generation'] > 0
                assert ts.metadata['SLiM']['generation'] >= np.max(ts.tables.nodes.time)

    def test_recover_metadata(self):
        # msprime <=0.7.5 discards metadata, but we can recover it from provenance
        for ts in self.get_slim_examples():
            t = ts.tables
            t.metadata_schema = tskit.MetadataSchema(None)
            t.metadata = b''
            new_ts = pyslim.load_tables(t)
            assert new_ts.metadata == ts.metadata


class TestDumpLoad(tests.PyslimTestCase):
    '''
    Test reading and writing.
    '''

    def setup_class(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="pyslim_ts_")
        os.close(fd)

    def teardown_class(self):
        os.unlink(self.temp_file)

    def verify_times(self, ts, slim_ts):
        gen = slim_ts.slim_generation
        assert ts.num_nodes == slim_ts.num_nodes
        # verify internal consistency
        for j in range(slim_ts.num_nodes):
            assert slim_ts.node(j).time == slim_ts.tables.nodes.time[j]
        # verify consistency between tree sequences
        for n1, n2 in zip(ts.nodes(), slim_ts.nodes()):
            assert n1.time == n2.time

    def verify_dump_equality(self, ts):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(self.temp_file)
        ts2 = pyslim.load(self.temp_file)
        assert ts.num_samples == ts2.num_samples
        assert ts.sequence_length == ts2.sequence_length
        assert ts.tables == ts2.tables
        assert ts.reference_sequence == ts2.reference_sequence

    def test_load_tables(self):
        for ts in self.get_slim_examples():
            assert isinstance(ts, pyslim.SlimTreeSequence)
            tables = ts.tables
            new_ts = pyslim.load_tables(tables)
            assert isinstance(new_ts, pyslim.SlimTreeSequence)
            new_tables = new_ts.tables
            assert tables == new_tables

    def test_load(self):
        for _, ex in self.get_slim_examples(return_info=True):
            fn = ex['basename'] + ".trees"
            # load in msprime then switch
            msp_ts = tskit.load(fn)
            assert isinstance(msp_ts, tskit.TreeSequence)
            # transfer tables
            msp_tables = msp_ts.tables
            new_ts = pyslim.load_tables(msp_tables)
            assert isinstance(new_ts, pyslim.SlimTreeSequence)
            self.verify_times(msp_ts, new_ts)
            new_tables = new_ts.tables
            self.assertTableCollectionsEqual(msp_tables, new_tables)
            # convert directly
            new_ts = pyslim.SlimTreeSequence(msp_ts)
            assert isinstance(new_ts, pyslim.SlimTreeSequence)
            self.verify_times(msp_ts, new_ts)
            new_tables = new_ts.tables
            self.assertTableCollectionsEqual(msp_tables, new_tables)
            # load to pyslim from file
            slim_ts = pyslim.load(fn)
            assert isinstance(slim_ts, pyslim.SlimTreeSequence)
            slim_tables = slim_ts.tables
            self.assertTableCollectionsEqual(msp_tables, slim_tables)
            assert slim_ts.slim_generation == new_ts.slim_generation

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
                assert u["nucleotide"] >= -1
                assert u["nucleotide"] <= 3

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
            with pytest.warns(FutureWarning):
                omd = decoder(b.metadata)
            if md is None:
                assert omd is None
            else:
                assert md == omd.asdict()

    def verify_mutation_decoding(self, t):
        ms = tskit.MetadataSchema(None)
        nt = t.copy()
        nt.metadata_schema = ms
        for a, b in zip(t, nt):
            md = a.metadata
            with pytest.warns(FutureWarning):
                omd = pyslim.decode_mutation(b.metadata)
            assert md == {"mutation_list": [u.asdict() for u in omd]}

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
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.slim_id
                    assert 'legacy' in str(exec_info)
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.selfing_fraction
                    assert 'legacy' in str(exec_info)
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.sex_ratio
                    assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.ping
                assert "has no attribute 'ping'" in str(exec_info)
                break
            break

    def test_individual_error(self):
        for ts in self.get_slim_examples():
            for x in ts.individuals():
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.pedigree_id
                assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.age
                assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.subpopulation
                assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.sex
                assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.flags
                assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.pong
                assert "has no attribute 'pong'" in str(exec_info)
                break
            break

    def test_node_error(self):
        for ts in self.get_slim_examples():
            for x in ts.nodes():
                if x.metadata is not None:
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.slim_id
                    assert 'legacy' in str(exec_info)
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.is_null
                    assert 'legacy' in str(exec_info)
                    with pytest.raises(AttributeError) as exec_info:
                        _ = x.metadata.genome_type
                    assert 'legacy' in str(exec_info)
                with pytest.raises(AttributeError) as exec_info:
                    _ = x.metadata.pang
                assert "has no attribute 'pang'" in str(exec_info)
                break
            break

    def test_mutation_error(self):
        for ts in self.get_slim_examples():
            for x in ts.mutations():
                with pytest.raises(KeyError) as exec_info:
                    _ = x.metadata[0]
                assert 'legacy' in str(exec_info)
                with pytest.raises(KeyError) as exec_info:
                    _ = x.metadata[999]
                assert 'legacy' in str(exec_info)
                with pytest.raises(KeyError):
                    _ = x.metadata['ping']
                break
            break
