"""
Test cases for the metadata reading/writing of pyslim.
"""
import os

import numpy as np
import pytest
import msprime
import tskit
import pyslim

import tests
from .recipe_specs import recipe_eq


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

    def test_slim_metadata(self, recipe):
        for ts in recipe["ts"].values():
            tables = ts.dump_tables()
            for t in (tables.populations, tables.individuals, tables.nodes, tables.edges,
                      tables.sites, tables.mutations, tables.migrations):
                self.validate_table_metadata(t)

    def test_default_metadata_errors(self):
        with pytest.raises(ValueError, match="Unknown metadata request"):
            _ = pyslim.default_slim_metadata("xxx")

    def test_default_metadata(self):
        for k in pyslim.slim_metadata_schemas:
            schema = pyslim.slim_metadata_schemas[k]
            entry = pyslim.default_slim_metadata(k)
            sd = schema.asdict()
            if sd is not None:
                for p in sd['properties']:
                    assert p in entry
            encoded = schema.validate_and_encode_row(entry)
            decoded = schema.decode_row(encoded)
            if entry is None:
                assert decoded is None
            else:
                assert entry == decoded
        schema = pyslim.slim_metadata_schemas["mutation"]
        entry = pyslim.default_slim_metadata("mutation")
        entry['mutation_list'].append(
            pyslim.default_slim_metadata("mutation_list_entry")
        )
        encoded = schema.validate_and_encode_row(entry)
        decoded = schema.decode_row(encoded)
        assert entry == decoded
        entry['mutation_list'].append(
            pyslim.default_slim_metadata("mutation_list_entry")
        )
        encoded = schema.validate_and_encode_row(entry)
        decoded = schema.decode_row(encoded)
        assert entry == decoded

    def test_slim_metadata_schema_equality(self, recipe):
        num_chromosomes = len(recipe["ts"])
        for ts in recipe["ts"].values():
            t = ts.dump_tables()
            assert t.metadata_schema == pyslim.slim_metadata_schemas['tree_sequence']
            assert t.edges.metadata_schema == pyslim.slim_metadata_schemas['edge']
            assert t.sites.metadata_schema == pyslim.slim_metadata_schemas['site']
            assert t.mutations.metadata_schema == pyslim.slim_metadata_schemas['mutation']
            node_schema = pyslim.slim_metadata_schemas['node'].asdict()
            node_schema['properties']['is_vacant']['length'] = int((num_chromosomes + 7)/8)
            assert t.nodes.metadata_schema.asdict() == node_schema
            assert t.individuals.metadata_schema == pyslim.slim_metadata_schemas['individual']
            assert t.populations.metadata_schema == pyslim.slim_metadata_schemas['population']

    def test_node_schema(self):
        s = pyslim.slim_node_metadata_schema()
        sd = s.asdict()
        assert sd['properties']['is_vacant']['length'] == 1
        s1 = pyslim.slim_node_metadata_schema(num_chromosomes=1)
        assert s == s1
        s8 = pyslim.slim_node_metadata_schema(num_chromosomes=8)
        assert s == s8
        for nc in [7, 12, 25, 32]:
            sx = pyslim.slim_node_metadata_schema(num_chromosomes=nc)
            sxd = sx.asdict()
            n = sxd['properties']['is_vacant']['length']
            assert n - 1 < nc / 8 and nc / 8 <= n
            sxd['properties']['is_vacant']['length']  = 1
            assert sxd == sd


class TestTreeSequenceMetadata(tests.PyslimTestCase):
    arbitrary_recipe = [next(recipe_eq())]  # for testing any one recipe

    def validate_slim_metadata(self, t):
        # t could be tables or a tree sequence
        schema = t.metadata_schema.schema
        assert 'SLiM' in schema['properties']
        assert 'SLiM' in t.metadata
        for k in pyslim.default_slim_metadata('tree_sequence')['SLiM']:
            assert k in schema['properties']['SLiM']['properties']
            assert k in t.metadata['SLiM']

    def validate_model_type(self, tsdict, model_type):
        for _, ts in tsdict.items():
            assert ts.metadata['SLiM']['file_version'] == pyslim.slim_file_version
            assert ts.metadata['SLiM']['model_type'] == model_type
            assert ts.metadata['SLiM']['tick'] > 0
            assert ts.metadata['SLiM']['tick'] >= np.max(ts.tables.nodes.time)

    @pytest.mark.parametrize('recipe', arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata_errors(self, recipe):
        for _, ts in recipe["ts"].items():
            tables = ts.dump_tables()
            tables.metadata_schema = tskit.MetadataSchema(None)
            assert len(tables.metadata) > 0
            with pytest.raises(ValueError):
                pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)

    @pytest.mark.parametrize('recipe', arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata_keeps(self, recipe):
        # make sure doesn't overwrite other stuff
        ts = list(recipe["ts"].values())[0]
        for x in [{}, { 'properties': { 'abc': { 'type': 'string' } } }]:
            schema_dict = {
                    'codec': 'json',
                    'type': 'object',
            }
            schema_dict.update(x)
            dummy_schema = tskit.MetadataSchema(schema_dict)
            dummy_metadata = { 'abc': 'foo' }
            tables = ts.dump_tables()
            tables.metadata_schema = dummy_schema
            tables.metadata = dummy_metadata
            pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)
            schema = tables.metadata_schema.schema
            for k in dummy_metadata:
                if len(x) > 0:
                    assert k in schema['properties']
                assert k in tables.metadata
                assert tables.metadata[k] == dummy_metadata[k]
            self.validate_slim_metadata(tables)
            assert tables.metadata['SLiM']['model_type'] == "nonWF"
            assert tables.metadata['SLiM']['tick'] == 0

    @pytest.mark.parametrize('recipe', arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata(self, recipe):
        ts = list(recipe["ts"].values())[0]
        tables = ts.dump_tables()
        pyslim.set_tree_sequence_metadata(
                tables,
                "WF",
                tick=99,
                cycle=40,
                stage="early",
                spatial_dimensionality='xy',
                spatial_periodicity='y',
                separate_sexes=False,
                nucleotide_based=True
        )
        self.validate_slim_metadata(tables)
        assert tables.metadata['SLiM']['model_type'] == "WF"
        assert tables.metadata['SLiM']['tick'] == 99
        assert tables.metadata['SLiM']['cycle'] == 40
        assert tables.metadata['SLiM']['stage'] == 'early'
        assert tables.metadata['SLiM']['spatial_dimensionality'] == 'xy'
        assert tables.metadata['SLiM']['spatial_periodicity'] == 'y'
        assert tables.metadata['SLiM']['separate_sexes'] == False
        assert tables.metadata['SLiM']['nucleotide_based'] == True

    @pytest.mark.parametrize('recipe', recipe_eq("WF"), indirect=True)
    def test_WF_model_type(self, recipe):
        self.validate_model_type(recipe["ts"], "WF")

    @pytest.mark.parametrize('recipe', recipe_eq("nonWF"), indirect=True)
    def test_nonWF_model_type(self, recipe):
        self.validate_model_type(recipe["ts"], "nonWF")

    @pytest.mark.parametrize(
        'recipe', recipe_eq(exclude=["user_metadata", "multichrom"]), indirect=True)
    def test_recover_metadata(self, recipe):
        # msprime <=0.7.5 discards metadata, but we can recover it from provenance
        # HOWEVER: multichromosome information is not saved
        # but this is not something we need to maintain any more.
        for _, ts in recipe["ts"].items():
            tables = ts.dump_tables()
            tables.metadata_schema = tskit.MetadataSchema(None)
            tables.metadata = b''
            pyslim.update_tables(tables)
            md = tables.metadata
            assert 'SLiM' in md
            for k in ts.metadata['SLiM']:
                if k in ("chromosomes", "this_chromosome"):
                    continue  # TODO: see https://github.com/MesserLab/SLiM/issues/520
                assert k in md['SLiM']
                # slim does not write out empty descriptions
                if k != 'description' or ts.metadata['SLiM'][k] != "":
                    assert ts.metadata['SLiM'][k] == md['SLiM'][k]

    @pytest.mark.parametrize('recipe', recipe_eq("recipe_with_metadata.slim"), indirect=True)
    def test_user_metadata(self, recipe):
        for _, ts in recipe["ts"].items():
            md = ts.metadata["SLiM"]
            assert "user_metadata" in md
            assert md['user_metadata'] == {
                    "hello" : ["world"],
                    "pi" : [3, 1, 4, 1, 5, 9]
                    }

    @pytest.mark.parametrize('recipe', recipe_eq("recipe_with_metadata.slim"), indirect=True)
    def test_population_names(self, recipe):
        for _, ts in recipe["ts"].items():
            md = ts.metadata["SLiM"]
            assert ts.num_populations == 4
            p = ts.population(1)
            assert p.metadata['name'] == "first_population"
            assert p.metadata['description'] == "i'm the first population"
            p = ts.population(3)
            assert p.metadata['name'] == "other_population"
            assert p.metadata['description'] == "i'm the other population"


class TestAlleles(tests.PyslimTestCase):
    '''
    Test nothing got messed up with haplotypes.
    '''

    def test_haplotypes(self, recipe):
        for _, slim_ts in recipe["ts"].items():
            tables = slim_ts.dump_tables()
            ts = tables.tree_sequence()
            self.verify_haplotype_equality(ts, slim_ts)


class TestNucleotides(tests.PyslimTestCase):
    '''
    Test nucleotide support
    '''

    def test_nucleotides(self, recipe):
        '''
        Check that nucleotides are all valid, i.e.,
        -1, 0, 1, 2, or 3.
        '''
        for _, ts in recipe["ts"].items():
            for mut in ts.mutations():
                for u in mut.metadata['mutation_list']:
                    assert u["nucleotide"] >= -1
                    assert u["nucleotide"] <= 3


class TestMultichrom(tests.PyslimTestCase):
    '''
    Test multichromosome metadata
    '''

    @pytest.mark.parametrize('recipe', recipe_eq("multichrom"), indirect=True)
    def test_chromosome_types(self, recipe):
        chroms = recipe["ts"]
        chrom_info = {}
        chrom_list = None
        # check 'chromosomes' identical for all
        for chrom, ts in chroms.items():
            md = ts.metadata['SLiM']
            chrom_info[chrom] = md['this_chromosome']
            # optional but should be there for all these examples
            assert 'chromosomes' in md
            if chrom_list is None:
                chrom_list = md['chromosomes']
            else:
                assert chrom_list == md['chromosomes']
        # check 'this_chromosome' matches entry in 'chromosomes'
        chrom_list_d = { f"chromosome_{x['symbol']}" : x for x in chrom_list }
        for chrom in chroms.keys():
            assert chrom in chrom_list_d
            assert chrom_list_d[chrom] == chrom_info[chrom]

