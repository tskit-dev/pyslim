"""
Test cases for the metadata reading/writing of pyslim.
"""

import sys

import numpy as np
import pytest
import tskit

import pyslim
import tests

from .recipe_specs import recipe_eq


def assert_nan_equal(a, b):
    if isinstance(a, float):
        assert (a == b) or (np.isnan(a) and np.isnan(b))
    elif isinstance(a, dict) and isinstance(b, dict):
        assert a.keys() == b.keys()
        for k in a:
            assert_nan_equal(a[k], b[k])
    elif isinstance(a, list) and isinstance(b, list):
        for x, y in zip(a, b, strict=True):
            assert_nan_equal(x, y)
    else:
        assert a == b


class TestMetadataSchemas(tests.PyslimTestCase):
    def validate_table_metadata(self, table):
        ms = table.metadata_schema
        for j, row in enumerate(table):
            a = table.metadata_offset[j]
            b = table.metadata_offset[j + 1]
            raw_md = table.metadata[a:b]
            # this checks to make sure metadata follows the schema
            enc_md = ms.validate_and_encode_row(row.metadata)
            assert bytes(raw_md) == enc_md

    def test_slim_metadata(self, recipe):
        for ts in recipe["ts"].values():
            tables = ts.dump_tables()
            for t in (
                tables.populations,
                tables.individuals,
                tables.nodes,
                tables.edges,
                tables.sites,
                tables.mutations,
                tables.migrations,
            ):
                self.validate_table_metadata(t)

    def test_default_metadata_errors(self):
        with pytest.raises(ValueError, match="Unknown metadata request"):
            _ = pyslim.default_slim_metadata("xxx")

    def test_default_metadata(self):
        for k in pyslim.slim_metadata_schemas:
            schema = pyslim.slim_metadata_schemas[k]
            entry = pyslim.default_slim_metadata(k)
            sd = schema.asdict()
            if k != "tree_sequence":
                if sd is not None:
                    for p in sd["properties"]:
                        assert p in entry
                encoded = schema.validate_and_encode_row(entry)
                decoded = schema.decode_row(encoded)
                if entry is None:
                    assert decoded is None
                else:
                    # some defaults have nans, which are not equal
                    assert_nan_equal(entry, decoded)
            else:
                assert k == "tree_sequence"
                for p in sd["json"]["properties"]:
                    assert p in entry
                encoded = schema.validate_and_encode_row(entry)
                decoded = schema.decode_row(encoded)
                assert entry == decoded
                entry["SLiM_mutation_list"].append(
                    pyslim.default_slim_metadata("mutation_list_entry")
                )
                encoded = schema.validate_and_encode_row(entry)
                decoded = schema.decode_row(encoded)
                assert entry == decoded
                entry["SLiM_mutation_list"].append(
                    pyslim.default_slim_metadata("mutation_list_entry")
                )
                encoded = schema.validate_and_encode_row(entry)
                decoded = schema.decode_row(encoded)
                assert entry == decoded

    @pytest.mark.parametrize("num_traits", [1, 5])
    def test_ts_metadata(self, num_traits):
        schema = pyslim.slim_tree_sequence_metadata_schema(num_traits=num_traits)
        entry = pyslim.default_slim_metadata("tree_sequence", num_traits=num_traits)
        assert len(entry["SLiM"]["traits"]) == num_traits
        trait_names = []
        for j, x in enumerate(entry["SLiM"]["traits"]):
            assert x["name"] not in trait_names
            trait_names.append(x["name"])
            assert x["index"] == j
        encoded = schema.validate_and_encode_row(entry)
        decoded = schema.decode_row(encoded)
        assert entry == decoded
        for _ in range(5):
            entry["SLiM_mutation_list"].append(
                pyslim.default_slim_metadata(
                    "mutation_list_entry", num_traits=num_traits
                )
            )
            encoded = schema.validate_and_encode_row(entry)
            decoded = schema.decode_row(encoded)
        for x in entry["SLiM_mutation_list"]:
            assert len(x["per_trait"]) == num_traits

    @pytest.mark.parametrize("num_traits", [1, 5])
    def test_ind_metadata(self, num_traits):
        schema = pyslim.slim_individual_metadata_schema(num_traits=num_traits)
        entry = pyslim.default_slim_metadata("individual", num_traits=num_traits)
        assert len(entry["per_trait"]) == num_traits
        encoded = schema.validate_and_encode_row(entry)
        decoded = schema.decode_row(encoded)
        assert_nan_equal(entry, decoded)

    @pytest.mark.parametrize("num_chroms", [1, 3, 50])
    def test_node_metadata(self, num_chroms):
        schema = pyslim.slim_node_metadata_schema(num_chromosomes=num_chroms)
        entry = pyslim.default_slim_metadata("node", num_chromosomes=num_chroms)
        assert len(entry["is_vacant"]) == pyslim.is_vacant_num_bytes(num_chroms)
        encoded = schema.validate_and_encode_row(entry)
        decoded = schema.decode_row(encoded)
        assert entry == decoded

    @pytest.mark.skipif(
        sys.platform.startswith("win"),
        reason="failing because of dict and OrderedDict comparison?",
    )
    @pytest.mark.parametrize("recipe", recipe_eq("minimal"), indirect=True)
    def test_slim_metadata_schema_equality(self, recipe):
        num_chromosomes = len(recipe["ts"])
        for ts in recipe["ts"].values():
            t = ts.dump_tables()
            num_traits = len(t.metadata["SLiM"]["traits"])
            ts_schema = pyslim.slim_metadata_schemas["tree_sequence"].asdict()
            ts_schema["struct"]["properties"]["SLiM_mutation_list"]["items"][
                "properties"
            ]["per_trait"]["length"] = num_traits
            assert t.metadata_schema.asdict() == ts_schema
            assert t.edges.metadata_schema == pyslim.slim_metadata_schemas["edge"]
            assert t.sites.metadata_schema == pyslim.slim_metadata_schemas["site"]
            assert (
                t.mutations.metadata_schema == pyslim.slim_metadata_schemas["mutation"]
            )
            node_schema = pyslim.slim_metadata_schemas["node"].asdict()
            node_schema["properties"]["is_vacant"]["length"] = int(
                (num_chromosomes + 7) / 8
            )
            assert t.nodes.metadata_schema.asdict() == node_schema
            ind_schema = pyslim.slim_metadata_schemas["individual"].asdict()
            ind_schema["properties"]["per_trait"]["length"] = num_traits
            assert t.individuals.metadata_schema.asdict() == ind_schema
            assert (
                t.populations.metadata_schema
                == pyslim.slim_metadata_schemas["population"]
            )

    def test_node_schema(self):
        s = pyslim.slim_node_metadata_schema()
        sd = s.asdict()
        assert sd["properties"]["is_vacant"]["length"] == 1
        s1 = pyslim.slim_node_metadata_schema(num_chromosomes=1)
        assert s == s1
        s8 = pyslim.slim_node_metadata_schema(num_chromosomes=8)
        assert s == s8
        for nc in [7, 12, 25, 32]:
            sx = pyslim.slim_node_metadata_schema(num_chromosomes=nc)
            sxd = sx.asdict()
            n = sxd["properties"]["is_vacant"]["length"]
            assert n - 1 < nc / 8 and nc / 8 <= n
            sxd["properties"]["is_vacant"]["length"] = 1
            assert sxd == sd


class TestTreeSequenceMetadata(tests.PyslimTestCase):
    arbitrary_recipe = [next(recipe_eq())]  # for testing any one recipe

    def validate_slim_metadata(self, t):
        # t could be tables or a tree sequence
        schema = t.metadata_schema.schema
        assert schema["codec"] == "json+struct"
        assert "SLiM" in schema["json"]["properties"]
        tmd = t.metadata
        assert "SLiM" in tmd
        for k in pyslim.default_slim_metadata("tree_sequence")["SLiM"]:
            assert k in schema["json"]["properties"]["SLiM"]["properties"]
            assert k in tmd["SLiM"]
        sml = schema["struct"]["properties"]
        assert "SLiM_mutation_list" in sml
        for k in pyslim.default_slim_metadata("mutation_list_entry"):
            assert k in sml["SLiM_mutation_list"]["items"]["properties"]

    def validate_model_type(self, tsdict, model_type):
        for _, ts in tsdict.items():
            md = ts.metadata
            assert md["SLiM"]["file_version"] == pyslim.slim_file_version
            assert md["SLiM"]["model_type"] == model_type
            assert md["SLiM"]["tick"] > 0
            assert md["SLiM"]["tick"] >= np.max(ts.tables.nodes.time)

    @pytest.mark.parametrize("recipe", arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata_errors(self, recipe):
        for _, ts in recipe["ts"].items():
            tables = ts.dump_tables()
            tables.metadata_schema = tskit.MetadataSchema(None)
            assert len(tables.metadata) > 0
            with pytest.raises(ValueError):
                pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)

    @pytest.mark.parametrize("recipe", arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata_keeps(self, recipe):
        # make sure doesn't overwrite other stuff
        ts = list(recipe["ts"].values())[0]
        for x in [{}, {"properties": {"abc": {"type": "string"}}}]:
            schema_dict = {
                "codec": "json",
                "type": "object",
            }
            schema_dict.update(x)
            dummy_schema = tskit.MetadataSchema(schema_dict)
            dummy_metadata = {"abc": "foo"}
            tables = ts.dump_tables()
            tables.metadata_schema = dummy_schema
            tables.metadata = dummy_metadata
            pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)
            schema = tables.metadata_schema.schema
            tmd = tables.metadata
            for k in dummy_metadata:
                if len(x) > 0:
                    assert k in schema["json"]["properties"]
                assert k in tmd
                assert tmd[k] == dummy_metadata[k]
            self.validate_slim_metadata(tables)
            assert tmd["SLiM"]["model_type"] == "nonWF"
            assert tmd["SLiM"]["tick"] == 0

    @pytest.mark.parametrize("recipe", arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata_keeps_struct(self, recipe):
        # make sure doesn't overwrite other stuff
        ts = list(recipe["ts"].values())[0]
        json_props = {"num": {"type": "number"}}
        struct_props = {"binnum": {"type": "integer", "binaryFormat": "i", "default": 0}}
        schema_dict = {
            "codec": "json+struct",
            "type": "object",
            "json": {"codec": "json", "type": "object", "properties": json_props},
            "struct": {"struct": "json", "type": "object", "properties": struct_props},
        }
        dummy_schema = tskit.MetadataSchema(schema_dict)
        dummy_metadata = {"num": 12, "binnum": 42}
        tables = ts.dump_tables()
        tables.metadata_schema = dummy_schema
        tables.metadata = dummy_metadata
        pyslim.set_tree_sequence_metadata(tables, "nonWF", 0)
        schema = tables.metadata_schema.schema
        tmd = tables.metadata
        for k in dummy_metadata:
            if k in json_props:
                assert schema["json"]["properties"][k] == json_props[k]
            else:
                assert schema["struct"]["properties"][k] == struct_props[k]
            assert k in tmd
            assert tmd[k] == dummy_metadata[k]
        self.validate_slim_metadata(tables)
        assert tmd["SLiM"]["model_type"] == "nonWF"
        assert tmd["SLiM"]["tick"] == 0

    @pytest.mark.parametrize("recipe", arbitrary_recipe, indirect=True)
    def test_set_tree_sequence_metadata(self, recipe):
        ts = list(recipe["ts"].values())[0]
        tables = ts.dump_tables()
        chroms = [
            {"id": 1, "name": "autosome_1", "symbol": "1", "type": "A", "index": 0},
            {"id": 35, "name": "mtDNA", "symbol": "MT", "type": "HF", "index": 1},
        ]
        traits = [
            {"index": 0, "name": "theTrait", "type": "additive"},
            {"index": 1, "name": "perfectness", "type": "multiplicative"},
        ]
        pyslim.set_tree_sequence_metadata(
            tables,
            "WF",
            tick=99,
            cycle=40,
            stage="early",
            spatial_dimensionality="xy",
            spatial_periodicity="y",
            separate_sexes=False,
            nucleotide_based=True,
            this_chromosome=chroms[1],
            chromosomes=chroms,
            traits=traits,
        )
        self.validate_slim_metadata(tables)
        tmd = tables.metadata
        assert tmd["SLiM"]["model_type"] == "WF"
        assert tmd["SLiM"]["tick"] == 99
        assert tmd["SLiM"]["cycle"] == 40
        assert tmd["SLiM"]["stage"] == "early"
        assert tmd["SLiM"]["spatial_dimensionality"] == "xy"
        assert tmd["SLiM"]["spatial_periodicity"] == "y"
        assert tmd["SLiM"]["separate_sexes"] == False
        assert tmd["SLiM"]["nucleotide_based"] == True
        assert tmd["SLiM"]["chromosomes"] == chroms
        assert tmd["SLiM"]["this_chromosome"] == chroms[1]
        assert tmd["SLiM"]["traits"] == traits

    @pytest.mark.parametrize("recipe", recipe_eq("WF"), indirect=True)
    def test_WF_model_type(self, recipe):
        self.validate_model_type(recipe["ts"], "WF")

    @pytest.mark.parametrize("recipe", recipe_eq("nonWF"), indirect=True)
    def test_nonWF_model_type(self, recipe):
        self.validate_model_type(recipe["ts"], "nonWF")

    @pytest.mark.parametrize(
        "recipe",
        recipe_eq(exclude=["user_metadata", "multichrom", "record_mutations"]),
        indirect=True,
    )
    def test_recover_metadata(self, recipe):
        # msprime <=0.7.5 discards metadata, but we can recover it from provenance
        for _, ts in recipe["ts"].items():
            tables = ts.dump_tables()
            tables.metadata_schema = tskit.MetadataSchema(None)
            tables.metadata = b""
            pyslim.update_tables(tables)
            md = tables.metadata
            assert "SLiM" in md
            tsmd = ts.metadata["SLiM"]
            for k in tsmd:
                assert k in md["SLiM"]
                # slim does not write out empty descriptions
                if k != "description" or tsmd[k] != "":
                    assert tsmd[k] == md["SLiM"][k]

    @pytest.mark.parametrize(
        "recipe", recipe_eq("recipe_with_metadata.slim"), indirect=True
    )
    def test_user_metadata(self, recipe):
        for _, ts in recipe["ts"].items():
            md = ts.metadata["SLiM"]
            assert "user_metadata" in md
            assert md["user_metadata"] == {"hello": ["world"], "pi": [3, 1, 4, 1, 5, 9]}

    @pytest.mark.parametrize(
        "recipe", recipe_eq("recipe_with_metadata.slim"), indirect=True
    )
    def test_population_names(self, recipe):
        for _, ts in recipe["ts"].items():
            md = ts.metadata["SLiM"]
            assert ts.num_populations == 4
            p = ts.population(1)
            assert p.metadata["name"] == "first_population"
            assert p.metadata["description"] == "i'm the first population"
            p = ts.population(3)
            assert p.metadata["name"] == "other_population"
            assert p.metadata["description"] == "i'm the other population"


class TestAlleles(tests.PyslimTestCase):
    """
    Test nothing got messed up with haplotypes.
    """

    def test_haplotypes(self, recipe):
        for _, slim_ts in recipe["ts"].items():
            tables = slim_ts.dump_tables()
            ts = tables.tree_sequence()
            self.verify_haplotype_equality(ts, slim_ts)


class TestNucleotides(tests.PyslimTestCase):
    """
    Test nucleotide support
    """

    def test_nucleotides(self, recipe):
        """
        Check that nucleotides are all valid, i.e.,
        -1, 0, 1, 2, or 3.
        """
        for _, ts in recipe["ts"].items():
            for u in ts.metadata["SLiM_mutation_list"]:
                assert u["nucleotide"] >= -1
                assert u["nucleotide"] <= 3


class TestMultichrom(tests.PyslimTestCase):
    """
    Test multichromosome metadata
    """

    @pytest.mark.parametrize("recipe", recipe_eq("multichrom"), indirect=True)
    def test_chromosome_types(self, recipe):
        chroms = recipe["ts"]
        chrom_info = {}
        chrom_list = None
        # check 'chromosomes' identical for all
        for chrom, ts in chroms.items():
            md = ts.metadata["SLiM"]
            chrom_info[chrom] = md["this_chromosome"]
            # optional but should be there for all these examples
            assert "chromosomes" in md
            if chrom_list is None:
                chrom_list = md["chromosomes"]
            else:
                assert chrom_list == md["chromosomes"]
        # check 'this_chromosome' matches entry in 'chromosomes'
        chrom_list_d = {f"chromosome_{x['symbol']}": x for x in chrom_list}
        for chrom in chroms.keys():
            assert chrom in chrom_list_d
            assert chrom_list_d[chrom] == chrom_info[chrom]
