"""
Test cases for the deprecated, legacy metadata representation of pyslim.
"""
import os
import warnings

import numpy as np
import pytest
import msprime
import tskit

import pyslim
import tests



class TestLegacyTypes(tests.PyslimTestCase):

    def test_test_warnings(self):
        # check that our checking for warnings works
        # (it didn't with DeprecationWarnings)
        with pytest.warns(FutureWarning):
            warnings.warn('hi', FutureWarning)

    def test_legacy_types(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        assert ts.legacy_metadata
        for pop in ts.populations():
            if pop.metadata is not None:
                assert isinstance(pop.metadata, pyslim.PopulationMetadata)
                break
        for ind in ts.individuals():
            assert isinstance(ind.metadata, pyslim.IndividualMetadata)
            break
        for n in ts.nodes():
            if n.metadata is not None:
                assert isinstance(n.metadata, pyslim.NodeMetadata)
                break
        for mut in ts.mutations():
            assert isinstance(mut.metadata, list)
            if len(mut.metadata) > 0:
                for u in mut.metadata:
                    assert isinstance(u, pyslim.MutationMetadata)
                break


class TestEncodeDecode(tests.PyslimTestCase):
    '''
    Tests for conversion to/from binary representations of metadata.
    '''

    def test_legacy_errors(self):
        defaults = pyslim.default_slim_metadata
        with pytest.raises(ValueError) as exec_info:
            pyslim.decode_mutation(defaults('mutation'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.decode_population(defaults('population'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.decode_individual(defaults('individual'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.decode_node(defaults('node'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.encode_mutation(defaults('mutation'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.encode_population(defaults('population'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.encode_individual(defaults('individual'))
        assert "legacy" in str(exec_info)
        with pytest.raises(ValueError) as exec_info:
            pyslim.encode_node(defaults('node'))
        assert "legacy" in str(exec_info)

    def test_decode_errors(self):
        with pytest.warns(FutureWarning):
            with pytest.raises(ValueError):
                pyslim.decode_mutation(2.0)
            with pytest.raises(ValueError):
                pyslim.decode_mutation([2.0, 3.0])
            with pytest.raises(ValueError):
                pyslim.decode_node(2.0)
            with pytest.raises(ValueError):
                pyslim.decode_node([1, 2])
            with pytest.raises(ValueError):
                pyslim.decode_individual(3.0)
            with pytest.raises(ValueError):
                pyslim.decode_individual([1, 2])
            with pytest.raises(ValueError):
                pyslim.decode_population(1.0)
            with pytest.raises(ValueError):
                pyslim.decode_population([2, 3])

    def test_decode_already_mutation(self):
        m = [pyslim.MutationMetadata(mutation_type = 0,
                                     selection_coeff = 0.2,
                                     population = k,
                                     slim_time = 130,
                                     nucleotide = 2) for k in range(4)]
        with pytest.warns(FutureWarning):
            dm = pyslim.decode_mutation(m)
        assert isinstance(dm, list)
        for a, b in zip(m, dm):
            assert a == b

    def test_decode_already_node(self):
        m = pyslim.NodeMetadata(slim_id=123, is_null=True, genome_type=0)
        with pytest.warns(FutureWarning):
            dm = pyslim.decode_node(m)
        assert m == dm

    def test_decode_already_population(self):
        m = pyslim.PopulationMetadata(slim_id=1, selfing_fraction=0.2,
                                      female_cloning_fraction=0.3,
                                      male_cloning_fraction=0.4,
                                      sex_ratio=0.5, bounds_x0=0, bounds_x1=10,
                                      bounds_y0=2, bounds_y1=20, bounds_z0=3,
                                      bounds_z1=30, migration_records=[])
        with pytest.warns(FutureWarning):
            dm = pyslim.decode_population(m)
        assert m == dm

    def test_decode_already_individual(self):
        m = pyslim.IndividualMetadata(pedigree_id=24, age=8, population=1,
                                      sex=1, flags=0)
        with pytest.warns(FutureWarning):
            dm = pyslim.decode_individual(m)
        assert m == dm

    def test_mutation_metadata(self):
        for md_length in [0, 1, 5]:
            md = [pyslim.MutationMetadata(
                     mutation_type=j, selection_coeff=0.5, population=j,
                     slim_time=10 + j, nucleotide=(j % 5) - 1)
                     for j in range(md_length)]
            with pytest.warns(FutureWarning):
                md_bytes = pyslim.encode_mutation(md)
            with pytest.warns(FutureWarning):
                new_md = pyslim.decode_mutation(md_bytes)
            assert len(md) == len(new_md)
            for x, y in zip(md, new_md):
                assert x == y

    def test_node_metadata(self):
        md = pyslim.NodeMetadata(slim_id=2, is_null=False,
                                 genome_type=pyslim.GENOME_TYPE_X)
        with pytest.warns(FutureWarning):
            md_bytes = pyslim.encode_node(md)
        with pytest.warns(FutureWarning):
            new_md = pyslim.decode_node(md_bytes)
        assert md == new_md

    def test_individual_metadata(self):
        md = pyslim.IndividualMetadata(
                age=2, pedigree_id=23, population=0,
                sex=pyslim.INDIVIDUAL_TYPE_MALE,
                flags=pyslim.INDIVIDUAL_FLAG_MIGRATED)
        with pytest.warns(FutureWarning):
            md_bytes = pyslim.encode_individual(md)
        with pytest.warns(FutureWarning):
            new_md = pyslim.decode_individual(md_bytes)
        assert md == new_md

    def test_population_metadata(self):
        mrs = [pyslim.PopulationMigrationMetadata(source_subpop=j, migration_rate=0.2)
               for j in range(3)]
        for mr_list in [[], mrs]:
            md = pyslim.PopulationMetadata(
                    slim_id=1, selfing_fraction=0.75, female_cloning_fraction=0.2,
                    male_cloning_fraction=0.8, sex_ratio=0.6, bounds_x0=-1.0,
                    bounds_x1=2.0, bounds_y0=0.0, bounds_y1=0.0, bounds_z0=0.0,
                    bounds_z1=1.0, migration_records=mr_list)
            with pytest.warns(FutureWarning):
                md_bytes = pyslim.encode_population(md)
            with pytest.warns(FutureWarning):
                new_md = pyslim.decode_population(md_bytes)
            assert md == new_md


class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests the table annotation methods.
    '''

    def test_annotate_mutations(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        tables = ts.dump_tables()
        new_tables = ts.dump_tables()
        metadata = []
        for md in tskit.unpack_bytes(tables.mutations.metadata,
                                     tables.mutations.metadata_offset):
            with pytest.warns(FutureWarning):
                dm = pyslim.decode_mutation(md)
            with pytest.warns(FutureWarning):
                edm = pyslim.encode_mutation(dm)
            assert md == edm
            metadata.append(dm)

        with pytest.warns(FutureWarning):
            pyslim.annotate_mutation_metadata(new_tables, metadata)
        self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_nodes(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        tables = ts.dump_tables()
        new_tables = ts.dump_tables()
        metadata = []
        for md in tskit.unpack_bytes(tables.nodes.metadata,
                                       tables.nodes.metadata_offset):
            with pytest.warns(FutureWarning):
                dm = pyslim.decode_node(md)
            with pytest.warns(FutureWarning):
                edm = pyslim.encode_node(dm)
            assert md == edm
            metadata.append(dm)

        with pytest.warns(FutureWarning):
            pyslim.annotate_node_metadata(new_tables, metadata)
        self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_individuals(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        tables = ts.dump_tables()
        new_tables = ts.dump_tables()
        metadata = []
        for md in tskit.unpack_bytes(tables.individuals.metadata,
                                       tables.individuals.metadata_offset):
            with pytest.warns(FutureWarning):
                dm = pyslim.decode_individual(md)
            with pytest.warns(FutureWarning):
                edm = pyslim.encode_individual(dm)
            assert md == edm
            metadata.append(dm)

        with pytest.warns(FutureWarning):
            pyslim.annotate_individual_metadata(new_tables, metadata)
        self.assertTableCollectionsEqual(tables, new_tables)

    def test_annotate_populations(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        tables = ts.dump_tables()
        new_tables = ts.dump_tables()
        metadata = []
        for md in tskit.unpack_bytes(tables.populations.metadata,
                                     tables.populations.metadata_offset):
            with pytest.warns(FutureWarning):
                dm = pyslim.decode_population(md)
            with pytest.warns(FutureWarning):
                edm = pyslim.encode_population(dm)
            assert md == edm
            metadata.append(dm)

        with pytest.warns(FutureWarning):
            pyslim.annotate_population_metadata(new_tables, metadata)
        self.assertTableCollectionsEqual(tables, new_tables)


class TestDumpLoad(tests.PyslimTestCase):
    '''
    Test reading and writing.
    '''

    def verify_times(self, ts, slim_ts):
        gen = slim_ts.slim_generation
        assert ts.num_nodes == slim_ts.num_nodes
        # verify internal consistency
        for j in range(slim_ts.num_nodes):
            assert slim_ts.node(j).time == slim_ts.tables.nodes.time[j]
        # verify consistency between tree sequences
        for n1, n2 in zip(ts.nodes(), slim_ts.nodes()):
            assert n1.time == n2.time

    def test_load_tables(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        assert isinstance(ts, pyslim.SlimTreeSequence)
        tables = ts.dump_tables()
        new_ts = pyslim.load_tables(tables, legacy_metadata=True)
        assert isinstance(new_ts, pyslim.SlimTreeSequence)
        new_tables = new_ts.dump_tables()
        assert tables == new_tables

    def test_load(self, recipe):
        assert recipe["path"]["ts"] == recipe["path"]["ts_legacy_metadata"]
        fn = recipe["path"]["ts"]
        # load in msprime then switch
        msp_ts = tskit.load(fn)
        assert isinstance(msp_ts, tskit.TreeSequence)
        # transfer tables
        msp_tables = msp_ts.dump_tables()
        new_ts = pyslim.load_tables(msp_tables, legacy_metadata=True)
        assert isinstance(new_ts, pyslim.SlimTreeSequence)
        self.verify_times(msp_ts, new_ts)
        new_tables = new_ts.dump_tables()
        self.assertTableCollectionsEqual(msp_tables, new_tables)
        # convert directly
        new_ts = pyslim.SlimTreeSequence(msp_ts)
        assert isinstance(new_ts, pyslim.SlimTreeSequence)
        self.verify_times(msp_ts, new_ts)
        new_tables = new_ts.dump_tables()
        self.assertTableCollectionsEqual(msp_tables, new_tables)
        # load to pyslim from file
        slim_ts = pyslim.load(fn, legacy_metadata=True)
        assert isinstance(slim_ts, pyslim.SlimTreeSequence)
        slim_tables = slim_ts.dump_tables()
        self.assertTableCollectionsEqual(msp_tables, slim_tables)
        assert slim_ts.slim_generation == new_ts.slim_generation

    def test_dump_equality(self, recipe, tmp_path):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        tmp_file = os.path.join(tmp_path, "test_dump.trees")
        ts = recipe["ts_legacy_metadata"]
        ts.dump(tmp_file)
        ts2 = pyslim.load(tmp_file, legacy_metadata=True)
        assert ts.num_samples == ts2.num_samples
        assert ts.sequence_length == ts2.sequence_length
        assert ts.tables == ts2.dump_tables()
        assert ts.reference_sequence == ts2.reference_sequence


class TestIndividualMetadata(tests.PyslimTestCase):
    # Tests for extra stuff related to Individuals.

    def test_individual_derived_info(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        for j, ind in enumerate(ts.individuals()):
            a = ts.tables.individuals.metadata_offset[j]
            b = ts.tables.individuals.metadata_offset[j+1]
            raw_md = ts.tables.individuals.metadata[a:b]
            with pytest.warns(FutureWarning):
                md = pyslim.decode_individual(raw_md)
            assert ind.metadata == md
            assert ts.individual(j).metadata == md
            for n in ind.nodes:
                assert ts.node(n).population == ind.population
                assert ts.node(n).time == ind.time


class TestNodeMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Nodes.
    '''

    def test_node_derived_info(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        for j, node in enumerate(ts.nodes()):
            a = ts.tables.nodes.metadata_offset[j]
            b = ts.tables.nodes.metadata_offset[j+1]
            raw_md = ts.tables.nodes.metadata[a:b]
            with pytest.warns(FutureWarning):
                md = pyslim.decode_node(raw_md)
            assert node.metadata == md
            assert ts.node(j).metadata == md


class TestMutationMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Mutations.
    '''

    def test_mutation_derived_info(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        for j, mut in enumerate(ts.mutations()):
            a = ts.tables.mutations.metadata_offset[j]
            b = ts.tables.mutations.metadata_offset[j+1]
            raw_md = ts.tables.mutations.metadata[a:b]
            with pytest.warns(FutureWarning):
                md = pyslim.decode_mutation(raw_md)
            assert mut.metadata == md
            assert ts.mutation(j).metadata == md


class TestPopulationMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Populations.
    '''

    def test_population_derived_info(self, recipe):
        ts = recipe["ts_legacy_metadata"]
        for j, pop in enumerate(ts.populations()):
            a = ts.tables.populations.metadata_offset[j]
            b = ts.tables.populations.metadata_offset[j+1]
            raw_md = ts.tables.populations.metadata[a:b]
            with pytest.warns(FutureWarning):
                md = pyslim.decode_population(raw_md)
            assert pop.metadata == md
            assert ts.population(j).metadata == md
