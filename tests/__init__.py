"""
Common code for the pyslim test cases.
"""
import os
import json
import random
import base64
import warnings

import pyslim
import tskit
import msprime
import pytest
import attr
import numpy as np


class PyslimTestCase:
    '''
    Base class for test cases in pyslim.
    '''

    def verify_haplotype_equality(self, ts, slim_ts):
        assert ts.num_sites == slim_ts.num_sites
        for j, v1, v2 in zip(range(ts.num_sites), ts.variants(),
                             slim_ts.variants()):
            g1 = [v1.alleles[x] for x in v1.genotypes]
            g2 = [v2.alleles[x] for x in v2.genotypes]
            assert np.array_equal(g1, g2)


    def assertTablesEqual(self, t1, t2, label=''):
        # make it easy to see what's wrong
        if hasattr(t1, "metadata_schema"):
            if t1.metadata_schema != t2.metadata_schema:
                print(f"{label} :::::::::: t1 ::::::::::::")
                print(t1.metadata_schema)
                print(f"{label} :::::::::: t2 ::::::::::::")
                print(t2.metadata_schema)
            assert t1.metadata_schema == t2.metadata_schema
        if t1.num_rows != t2.num_rows:
            print(f"{label}: t1.num_rows {t1.num_rows} != {t2.num_rows} t2.num_rows")
        for k, (e1, e2) in enumerate(zip(t1, t2)):
            if e1 != e2:
                print(f"{label} :::::::::: t1 ({k}) ::::::::::::")
                print(e1)
                print(f"{label} :::::::::: t2 ({k}) ::::::::::::")
                print(e2)
            assert e1 == e2
        assert t1.num_rows == t2.num_rows
        assert t1 == t2

    def assertMetadataEqual(self, t1, t2):
        # check top-level metadata, first the parsed version:
        assert t1.metadata_schema == t2.metadata_schema
        assert t1.metadata == t2.metadata
        # and now check the underlying bytes
        # TODO: use the public interface if https://github.com/tskit-dev/tskit/issues/832 happens
        md1 = t1._ll_tables.metadata
        md2 = t2._ll_tables.metadata
        assert md1 == md2

    def verify_trees_equal(self, ts1, ts2):
        # check that trees are equal by checking MRCAs between randomly
        # chosen nodes with matching slim_ids
        random.seed(23)
        assert ts1.sequence_length == ts2.sequence_length
        if isinstance(ts1, tskit.TableCollection):
            ts1 = ts1.tree_sequence()
        if isinstance(ts2, tskit.TableCollection):
            ts2 = ts2.tree_sequence()
        map1 = {}
        for j, n in enumerate(ts1.nodes()):
            if n.metadata is not None:
                map1[n.metadata['slim_id']] = j
        map2 = {}
        for j, n in enumerate(ts2.nodes()):
            if n.metadata is not None:
                map2[n.metadata['slim_id']] = j
        assert set(map1.keys()) == set(map2.keys())
        print(ts1)
        print(map1)
        print(ts2)
        print(map2)
        sids = list(map1.keys())
        for sid in sids:
            n1 = ts1.node(map1[sid])
            n2 = ts2.node(map2[sid])
            assert n1.time == n2.time
            assert n1.metadata == n2.metadata
            i1 = ts1.individual(n1.individual)
            i2 = ts2.individual(n2.individual)
            if i1.metadata != i2.metadata:
                print("i1: ", i1.metadata)
                print("i2: ", i2.metadata)
            assert i1.metadata == i2.metadata
        for _ in range(10):
            pos = random.uniform(0, ts1.sequence_length)
            t1 = ts1.at(pos)
            t2 = ts2.at(pos)
            for _ in range(10):
                a, b = random.choices(sids, k=2)
                print(a, b, map1[a], map1[b], map2[a], map2[b])
                print(t1)
                print('a', t1.time(map1[a]))
                print('b', t1.time(map1[b]))
                print('1', t1.tmrca(map1[a], map1[b]))
                print('2', t2.tmrca(map2[a], map2[b]))
                assert t1.tmrca(map1[a], map1[b]) == t2.tmrca(map2[a], map2[b])

    def assertTableCollectionsEqual(self, t1, t2,
            skip_provenance=False, check_metadata_schema=True,
            reordered_individuals=False):
        if isinstance(t1, tskit.TreeSequence):
            t1 = t1.dump_tables()
        if isinstance(t2, tskit.TreeSequence):
            t2 = t2.dump_tables()
        t1_samples = [(n.metadata['slim_id'], j) for j, n in enumerate(t1.nodes) if (n.flags & tskit.NODE_IS_SAMPLE)]
        t1_samples.sort()
        t2_samples = [(n.metadata['slim_id'], j) for j, n in enumerate(t2.nodes) if (n.flags & tskit.NODE_IS_SAMPLE)]
        t2_samples.sort()
        t1.simplify([j for (_, j) in t1_samples], record_provenance=False)
        t2.simplify([j for (_, j) in t2_samples], record_provenance=False)
        if skip_provenance is True:
            t1.provenances.clear()
            t2.provenances.clear()
        if skip_provenance == -1:
            assert t1.provenances.num_rows + 1 == t2.provenances.num_rows
            t2.provenances.truncate(t1.provenances.num_rows)
            assert t1.provenances.num_rows == t2.provenances.num_rows
        if check_metadata_schema:
            # this is redundant now, but will help diagnose if things go wrong
            assert t1.metadata_schema.schema == t2.metadata_schema.schema
            assert t1.populations.metadata_schema.schema == t2.populations.metadata_schema.schema
            assert t1.individuals.metadata_schema.schema == t2.individuals.metadata_schema.schema
            assert t1.nodes.metadata_schema.schema == t2.nodes.metadata_schema.schema
            assert t1.edges.metadata_schema.schema == t2.edges.metadata_schema.schema
            assert t1.sites.metadata_schema.schema == t2.sites.metadata_schema.schema
            assert t1.mutations.metadata_schema.schema == t2.mutations.metadata_schema.schema
            assert t1.migrations.metadata_schema.schema == t2.migrations.metadata_schema.schema
        if not check_metadata_schema:
            # need to pull out metadata to compare as dicts before zeroing the schema
            m1 = t1.metadata
            m2 = t2.metadata
            ms = tskit.MetadataSchema(None)
            for t in (t1, t2):
                t.metadata_schema = ms
                t.populations.metadata_schema = ms
                t.individuals.metadata_schema = ms
                t.nodes.metadata_schema = ms
                t.edges.metadata_schema = ms
                t.sites.metadata_schema = ms
                t.mutations.metadata_schema = ms
                t.migrations.metadata_schema = ms
            t1.metadata = b''
            t2.metadata = b''
            assert m1 == m2
        if reordered_individuals:
            ind1 = {i.metadata['pedigree_id']: j for j, i in enumerate(t1.individuals)}
            ind2 = {i.metadata['pedigree_id']: j for j, i in enumerate(t2.individuals)}
            for pid in ind1:
                if not pid in ind2:
                    print("not in t2:", ind1[pid])
                assert pid in ind2
                if t1.individuals[ind1[pid]] != t2.individuals[ind2[pid]]:
                    print("t1:", t1.individuals[ind1[pid]])
                    print("t2:", t2.individuals[ind2[pid]])
                assert t1.individuals[ind1[pid]] == t2.individuals[ind2[pid]]
            for pid in ind2:
                if not pid in ind1:
                    print("not in t1:", ind2[pid])
                assert pid in ind1
            t1.individuals.clear()
            t2.individuals.clear()
        # go through one-by-one so we know which fails
        self.assertTablesEqual(t1.populations, t2.populations, "populations")
        self.assertTablesEqual(t1.individuals, t2.individuals, "individuals")
        self.assertTablesEqual(t1.nodes, t2.nodes, "nodes")
        self.assertTablesEqual(t1.edges, t2.edges, "edges")
        self.assertTablesEqual(t1.sites, t2.sites, "sites")
        self.assertTablesEqual(t1.mutations, t2.mutations, "mutations")
        self.assertTablesEqual(t1.migrations, t2.migrations, "migrations")
        self.assertTablesEqual(t1.provenances, t2.provenances, "provenances")
        self.assertMetadataEqual(t1, t2)
        assert t1.sequence_length == t2.sequence_length
        if t1.reference_sequence.data != t2.reference_sequence.data:
            print(t1.reference_sequence.data, " != ", t2.reference_sequence.data)
        assert t1.reference_sequence.data == t2.reference_sequence.data
        # time units are not passed on: https://github.com/tskit-dev/msprime/issues/1951
        # if t1.time_units != t2.time_units:
        #     print(t1.time_units, " != ", t2.time_units)
        # assert t1.time_units == t2.time_units
        # assert t1 == t2

    def do_recapitate(self, ts, *args, **kwargs):
        if ts.metadata['SLiM']['model_type'] != "WF":
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
                recap = pyslim.recapitate(ts, *args, **kwargs)
        else:
            recap = pyslim.recapitate(ts, *args, **kwargs)
        return recap

