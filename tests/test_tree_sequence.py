"""
Test cases for tree sequences.
"""
from __future__ import print_function
from __future__ import division

import tests
import unittest
import random
import os
import numpy as np
import pytest
import tskit
import msprime
import pyslim


class TestSlimTreeSequence(tests.PyslimTestCase):

    def clean_example(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.populations.add_row()
        tables.individuals.add_row()
        tables.nodes.add_row(time=0, flags=tskit.NODE_IS_SAMPLE, population=1, individual=0)
        tables.nodes.add_row(time=0, flags=tskit.NODE_IS_SAMPLE, population=1, individual=0)
        pyslim.annotate_defaults_tables(tables, model_type='nonWF', slim_generation=1)
        return tables

    def test_inconsistent_nodes(self):
        clean_tables = self.clean_example()
        tables = clean_tables.copy()
        tables.nodes.clear()
        for j, n in enumerate(clean_tables.nodes):
            tables.nodes.add_row(
                    time=n.time, flags=n.flags,
                    population=j,
                    individual=n.individual,
                    metadata=n.metadata)
        with pytest.raises(ValueError):
            pyslim.annotate_defaults_tables(tables, model_type='nonWF', slim_generation=1)
        ts = tables.tree_sequence()
        with pytest.raises(ValueError):
            _ = pyslim.SlimTreeSequence(ts)

    def test_inconsistent_times(self):
        clean_tables = self.clean_example()
        tables = clean_tables.copy()
        tables.nodes.clear()
        for j, n in enumerate(clean_tables.nodes):
            tables.nodes.add_row(time=j, flags=tskit.NODE_IS_SAMPLE, population=1, individual=0)
        ts = tables.tree_sequence()
        with pytest.raises(ValueError):
            _ = pyslim.SlimTreeSequence(ts)

    def test_bad_metadata(self):
        clean_tables = self.clean_example()
        tables = clean_tables.copy()
        tables.metadata_schema = tskit.MetadataSchema({"type": "object", "codec": "json"})
        tables.metadata = {}
        ts = tables.tree_sequence()
        with pytest.raises(ValueError):
            _ = pyslim.SlimTreeSequence(ts)

    # tmp_path is a pytest fixture, and is a pathlib.Path object
    def test_slim_generation(self, tmp_path):
        # tests around awkward backwards-compatible patch for setting slim_generation
        ts = self.get_slim_example(name="recipe_nonWF")
        assert ts.slim_generation == ts.metadata['SLiM']['generation']
        new_sg = 12345
        ts.slim_generation = new_sg
        assert ts.slim_generation == new_sg
        # check persists through dump/load
        temp_file = tmp_path / "temp.trees"
        ts.dump(temp_file.name)
        loaded_ts = pyslim.load(temp_file.name)
        assert loaded_ts.slim_generation == new_sg
        assert loaded_ts.metadata['SLiM']['generation'] == new_sg
        # check persists through recapitate
        recap = ts.recapitate(recombination_rate=1e-8)
        assert recap.slim_generation == new_sg
        # check persists through simplify
        simp = ts.simplify(ts.samples())
        assert simp.slim_generation == new_sg


class TestSlimTime(tests.PyslimTestCase):
    # Tests for slim_time()

    def test_slim_time(self):
        for ts in self.get_slim_examples(init_mutated=False):
            for mut in ts.mutations():
                mut_time = max([x['slim_time'] for x in mut.metadata['mutation_list']])
                assert mut_time == ts.slim_time(mut.time)
        # the mutations in "init_mutated" examples have mutations that are *added*
        # in *early*, and so their times match in that stage.
        for ts in self.get_slim_examples(init_mutated=True):
            for mut in ts.mutations():
                mut_time = max([x['slim_time'] for x in mut.metadata['mutation_list']])
                assert mut_time == ts.slim_time(mut.time, stage="early")


class TestMutate(tests.PyslimTestCase):
    # Tests for making a tree sequence a SlimTreeSequence
    # again after msprime.mutate.

    def test_mutate(self):
        for ts in self.get_slim_examples():
            mts = msprime.mutate(ts, rate=1e-8, random_seed=5)
            pts = pyslim.SlimTreeSequence(mts)
            assert ts.metadata == pts.metadata


class TestRecapitate(tests.PyslimTestCase):
    '''
    Tests for recapitation.
    '''

    def check_recap_consistency(self, ts, recap):
        assert ts.slim_generation == recap.slim_generation
        assert all(tree.num_roots == 1 for tree in recap.trees())

        ts_samples = list(ts.samples())
        for u in recap.samples():
            n1 = recap.node(u)
            assert n1.individual >= 0
            i1 = recap.individual(n1.individual)
            remembered = ((pyslim.INDIVIDUAL_REMEMBERED & i1.flags) > 0)
            alive = ((pyslim.INDIVIDUAL_ALIVE & i1.flags) > 0)
            assert alive or remembered
            assert u in ts_samples
            n2 = ts.node(u)
            assert n1.time == n2.time
            assert n1.individual == n2.individual
            assert n1.flags == n2.flags
            assert n1.metadata == n2.metadata
            assert n1.population == n2.population
        assert ts.num_populations <= recap.num_populations
        for k in range(ts.num_populations):
            p1 = ts.population(k)
            p2 = recap.population(k)
            assert p1.metadata == p2.metadata

    def test_recapitate_errors(self):
        ts = next(self.get_slim_examples())
        with pytest.raises(ValueError):
            _ = ts.recapitate(
                        recombination_rate=0.0,
                        keep_first_generation=True)

    def test_recapitation(self):
        for ts in self.get_slim_examples():
            if ts.num_populations <= 2:
                # if not we need migration rates
                recomb_rate = 1.0 / ts.sequence_length
                recap = ts.recapitate(recombination_rate=recomb_rate)
                # there should be no new mutations
                assert ts.num_mutations == recap.num_mutations
                assert ts.num_sites == recap.num_sites
                assert list(ts.tables.sites.position) == list(recap.tables.sites.position)
                self.check_recap_consistency(ts, recap)
                for t in recap.trees():
                    assert t.num_roots == 1

                recap = ts.recapitate(recombination_rate=recomb_rate, Ne=1e-6)
                self.check_recap_consistency(ts, recap)
                if ts.slim_generation < 200:
                    for t in recap.trees():
                        assert t.num_roots == 1
                        assert abs(recap.node(t.root).time - recap.slim_generation) < 1e-4


class TestIndividualMetadata(tests.PyslimTestCase):
    # Tests for extra stuff related to Individuals.

    def test_individual_derived_info(self):
        for ts in self.get_slim_examples():
            for ind in ts.individuals():
                for n in ind.nodes:
                    assert ts.node(n).population == ind.population
                    assert ts.node(n).time == ind.time

    def test_individual_embellishments(self):
        # Test the individual additional information.
        for ts in self.get_slim_examples():
            is_wf = (ts.metadata["SLiM"]["model_type"] == "WF")
            for j, ind in enumerate(ts.individuals()):
                assert ts.individual_times[j] == ind.time
                if is_wf:
                    assert ts.individual_ages[j] == 0
                else:
                    assert ts.individual_ages[j] == ind.metadata["age"]
                assert ts.individual_populations[j] == ind.population
                assert np.array_equal(ts.individual_locations[j], ind.location)

    def test_first_gen_nodes(self):
        # check that all the roots of the trees are present
        for ts in self.get_slim_examples():
            root_time = ts.slim_generation
            if (ts.metadata['SLiM']['stage'] == 'early'
                    or ts.metadata['SLiM']['model_type'] == 'nonWF'):
                root_time -= 1
            for t in ts.trees():
                for u in t.roots:
                    assert ts.node(u).time == root_time


class TestMutationMetadata(tests.PyslimTestCase):
    '''
    Tests for extra stuff related to Mutations.
    '''

    def test_slim_time(self):
        # check that slim_times make sense
        for ts in self.get_slim_examples(init_mutated=False):
            # Mutation's slim_times are one less than the corresponding node's slim times
            # in WF models, but not in WF models, for some reason.
            is_wf = (ts.metadata["SLiM"]["model_type"] == "WF")
            for mut in ts.mutations():
                node_slim_time = ts.slim_generation - ts.node(mut.node).time
                mut_slim_time = max([u["slim_time"] for u in mut.metadata["mutation_list"]])
                assert node_slim_time >= mut_slim_time


class TestIndividualAges(tests.PyslimTestCase):
    # tests for individuals_alive_at and individual_ages_at

    def test_errors(self):
        ts = next(self.get_slim_examples(everyone=True))
        for stage in ['abcd', 10, []]:
            with pytest.raises(ValueError):
                ts.individuals_alive_at(0, stage=stage)
            with pytest.raises(ValueError):
                ts.individuals_alive_at(0, remembered_stage=stage)
            with pytest.raises(ValueError):
                ts.individual_ages_at(0, stage=stage)

    def test_mismatched_remembered_stage(self):
        for ts, ex in self.get_slim_examples(pedigree=True, WF=True, return_info=True):
            info = ex['info']
            if "remembered_early" in ex:
                with pytest.warns(UserWarning):
                    ts.individuals_alive_at(0, remembered_stage="late")
            else:
                with pytest.warns(UserWarning):
                    ts.individuals_alive_at(0, remembered_stage="early")

    def test_population(self):
        for ts in self.get_slim_examples(multipop=True, remembered_early=False):
            all_inds = ts.individuals_alive_at(0)
            for p in range(ts.num_populations):
                sub_inds = ts.individuals_alive_at(0, population=p)
                assert set(sub_inds) == set(all_inds[ts.individual_populations == p])
                sub_inds = ts.individuals_alive_at(0, population=[p])
                assert set(sub_inds) == set(all_inds[ts.individual_populations == p])
            sub_inds = ts.individuals_alive_at(0, population=np.arange(p))
            assert set(sub_inds) == set(all_inds[ts.individual_populations != p])

    def test_samples_only(self):
        for ts in self.get_slim_examples(nonWF=True, remembered_early=False):
            all_inds = ts.individuals_alive_at(0)
            assert set(all_inds) == set(ts.individuals_alive_at(0, samples_only=False))
            sub_inds = np.random.choice(all_inds, size=min(len(all_inds), 4), replace=False)
            flags = np.array([n.flags & (tskit.NODE_IS_SAMPLE * n.individual in sub_inds)
                              for n in ts.nodes()], dtype=np.uint32)
            tables = ts.tables
            tables.nodes.flags = flags
            new_ts = pyslim.SlimTreeSequence(tables.tree_sequence())
            assert set(sub_inds) == set(new_ts.individuals_alive_at(0, samples_only=True))

    def test_after_simplify(self):
        for ts in self.get_slim_examples(remembered_early=False):
            sts = ts.simplify()
            orig_inds = ts.individuals_alive_at(0)
            simp_inds = sts.individuals_alive_at(0)
            odict = {ts.individual(i).metadata["pedigree_id"]: i for i in orig_inds}
            sdict = {sts.individual(i).metadata["pedigree_id"]: i for i in simp_inds}
            for slim_id in odict:
                i = odict[slim_id]
                ind = ts.individual(i)
                n = ts.node(ind.nodes[0])
                if n.flags & tskit.NODE_IS_SAMPLE:
                    assert slim_id in sdict

    def test_ages(self):
        for ts, ex in self.get_slim_examples(pedigree=True, return_info=True):
            info = ex['info']
            remembered_stage = 'early' if 'remembered_early' in ex else 'late'
            assert remembered_stage == ts.metadata['SLiM']['stage']
            max_time_ago = ts.slim_generation
            if remembered_stage == 'early':
                max_time_ago -= 1
            for time in range(0, max_time_ago):
                # if written out during 'early' in a WF model,
                # tskit time 0 will be the SLiM time step *before* slim_generation
                slim_time = ts.slim_generation - time
                if remembered_stage == 'early' and ts.metadata["SLiM"]["model_type"] == "WF":
                    slim_time -= 1
                if remembered_stage == 'early' and time == 0:
                    # if we remember in early we don't know who's still there
                    # in late of the last time step
                    check_stages = ('early',)
                else:
                    check_stages = ('early', 'late')
                for stage in check_stages:
                    alive = ts.individuals_alive_at(
                                time,
                                stage=stage,
                                remembered_stage=remembered_stage)
                    ages = ts.individual_ages_at(
                                time,
                                stage=stage,
                                remembered_stage=remembered_stage)
                    for ind in ts.individuals():
                        if 'everyone' in ex or ind.time == 0:
                            slim_id = ind.metadata["pedigree_id"]
                            assert slim_id in info
                            slim_alive = (slim_time, stage) in info[slim_id]['age']
                            pyslim_alive = ind.id in alive
                            print(time, (slim_time, stage))
                            print(ind)
                            print(info[slim_id])
                            print(slim_alive, pyslim_alive)
                            assert slim_alive == pyslim_alive
                            if slim_alive:
                                slim_age = info[slim_id]['age'][(slim_time, stage)]
                                if ts.metadata["SLiM"]["model_type"] == "WF":
                                    # SLiM records -1 but we return 0 in late and 1 in early
                                    slim_age = 0 + (stage == 'early')
                                print('age:', ages[ind.id], slim_age)
                                assert ages[ind.id] == slim_age
                            else:
                                assert np.isnan(ages[ind.id])


class TestHasIndividualParents(tests.PyslimTestCase):

    def verify_has_parents(self, ts):
        right_answer = np.repeat(True, ts.num_individuals)
        node_indivs = ts.tables.nodes.individual
        parent_ids = [set() for _ in ts.individuals()]
        node_parent_ids = [set() for _ in ts.nodes()]
        for t in ts.trees():
            for i in ts.individuals():
                if len(i.nodes) != 2:
                    right_answer[i.id] = False
                for n in i.nodes:
                    pn = t.parent(n)
                    if pn == tskit.NULL:
                        right_answer[i.id] = False
                    else:
                        p = node_indivs[t.parent(n)]
                        if p == tskit.NULL:
                            right_answer[i.id] = False
                        else:
                            ptime = ts.individual_times[p]
                            if ts.metadata["SLiM"]["model_type"] == "WF":
                                if i.time + 1 != ptime:
                                    right_answer[i.id] = False
                            else:
                                pdeath = ptime - ts.individual_ages[p]
                                if i.time + 1 < pdeath:
                                    right_answer[i.id] = False
                            parent_ids[i.id].add(p)
                            node_parent_ids[n].add(p)
        for j, p in enumerate(parent_ids):
            if len(p) == 0:
                right_answer[j] = False
        for j, p in enumerate(node_parent_ids):
            if len(p) != 1:
                ind = ts.node(j).individual
                if ind != tskit.NULL:
                    right_answer[ts.node(j).individual] = False
        has_parents = ts.has_individual_parents()
        for j, (a, b) in enumerate(zip(right_answer, has_parents)):
            if a != b:
                print("------------")
                print(j, a, b)
                print(parent_ids[j])
                for n in ts.individual(j).nodes:
                    print(n, node_parent_ids[n])
                print("------------")
        assert np.array_equal(right_answer, has_parents)

    def get_first_gen(self, ts):
        root_time = ts.metadata["SLiM"]["generation"]
        if ts.metadata['SLiM']['model_type'] != 'WF' or ts.metadata['SLiM']['stage'] != 'late':
            root_time -= 1
        first_gen = set(ts.tables.nodes.individual[ts.tables.nodes.time == root_time])
        first_gen.discard(tskit.NULL)
        return np.array(list(first_gen), dtype='int')

    def test_everyone(self):
        # since everyone is recorded, only the initial individuals should
        # not have parents
        for ts in self.get_slim_examples(everyone=True):
            right_answer = np.repeat(True, ts.num_individuals)
            first_gen = self.get_first_gen(ts)
            right_answer[first_gen] = False
            has_parents = ts.has_individual_parents()
            assert np.array_equal(right_answer, has_parents)
            self.verify_has_parents(ts)

    def test_post_recap(self):
        # the same should be true after recapitation
        for ts in self.get_slim_examples(everyone=True):
            right_answer = np.repeat(True, ts.num_individuals)
            first_gen = self.get_first_gen(ts)
            right_answer[first_gen] = False
            assert(ts.num_populations <= 2)
            ts = ts.recapitate(recombination_rate=0.01)
            assert(ts.num_individuals == ts.num_individuals)
            has_parents = ts.has_individual_parents()
            assert np.array_equal(right_answer, has_parents)
            self.verify_has_parents(ts)

    def test_post_simplify(self):
        for ts in self.get_slim_examples(everyone=True):
            keep_indivs = np.random.choice(
                    np.where(ts.individual_times < ts.slim_generation - 1)[0],
                    size=30, replace=False)
            keep_nodes = []
            for i in keep_indivs:
                keep_nodes.extend(ts.individual(i).nodes)
            ts = ts.simplify(samples=keep_nodes, filter_individuals=True)
            assert(ts.num_populations <= 2)
            ts = ts.recapitate(recombination_rate=0.01)
            has_parents = ts.has_individual_parents()
            assert sum(has_parents) > 0
            self.verify_has_parents(ts)

    def test_pedigree(self):
        for ts, ex in self.get_slim_examples(pedigree=True, return_info=True):
            has_parents = ts.has_individual_parents()
            info = ex['info']
            slim_map = {}
            for ind in ts.individuals():
                slim_map[ind.metadata["pedigree_id"]] = ind.id
            for hasp, ind in zip(has_parents, ts.individuals()):
                slim_parents = info[ind.metadata["pedigree_id"]]['parents']
                slim_hasp = len(slim_parents) > 0
                for p in slim_parents:
                    if p not in slim_map:
                        slim_hasp = False
                assert hasp == slim_hasp


class TestSimplify(tests.PyslimTestCase):
    '''
    Our simplify() is just a wrapper around the tskit simplify.
    '''

    def test_simplify(self):
        for ts in self.get_slim_examples():
            sts = ts.simplify(map_nodes=False)
            assert ts.sequence_length == sts.sequence_length
            assert type(ts) == type(sts)
            assert sts.samples()[0] == 0
            sts, _ = ts.simplify(map_nodes=True)
            assert ts.sequence_length == sts.sequence_length
            assert type(ts) == type(sts)
            assert sts.samples()[0] == 0


class TestReferenceSequence(tests.PyslimTestCase):
    '''
    Test for operations involving the reference sequence
    '''

    def test_reference_sequence(self):
        for ts in self.get_slim_examples():
            if ts.num_mutations > 0:
                mut_md = ts.mutation(0).metadata
                has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
                if not has_nucleotides:
                    assert ts.reference_sequence == None
                else:
                    assert type(ts.reference_sequence) == type('')
                    assert len(ts.reference_sequence) == ts.sequence_length
                    for u in ts.reference_sequence:
                        assert u in pyslim.NUCLEOTIDES
                sts = ts.simplify(ts.samples()[:2])
                assert sts.reference_sequence == ts.reference_sequence

    def test_mutation_at_errors(self):
        for ts in self.get_slim_examples():
            u = ts.samples()[0]
            with pytest.raises(ValueError):
                ts.mutation_at(-2, 3)
            with pytest.raises(ValueError):
                ts.mutation_at(u, -3)
            with pytest.raises(ValueError):
                ts.mutation_at(ts.num_nodes + 2, 3)
            with pytest.raises(ValueError):
                ts.mutation_at(u, ts.sequence_length)

    def test_nucleotide_at_errors(self):
        for ts in self.get_slim_examples():
            u = ts.samples()[0]
            if ts.num_mutations > 0:
                mut_md = ts.mutation(0).metadata
                has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
                if not has_nucleotides:
                    with pytest.raises(ValueError):
                        ts.nucleotide_at(u, 3)

    def test_mutation_at(self):
        random.seed(42)
        for ts in self.get_slim_examples():
            for _ in range(100):
                node = random.randint(0, ts.num_nodes - 1)
                pos = random.randint(0, ts.sequence_length - 1)
                tree = ts.at(pos)
                parent = tree.parent(node)
                a = ts.mutation_at(node, pos)
                if parent == tskit.NULL:
                    assert a == tskit.NULL
                else:
                    b = ts.mutation_at(parent, pos)
                    c = ts.mutation_at(node, pos, ts.node(parent).time)
                    assert b == c
                    for k in np.where(node == ts.tables.mutations.node)[0]:
                        mut = ts.mutation(k)
                        if ts.site(mut.site).position == pos:
                            b = mut.id
                    assert a == b

    def test_nucleotide_at(self):
        random.seed(42)
        for ts in self.get_slim_examples():
            if ts.num_mutations > 0:
                mut_md = ts.mutation(0).metadata
                has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
                if has_nucleotides:
                    for _ in range(100):
                        node = random.randint(0, ts.num_nodes - 1)
                        pos = random.randint(0, ts.sequence_length - 1)
                        tree = ts.at(pos)
                        parent = tree.parent(node)
                        a = ts.nucleotide_at(node, pos)
                        if parent == tskit.NULL:
                            nuc = ts.reference_sequence[int(pos)]
                            assert a == pyslim.NUCLEOTIDES.index(nuc)
                        else:
                            b = ts.nucleotide_at(parent, pos)
                            c = ts.nucleotide_at(node, pos, ts.node(parent).time)
                            assert b == c
                            for k in np.where(node == ts.tables.mutations.node)[0]:
                                mut = ts.mutation(k)
                                if ts.site(mut.site).position == pos:
                                    b = mut.metadata["mutation_list"][0]["nucleotide"]
                            assert a == b

class TestDeprecations(tests.PyslimTestCase):

    def test_first_gen(self):
        ts = next(self.get_slim_examples())
        with pytest.warns(FutureWarning):
            _ = ts.first_generation_individuals()
