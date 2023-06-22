"""
Test cases for tree sequences.
"""
import random
import numpy as np
import os
import json

import pytest
import tskit
import msprime
import pyslim

import tests

from .recipe_specs import recipe_eq

class TestSlimTime(tests.PyslimTestCase):
    # Tests for slim_time()

    @pytest.mark.parametrize('recipe', recipe_eq(exclude="long"), indirect=True)
    def test_slim_time(self, recipe):
        ts = recipe["ts"]
        if "init_mutated" not in recipe:
            for mut in ts.mutations():
                mut_time = max([x['slim_time'] for x in mut.metadata['mutation_list']])
                assert mut_time == pyslim.slim_time(ts, mut.time)
        # the mutations in "init_mutated" examples have mutations that are *added*
        # in *early*, and so their times match in that stage.
        else:
            for mut in ts.mutations():
                mut_time = max([x['slim_time'] for x in mut.metadata['mutation_list']])
                assert mut_time == pyslim.slim_time(ts, mut.time, stage="early")

class TestNextMutationID(tests.PyslimTestCase):
    '''
    Tests for the function that returns the largest SLiM mutation ID.
    '''
    def test_next_id(self, recipe):
        ts = recipe["ts"]
        mt_ids_str = ','.join(tskit.unpack_strings(ts.tables.mutations.derived_state, ts.tables.mutations.derived_state_offset))
        mt_ids = [int(i or 0) for i in mt_ids_str.split(',')]
        max_mt_id = max(mt_ids)
        assert max_mt_id + 1 == pyslim.next_slim_mutation_id(ts)

    @pytest.mark.parametrize('recipe', recipe_eq("adds_mutations"), indirect=True)
    def test_reload_slim(self, recipe, helper_functions, tmp_path):
        ts = recipe["ts"]
        rts = self.do_recapitate(
                    ts,
                    recombination_rate=1e-8,
                    ancestral_Ne=100,
                    random_seed=875,
        )
        next_id = pyslim.next_slim_mutation_id(rts)
        mts = msprime.sim_mutations(
                rts,
                rate=1e-4,
                keep=True,
                model=msprime.SLiMMutationModel(type=1, next_id=next_id),
        )
        assert mts.num_mutations > rts.num_mutations
        rrts = helper_functions.run_slim_restart(
                mts,
                "restart_nucleotides_WF.slim",
                tmp_path,
                WF=False,
        )
        # nothing should change
        assert pyslim.next_slim_mutation_id(mts) == pyslim.next_slim_mutation_id(rrts)
        assert rrts.num_mutations == mts.num_mutations

        def test_invalid_derived_state(self):
            ts = msprime.sim_ancestry(
                    4,
                    sequence_length=10,
                    population_size=10,
                    random_seed=10,
            )
            mts = msprime.sim_mutations(ts,
                    model="jc69",
                    rate=0.5,
                    random_seed=23)
            with pytest.raises(ValueError, match="need to be coercible to int"):
                pyslim.next_slim_mutation_id(mts)

class TestRecapitate(tests.PyslimTestCase):
    '''
    Tests for recapitation.
    '''

    def check_recap_consistency(self, ts, recap, with_ancestral_Ne=True):
        assert ts.metadata['SLiM']['tick'] == recap.metadata['SLiM']['tick']
        assert ts.metadata['SLiM']['cycle'] == recap.metadata['SLiM']['cycle']
        assert ts.metadata['SLiM']['stage'] == recap.metadata['SLiM']['stage']
        assert ts.metadata['SLiM']['name'] == recap.metadata['SLiM']['name']
        assert all(tree.num_roots == 1 for tree in recap.trees())
        assert ts.has_reference_sequence() == recap.has_reference_sequence()
        if ts.has_reference_sequence():
            assert ts.reference_sequence.data == recap.reference_sequence.data

        root_times = list(set([ts.node(n).time for t in ts.trees() for n in t.roots]))
        assert len(root_times) == 1
        if with_ancestral_Ne:
            # check that time recorded in provenance is correct
            assert recap.num_provenances == 2
            recap_prov = json.loads(recap.provenance(1).record)
            recap_events = recap_prov['parameters']['demography']['events']
            assert len(recap_events) == 1
            recap_time = recap_events[0]['time']
            assert np.allclose(recap_time, root_times[0])
            # the oldest nodes in all trees with SLiM provenance should be at that time
            if recap.num_nodes <= 500: # takes a long time otherwise
                for t in recap.trees():
                    for n in t.nodes():
                        rn = recap.node(n)
                        if rn.metadata is not None:
                            p = t.parent(n)
                            if p == tskit.NULL or recap.node(p).metadata is None:
                                assert rn.time == root_times[0]

        ts_samples = list(ts.samples())
        for u in recap.samples():
            n1 = recap.node(u)
            assert n1.individual >= 0
            i1 = recap.individual(n1.individual)
            remembered = ((pyslim.INDIVIDUAL_REMEMBERED & i1.flags) > 0)
            retained = ((pyslim.INDIVIDUAL_RETAINED & i1.flags) > 0)
            alive = ((pyslim.INDIVIDUAL_ALIVE & i1.flags) > 0)
            assert alive or remembered or retained
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
        # find ancestral pop in which recapitation happens
        tables = ts.tables
        anc_nodes = np.where(tables.nodes.time > ts.metadata['SLiM']['tick'])[0]
        if len(anc_nodes) > 0:
            for pop in ts.populations():
                if pop.metadata is not None and pop.metadata['name'] == "ancestral":
                    break
            assert pop.metadata['name'] == "ancestral"
            assert np.all(tables.nodes.population[anc_nodes] == pop.id)

    # Just test on the first recipe
    @pytest.mark.parametrize('recipe', [next(recipe_eq())], indirect=True)
    def test_recapitate_errors(self, recipe):
        ts = recipe["ts"]
        with pytest.raises(ValueError, match="cannot specify both `demography` and `ancestral_Ne`"):
            _ = self.do_recapitate(
                        ts,
                        recombination_rate=0.0,
                        demography=msprime.Demography.from_tree_sequence(ts),
                        ancestral_Ne=10,
                        random_seed=123,
            )

    def test_root_mismatch_error(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, random_seed=12)
        recap_time = 100
        ts = pyslim.annotate(ts, model_type="nonWF", tick=recap_time)
        assert ts.node(ts.first().roots[0]).time < recap_time
        with pytest.raises(ValueError, match="at the time expected by recapitate"):
            rts = self.do_recapitate(ts, ancestral_Ne=10)

    def test_unique_names(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, random_seed=12, end_time=1.0)
        # adjust seed if not
        assert min([t.num_roots for t in ts.trees()]) > 1
        ts = pyslim.annotate(ts, model_type="nonWF", tick=1)
        t = ts.dump_tables()
        md = t.populations[0].metadata
        md.update({"name": "ancestral"})
        t.populations[0] = t.populations[0].replace(metadata=md)
        md.update({"name": "ancestral_ancestral", "slim_id": 1})
        t.populations.add_row(metadata=md)
        ts = t.tree_sequence()
        rts = self.do_recapitate(ts, ancestral_Ne=10)
        names = [pop.metadata['name'] for pop in rts.populations()]
        assert len(set(names)) == len(names)
        assert names[0] == "ancestral"
        assert names[-2] == "ancestral_ancestral"

    def test_recapitation(self, recipe):
        ts = recipe["ts"]
        recomb_rate = 1.0 / ts.sequence_length
        recap = self.do_recapitate(
                    ts,
                    recombination_rate=recomb_rate,
                    ancestral_Ne=10,
                    random_seed=5
        )
        # there should be no new mutations
        assert ts.num_mutations == recap.num_mutations
        assert ts.num_sites == recap.num_sites
        assert list(ts.tables.sites.position) == list(recap.tables.sites.position)
        self.check_recap_consistency(ts, recap)

        if ts.metadata['SLiM']['tick'] < 200:
            old_root_time = np.max(ts.tables.nodes.time)
            for t in recap.trees():
                assert t.num_roots == 1
                assert recap.node(t.root).time >= old_root_time

    @pytest.mark.parametrize('recipe', recipe_eq(exclude="long"), indirect=True)
    def test_with_recomb_map(self, recipe):
        ts = recipe["ts"]
        recomb_rate = 1.0 / ts.sequence_length
        recombination_map = msprime.RateMap(
                   position = [0.0, ts.sequence_length],
                   rate = [recomb_rate]
        )
        recap = self.do_recapitate(
                    ts,
                    recombination_rate=recombination_map,
                    ancestral_Ne=1e-6,
                    random_seed=875,
        )
        self.check_recap_consistency(ts, recap)

    @pytest.mark.parametrize('recipe', recipe_eq("multipop"), indirect=True)
    def test_with_demography(self, recipe):
        ts = recipe["ts"]
        recomb_rate = 1.0 / ts.sequence_length
        demography = msprime.Demography.from_tree_sequence(ts)
        for pop in demography.populations:
            pop.initial_size=100.0
        demography.add_population(
                initial_size=10,
                name='ancestral',
                extra_metadata={"slim_id": ts.num_populations},
        )
        demography.add_population_split(
                time=ts.metadata["SLiM"]["tick"] + 20.0,
                derived=[p.name for p in demography.populations if p.name != "ancestral"],
                ancestral="ancestral",
        )
        recap = self.do_recapitate(
                ts,
                demography=demography,
                recombination_rate=recomb_rate,
                random_seed=333,
        )
        self.check_recap_consistency(ts, recap, with_ancestral_Ne=False)

    @pytest.mark.parametrize('recipe', recipe_eq(exclude="long"), indirect=True)
    def test_first_gen_nodes(self, recipe):
        # check that all the roots of the trees are present
        # (note this will fail if some populations were started at different
        # times than others)
        ts = recipe["ts"]
        root_time = ts.metadata['SLiM']['tick']
        is_wf = (ts.metadata['SLiM']['model_type'] == 'WF')
        remembered_stage = ts.metadata['SLiM']['stage']
        if (not is_wf) or (remembered_stage != 'late'):
            root_time -= 1
        if (not is_wf) and ("begun_first" in recipe):
            root_time += 1
        if (not is_wf) and ("remembered_first" in recipe):
            root_time -= 1
        if is_wf and ("begun_late" in recipe):
            root_time -= 1
        for t in ts.trees():
            for u in t.roots:
                assert ts.node(u).time == root_time


class TestIndividualAges(tests.PyslimTestCase):
    # tests for individuals_alive_at and individual_ages_at

    @pytest.mark.parametrize('recipe', [next(recipe_eq("everyone"))], indirect=True)
    def test_errors(self, recipe):
        ts = recipe["ts"]
        for stage in ['abcd', 10, []]:
            with pytest.raises(ValueError):
                pyslim.individuals_alive_at(ts, 0, stage=stage)
            with pytest.raises(ValueError):
                pyslim.individuals_alive_at(ts, 0, remembered_stage=stage)
            with pytest.raises(ValueError):
                pyslim.individual_ages_at(ts, 0, stage=stage)

    @pytest.mark.parametrize('recipe', [next(recipe_eq("pedigree", "WF"))], indirect=True)
    def test_mismatched_remembered_stage(self, recipe):
        ts = recipe["ts"]
        if "remembered_early" in recipe:
            with pytest.warns(UserWarning):
                pyslim.individuals_alive_at(ts, 0, remembered_stage="late")
        else:
            with pytest.warns(UserWarning):
                pyslim.individuals_alive_at(ts, 0, remembered_stage="early")

    @pytest.mark.parametrize('recipe', recipe_eq("multipop", exclude="remembered_early"), indirect=True)
    def test_population(self, recipe):
        ts = recipe["ts"]
        individual_populations = ts.individuals_population
        all_inds = pyslim.individuals_alive_at(ts, 0)
        assert len(all_inds) > 0
        for p in range(ts.num_populations):
            sub_inds = pyslim.individuals_alive_at(ts, 0, population=p)
            assert set(sub_inds) == set(all_inds[individual_populations == p])
            sub_inds = pyslim.individuals_alive_at(ts, 0, population=[p])
            assert set(sub_inds) == set(all_inds[individual_populations == p])
        sub_inds = pyslim.individuals_alive_at(ts, 0, population=np.arange(p))
        assert set(sub_inds) == set(all_inds[individual_populations != p])

    @pytest.mark.parametrize('recipe', recipe_eq("nonWF", exclude="remembered_early"), indirect=True)
    def test_samples_only(self, recipe):
        ts = recipe["ts"]
        all_inds = pyslim.individuals_alive_at(ts, 0)
        assert set(all_inds) == set(pyslim.individuals_alive_at(ts, 0, samples_only=False))
        sub_inds = np.random.choice(all_inds, size=min(len(all_inds), 4), replace=False)
        flags = np.array([n.flags & (tskit.NODE_IS_SAMPLE * n.individual in sub_inds)
                          for n in ts.nodes()], dtype=np.uint32)
        tables = ts.dump_tables()
        tables.nodes.flags = flags
        new_ts = tables.tree_sequence()
        assert set(sub_inds) == set(pyslim.individuals_alive_at(new_ts, 0, samples_only=True))

    @pytest.mark.parametrize('recipe', recipe_eq(exclude=("remembered_early", "long")), indirect=True)
    def test_after_simplify(self, recipe):
        ts = recipe["ts"]
        sts = ts.simplify()
        orig_inds = pyslim.individuals_alive_at(ts, 0)
        simp_inds = pyslim.individuals_alive_at(sts, 0)
        odict = {ts.individual(i).metadata["pedigree_id"]: i for i in orig_inds}
        sdict = {sts.individual(i).metadata["pedigree_id"]: i for i in simp_inds}
        for slim_id in odict:
            i = odict[slim_id]
            ind = ts.individual(i)
            n = ts.node(ind.nodes[0])
            if n.flags & tskit.NODE_IS_SAMPLE:
                assert slim_id in sdict

    @pytest.mark.parametrize('recipe', recipe_eq("pedigree"), indirect=True)
    def test_ages(self, recipe):
        ts = recipe["ts"]
        info = recipe["info"]
        remembered_stage = 'late'
        if 'remembered_first' in recipe:
            remembered_stage = 'first' 
        elif 'remembered_early' in recipe:
            remembered_stage = 'early' 
        assert remembered_stage == ts.metadata['SLiM']['stage']
        max_time_ago = ts.metadata['SLiM']['tick']
        if remembered_stage in ('first', 'early'):
            max_time_ago -= 1
        for time in range(0, max_time_ago):
            slim_tick = ts.metadata['SLiM']['tick'] - time
            check_stages = ('first', 'early', 'late')
            if time == 0:
                if remembered_stage == 'first':
                    # if we remember in first we don't know who's still there
                    # in later stages of that time step
                    check_stages = ('first', )
                elif remembered_stage == 'early':
                    check_stages = ('first', 'early')
            if time == max_time_ago:
                if remembered_stage == 'early':
                    # if we set up the population in early there aren't individuals in first
                    # of the very first time step
                    check_stages = ('early', 'late')
                elif remembered_stage == 'late':
                    # similarly for late
                    check_stages = ('late', )
            for stage in check_stages:
                alive = pyslim.individuals_alive_at(
                            ts,
                            time,
                            stage=stage,
                            remembered_stage=remembered_stage
                )
                ages = pyslim.individual_ages_at(
                            ts,
                            time,
                            stage=stage,
                            remembered_stage=remembered_stage
                )
                for ind in ts.individuals():
                    ind_time = ts.node(ind.nodes[0]).time
                    # bad things can happen for the very first individuals
                    # depending on when the subpops are created
                    if (('everyone' in recipe or ind_time == 0) and
                            (remembered_stage == "early" or ind_time < max_time_ago)):
                        slim_id = ind.metadata["pedigree_id"]
                        assert slim_id in info
                        slim_alive = (slim_tick, stage) in info[slim_id]['age']
                        pyslim_alive = ind.id in alive
                        assert slim_alive == pyslim_alive
                        if slim_alive:
                            slim_age = info[slim_id]['age'][(slim_tick, stage)]
                            if ts.metadata["SLiM"]["model_type"] == "WF":
                                # SLiM records -1 but we return 0 in late and 1 in early
                                slim_age = 0 + (stage in ('first', 'early'))
                            assert ages[ind.id] == slim_age
                        else:
                            assert np.isnan(ages[ind.id])


class TestHasIndividualParents(tests.PyslimTestCase):

    def verify_has_parents(self, ts):
        right_answer = np.repeat(True, ts.num_individuals)
        node_indivs = ts.tables.nodes.individual
        parent_ids = [set() for _ in ts.individuals()]
        node_parent_ids = [set() for _ in ts.nodes()]
        individual_times = ts.individuals_time
        individual_ages = pyslim.individual_ages(ts)
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
                            ptime = individual_times[p]
                            parent_alive = True
                            if ts.metadata["SLiM"]["model_type"] == "WF":
                                if individual_times[i.id] + 1 != ptime:
                                    parent_alive = False
                            else:
                                pdeath = ptime - individual_ages[p]
                                if individual_times[i.id] + 1 < pdeath:
                                    parent_alive = False
                            if not parent_alive:
                                right_answer[i.id] = False
                            else:
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
        right_parents = []
        for j, p in enumerate(parent_ids):
            if right_answer[j]:
                for pp in p:
                    right_parents.append([pp, j])
        has_parents = pyslim.has_individual_parents(ts)
        right_parents = np.sort(np.array(right_parents), axis=0)
        parents = np.sort(pyslim.individual_parents(ts), axis=0)
        assert np.array_equal(right_answer, has_parents)
        assert np.array_equal(right_parents, parents)

    def get_first_gen(self, ts):
        # root_time = ts.metadata["SLiM"]["tick"]
        # if ts.metadata['SLiM']['model_type'] != 'WF' or ts.metadata['SLiM']['stage'] != 'late':
        #     root_time -= 1
        nodes = ts.tables.nodes
        root_time = np.max(nodes.time)
        first_gen = set(nodes.individual[nodes.time == root_time])
        first_gen.discard(tskit.NULL)
        return np.array(list(first_gen), dtype='int')

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_everyone(self, recipe):
        # since everyone is recorded, only the initial individuals should
        # not have parents
        ts = recipe["ts"]
        right_answer = np.repeat(True, ts.num_individuals)
        first_gen = self.get_first_gen(ts)
        assert len(first_gen) > 0
        right_answer[first_gen] = False
        has_parents = pyslim.has_individual_parents(ts)
        assert np.array_equal(right_answer, has_parents)
        self.verify_has_parents(ts)

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_post_recap(self, recipe):
        # the same should be true after recapitation
        ts = recipe["ts"]
        right_answer = np.repeat(True, ts.num_individuals)
        first_gen = self.get_first_gen(ts)
        right_answer[first_gen] = False
        assert(ts.num_populations <= 2)
        ts = self.do_recapitate(ts, recombination_rate=0.01, ancestral_Ne=10, random_seed=11)
        has_parents = pyslim.has_individual_parents(ts)
        assert np.array_equal(right_answer, has_parents)
        self.verify_has_parents(ts)

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_post_simplify(self, recipe):
        ts = recipe["ts"]
        rng = np.random.default_rng(seed=3)
        individual_times = ts.individuals_time
        keep_indivs = rng.choice(
                np.where(individual_times < ts.metadata['SLiM']['tick'] - 1)[0],
                size=30,
                replace=False
        )
        keep_nodes = []
        for i in keep_indivs:
            keep_nodes.extend(ts.individual(i).nodes)
        ts = ts.simplify(samples=keep_nodes, filter_individuals=True, keep_input_roots=True)
        assert(ts.num_populations <= 2)
        ts = self.do_recapitate(ts, recombination_rate=0.01, ancestral_Ne=10)
        has_parents = pyslim.has_individual_parents(ts)
        assert sum(has_parents) > 0
        self.verify_has_parents(ts)

    @pytest.mark.parametrize('recipe', recipe_eq("everyone"), indirect=True)
    def test_pedigree_parents_everyone(self, recipe):
        # We can only guarantee to correctly reconstruct parents when everyone is remembered:
        # for instance, A selfs to produce B who selfs to produce C; if A and C are present
        # but B is not, and A is still alive, we will think that A is C's parent.
        # Or, suppose that X is the parent of Y, but we did not remember X at or after
        # the time that Y was born, so that although X is alive at Y's birth, we don't know it,
        # and so the parentage would not be reported by `individual_parents()`.
        ts = recipe["ts"]
        info = recipe["info"]
        has_parents = pyslim.has_individual_parents(ts)
        parents = pyslim.individual_parents(ts)
        slim_map = {}
        for ind in ts.individuals():
            slim_map[ind.metadata["pedigree_id"]] = ind.id
        ts_to_slim = {sid: [] for sid in slim_map}
        for (pa, ch) in parents:
            assert pa >= 0 and pa < ts.num_individuals
            assert ch >= 0 and pa < ts.num_individuals
            pa_ind = ts.individual(pa).metadata["pedigree_id"]
            ch_ind = ts.individual(ch).metadata["pedigree_id"]
            ts_to_slim[ch_ind].append(pa_ind)
        for hasp, ind in zip(has_parents, ts.individuals()):
            for n in ind.nodes:
                assert ts.node(n).is_sample()
            assert len(ind.nodes) == 2
            sid = ind.metadata["pedigree_id"]
            ts_p = ts_to_slim[sid]
            assert hasp == (len(ts_p) > 0)
            # parents, as recorded by SLiM, that we know about:
            slim_p = [x for x in info[sid]["parents"] if x in slim_map]
            assert set(slim_p) == set(ts_p)

    @pytest.mark.parametrize('recipe', recipe_eq("pedigree"), indirect=True)
    def test_pedigree_parents(self, recipe):
        # Less strict test for consistency only: see caveats above in test_pedigree_parents_everyone.
        # In particular, we are only guaranteed to have whole genomes for ALIVE
        # or REMEMBERED individuals (not RETAINED), and the same for parental genomes.
        ts = recipe["ts"]
        info = recipe["info"]
        has_parents = pyslim.has_individual_parents(ts)
        parents = pyslim.individual_parents(ts)
        slim_map = {}
        for ind in ts.individuals():
            slim_map[ind.metadata["pedigree_id"]] = ind.id
        ts_to_slim = {sid: [] for sid in slim_map}
        for (pa, ch) in parents:
            assert pa >= 0 and pa < ts.num_individuals
            assert ch >= 0 and pa < ts.num_individuals
            pa_ind = ts.individual(pa).metadata["pedigree_id"]
            ch_ind = ts.individual(ch).metadata["pedigree_id"]
            ts_to_slim[ch_ind].append(pa_ind)
        for hasp, ind in zip(has_parents, ts.individuals()):
            all_there = (ind.flags & (pyslim.INDIVIDUAL_ALIVE | pyslim.INDIVIDUAL_REMEMBERED) > 0)
            if all_there:
                for n in ind.nodes:
                    assert ts.node(n).is_sample()
                assert len(ind.nodes) == 2
            sid = ind.metadata["pedigree_id"]
            ts_p = ts_to_slim[sid]
            assert hasp == (len(ts_p) > 0)
            # parents, as recorded by SLiM, that we know about:
            slim_p = [x for x in info[sid]["parents"] if x in slim_map]
            # all pyslim parents should be legit (so set(ts_p) - set(slim_p) should usually be empty)
            # BUT sometimes we can mistake a (great)^n-grandparent for a parent
            gfolks = []
            for a in set(info[sid]["parents"]) - set(ts_p):
                gfolks.extend(info[a]["parents"])
            for a in set(ts_p) - set(slim_p):
                assert a in gfolks


class TestReferenceSequence(tests.PyslimTestCase):
    '''
    Test for operations involving the reference sequence
    '''

    @pytest.mark.parametrize('recipe', recipe_eq(exclude="long"), indirect=True)
    def test_reference_sequence(self, recipe):
        ts = recipe["ts"]
        if ts.num_mutations > 0:
            mut_md = ts.mutation(0).metadata
            has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
            if not has_nucleotides:
                assert not ts.has_reference_sequence()
            else:
                assert type(ts.reference_sequence.data) == type('')
                assert len(ts.reference_sequence.data) == ts.sequence_length
                for u in ts.reference_sequence.data:
                    assert u in pyslim.NUCLEOTIDES
            sts = ts.simplify(ts.samples()[:2])
            assert sts.has_reference_sequence() == ts.has_reference_sequence()
            if sts.has_reference_sequence():
                assert sts.reference_sequence.data == ts.reference_sequence.data

    def test_mutation_at_errors(self, recipe):
        ts = recipe["ts"]
        u = ts.samples()[0]
        with pytest.raises(ValueError):
            pyslim.mutation_at(ts, -2, 3)
        with pytest.raises(ValueError):
            pyslim.mutation_at(ts, u, -3)
        with pytest.raises(ValueError):
            pyslim.mutation_at(ts, ts.num_nodes + 2, 3)
        with pytest.raises(ValueError):
            pyslim.mutation_at(ts, u, ts.sequence_length)

    def test_nucleotide_at_errors(self, recipe):
        ts = recipe["ts"]
        u = ts.samples()[0]
        if ts.num_mutations > 0:
            mut_md = ts.mutation(0).metadata
            has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
            if not has_nucleotides:
                with pytest.raises(ValueError):
                    pyslim.nucleotide_at(ts, u, 3)

    def test_mutation_at(self, recipe):
        random.seed(42)
        ts = recipe["ts"]
        for _ in range(100):
            node = random.randint(0, ts.num_nodes - 1)
            pos = random.randint(0, ts.sequence_length - 1)
            tree = ts.at(pos)
            parent = tree.parent(node)
            a = pyslim.mutation_at(ts, node, pos)
            if parent == tskit.NULL:
                assert a == tskit.NULL
            else:
                b = pyslim.mutation_at(ts, parent, pos)
                c = pyslim.mutation_at(ts, node, pos, ts.node(parent).time)
                assert b == c
                for k in np.where(node == ts.tables.mutations.node)[0]:
                    mut = ts.mutation(k)
                    if ts.site(mut.site).position == pos:
                        b = mut.id
                assert a == b

    def test_nucleotide_at(self, recipe):
        random.seed(42)
        ts = recipe["ts"]
        if ts.num_mutations > 0:
            mut_md = ts.mutation(0).metadata
            has_nucleotides = (mut_md["mutation_list"][0]["nucleotide"] >= 0)
            if has_nucleotides:
                assert ts.has_reference_sequence()
                assert len(ts.reference_sequence.data) == ts.sequence_length
                for _ in range(100):
                    node = random.randint(0, ts.num_nodes - 1)
                    pos = random.randint(0, ts.sequence_length - 1)
                    tree = ts.at(pos)
                    parent = tree.parent(node)
                    a = pyslim.nucleotide_at(ts, node, pos)
                    if parent == tskit.NULL:
                        nuc = ts.reference_sequence.data[int(pos)]
                        assert a == pyslim.NUCLEOTIDES.index(nuc)
                    else:
                        b = pyslim.nucleotide_at(ts, parent, pos)
                        c = pyslim.nucleotide_at(ts, node, pos, ts.node(parent).time)
                        assert b == c
                        for k in np.where(node == ts.tables.mutations.node)[0]:
                            mut = ts.mutation(k)
                            if ts.site(mut.site).position == pos:
                                b = mut.metadata["mutation_list"][0]["nucleotide"]
                        assert a == b

    @pytest.mark.parametrize('recipe', recipe_eq("mutation_spectrum"), indirect=True)
    def test_nucleotide_spectrum(self, recipe):
        # this is modified from Recipe 18.13
        # Also note that we are comparing to "truth" in a recipe where
        # truth is recorded during a mutation callback, in which we have
        # access to the parental genome, so if two adjacent mutations
        # occur in the same meiosis then each will not know about the other.
        ts = recipe["ts"]
        mutation_spectrum = recipe["mutation_info"]
        M = {
                a + b + c + "," + d : 0
                for a in pyslim.NUCLEOTIDES
                for b in pyslim.NUCLEOTIDES
                for c in pyslim.NUCLEOTIDES
                for d in pyslim.NUCLEOTIDES
        }
        nmuts = 0
        for mut in ts.mutations():
            pos = ts.site(mut.site).position 
            if pos > 0 and pos < ts.sequence_length - 1:
                nmuts += 1
                mut_list = mut.metadata["mutation_list"]
                k = np.argmax([u["slim_time"] for u in mut_list])
                derived_nuc = mut_list[k]["nucleotide"]
                left_nuc = pyslim.nucleotide_at(ts, mut.node, pos - 1, time = mut.time + 1.0)
                right_nuc = pyslim.nucleotide_at(ts, mut.node, pos + 1, time = mut.time + 1.0)
                parent_nuc = pyslim.nucleotide_at(ts, mut.node, pos, time = mut.time + 1.0)
                context = "".join([pyslim.NUCLEOTIDES[k] for k in (left_nuc, parent_nuc, right_nuc)])
                key = context + "," + pyslim.NUCLEOTIDES[derived_nuc]
                M[key] += 1
                if key == "ACA,T" or key == "CCA,T":
                    print(key, pos, mut.node, mut.time)
        assert sum([M[k] for k in M]) == nmuts
        assert sum([mutation_spectrum[k][0] for k in mutation_spectrum]) == nmuts
        for k in M:
            if M[k] != mutation_spectrum[k][0]:
                print(k, M[k], mutation_spectrum[k])
        for k in M:
            assert len(mutation_spectrum[k]) == 1
            assert M[k] == mutation_spectrum[k][0]


class TestConvertNucleotides(tests.PyslimTestCase):
    '''
    Test for operations involving the converting and generating nucleotides
    '''

    def last_slim_mutations(self, ts):
        # iterator over mutations, returning for each mutation in ts a tuple
        # (slim id, slim mutation metadata) of the slim mutation that is the
        # *most recent* one of any possibly stacked mutations. Note that it
        # is possible that this is ambiguous.
        for mut in ts.mutations():
            slim_muts = {
                    k : v for k, v in zip(
                        mut.derived_state.split(","),
                        mut.metadata['mutation_list']
                    )
            }
            if mut.parent == tskit.NULL:
                parent_slim_ids = []
            else:
                parent_mut = ts.mutation(mut.parent)
                parent_slim_ids = parent_mut.derived_state.split(",")
            max_time = max([md['slim_time'] for md in slim_muts.values()])
            any_new = any([k not in parent_slim_ids for k in slim_muts.keys()
                           if slim_muts[k]['slim_time'] == max_time])
            maybe_these = [k for k in slim_muts.keys()
                           if slim_muts[k]["slim_time"] == max_time
                           and ((k not in parent_slim_ids)
                                or (not any_new))
                           ]
            k = max(maybe_these)
            yield k, slim_muts[k]

    def verify_converted_nucleotides(self, ts, cts):
        assert ts.has_reference_sequence() == cts.has_reference_sequence()
        if ts.has_reference_sequence():
            assert ts.reference_sequence.data == cts.reference_sequence.data
        assert ts.num_sites == cts.num_sites
        for k, (s, ns) in enumerate(zip(ts.sites(), cts.sites())):
            assert s.position == ns.position
            assert s.metadata == ns.metadata
            assert ns.ancestral_state == ts.reference_sequence.data[int(s.position)]
        for m, cm, (_, sm) in zip(ts.mutations(), cts.mutations(), self.last_slim_mutations(ts)):
            assert m.site == cm.site
            assert m.node == cm.node
            assert m.parent == cm.parent
            assert m.time == cm.time
            assert m.metadata == cm.metadata
            nuc = sm['nucleotide']
            assert nuc in [0, 1, 2, 3]
            assert cm.derived_state == pyslim.NUCLEOTIDES[nuc]
        # should not have changed anything else
        tc = ts.dump_tables()
        ntc = cts.dump_tables()
        tc.sites.clear()
        ntc.sites.clear()
        tc.mutations.clear()
        ntc.mutations.clear()
        tc.provenances.clear()
        ntc.provenances.clear()
        assert tc == ntc

    def scramble_mutations(self, ts):
        # scramble order of mutations so that the most recent is not always first,
        # since we don't have a reliable way to get that out of SLiM
        rng = np.random.default_rng(123)
        t = ts.dump_tables()
        t.mutations.clear()
        for m in ts.mutations():
            a = np.array(m.derived_state.split(","))
            ii = rng.permutation(len(a))
            ml = [m.metadata['mutation_list'][i] for i in ii]
            t.mutations.append(
                    m.replace(
                        derived_state=",".join(a[ii]),
                        metadata={'mutation_list': ml}
                    )
            )
        t.compute_mutation_parents()
        return t.tree_sequence()

    def test_convert_alleles_errors(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, population_size=10)
        with pytest.raises(ValueError, match="must have a valid reference sequence"):
            _ = pyslim.convert_alleles(ts)
        ts = pyslim.annotate(ts, model_type="nonWF", tick=1)
        with pytest.raises(ValueError, match="must have a valid reference sequence"):
            _ = pyslim.convert_alleles(ts)
        mts = msprime.sim_mutations(ts,
                model=msprime.SLiMMutationModel(type=1),
                rate=0.1,
                random_seed=23)
        assert mts.num_mutations > 0
        mtt = mts.dump_tables()
        mtt.reference_sequence.data = 'A' * int(mts.sequence_length)
        mts = mtt.tree_sequence()
        with pytest.raises(ValueError, match="must be nucleotide mutations"):
            _ = pyslim.convert_alleles(mts)

    @pytest.mark.parametrize(
            'recipe', recipe_eq("nucleotides", exclude="non-nucleotides"), indirect=True
    )
    def test_convert_alleles(self, recipe):
        ts = recipe["ts"]
        cts = pyslim.convert_alleles(ts)
        self.verify_converted_nucleotides(ts, cts)

        # get some weirder situations in there
        t = ts.dump_tables()
        t.mutations.clear()
        for j, mut in enumerate(ts.mutations()):
            if j % 2 == 0:
                t.mutations.append(mut)
        t.compute_mutation_parents()
        ts = t.tree_sequence()
        cts = pyslim.convert_alleles(ts)
        self.verify_converted_nucleotides(ts, cts)

    @pytest.mark.parametrize(
            'recipe', recipe_eq("nucleotides", exclude="non-nucleotides"), indirect=True
    )
    def test_convert_alleles_scrambled(self, recipe):
        ts = self.scramble_mutations(recipe["ts"])
        cts = pyslim.convert_alleles(ts)
        self.verify_converted_nucleotides(ts, cts)

    @pytest.mark.parametrize(
            'recipe', recipe_eq("nucleotides", exclude="non-nucleotides"), indirect=True
    )
    def test_keeps_reference_sequence(self, recipe):
        ts = recipe["ts"]
        assert ts.has_reference_sequence()
        nts = pyslim.generate_nucleotides(ts, seed=123)
        assert nts.has_reference_sequence()
        assert ts.reference_sequence == nts.reference_sequence

    def test_generate_nucleotides_errors(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, population_size=10, random_seed=777)
        with pytest.raises(ValueError, match="must have length equal"):
            _ = pyslim.generate_nucleotides(ts, reference_sequence="AAA")
        with pytest.raises(ValueError, match="must have length equal"):
            _ = pyslim.generate_nucleotides(ts, reference_sequence=[1, 2, 3])
        with pytest.raises(ValueError, match="must be a string of"):
            _ = pyslim.generate_nucleotides(
                    ts,
                    reference_sequence="X" * int(ts.sequence_length)
            )
        with pytest.raises(ValueError, match="must be a string of"):
            _ = pyslim.generate_nucleotides(
                    ts,
                    reference_sequence=np.arange(int(ts.sequence_length)),
            )

    def verify_generate_nucleotides(self, ts, check_transitions=False):
        # if check_transitions is True, verify that derived states differ
        # from parental states - which we try to do but is not guaranteed,
        # for instance, if keep=True or in other weird situations.
        assert len(ts.reference_sequence.data) == ts.sequence_length
        muts = {}
        ts_muts = {
            j : v['nucleotide']
            for j, (_, v) in enumerate(self.last_slim_mutations(ts))
        }
        for mut in ts.mutations():
            aa = ts.reference_sequence.data[int(ts.site(mut.site).position)]
            for i, md in zip(
                    mut.derived_state.split(","),
                    mut.metadata['mutation_list']
                 ):
                nuc = md['nucleotide']
                assert nuc in [0, 1, 2, 3]
                if i in muts:
                    assert muts[i] == nuc
                muts[i] = nuc
            if check_transitions:
                if mut.parent == tskit.NULL:
                    assert pyslim.NUCLEOTIDES[nuc] != aa
                else:
                    if ts.mutation(mut.parent).derived_state != mut.derived_state:
                        assert ts_muts[mut.parent] != ts_muts[mut.id]

    def test_generate_nucleotides(self, recipe):
        ts = recipe["ts"]
        nts = pyslim.generate_nucleotides(ts, keep=False, seed=5)
        self.verify_generate_nucleotides(
                nts,
                check_transitions=("adds_mutations" not in recipe),
        )

    def test_generate_nucleotides_refseq(self):
        ts = msprime.sim_ancestry(
                4,
                sequence_length=10,
                population_size=10,
                random_seed=10,
        )
        ts = pyslim.annotate(ts, model_type='nonWF', tick=1)
        mts = msprime.sim_mutations(ts,
                model=msprime.SLiMMutationModel(type=1),
                rate=0.5,
                random_seed=23)
        refseq = "A" * int(mts.sequence_length)
        nts = pyslim.generate_nucleotides(mts, reference_sequence=refseq, seed=6)
        self.verify_generate_nucleotides(nts, check_transitions=True)
        assert nts.reference_sequence.data == refseq

    def test_generate_nucleotides_keep(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, population_size=10)
        ts = pyslim.annotate(ts, model_type='nonWF', tick=1)
        mts1 = msprime.sim_mutations(ts,
                model=msprime.SLiMMutationModel(type=1),
                rate=0.1,
                random_seed=23)
        mts1.dump("out.trees")
        nts1 = pyslim.generate_nucleotides(mts1, seed=10, keep=False)
        assert nts1.num_mutations > 0
        self.verify_generate_nucleotides(nts1, check_transitions=False)
        mts2 = msprime.sim_mutations(nts1,
                model=msprime.SLiMMutationModel(
                    type=2,
                    next_id=nts1.num_mutations,
                ),
                rate=0.1,
                random_seed=24,
        )
        # keep defaults to True
        nts2 = pyslim.generate_nucleotides(mts2, seed=12)
        assert nts2.num_mutations > nts1.num_mutations
        muts1 = {}
        for mut in nts1.mutations():
            for i, md in zip(mut.derived_state.split(","), mut.metadata['mutation_list']):
                muts1[i] = md['nucleotide']
        for mut in nts2.mutations():
            for i, md in zip(mut.derived_state.split(","), mut.metadata['mutation_list']):
                if md['mutation_type'] == 1:
                    assert i in muts1
                    assert muts1[i] == md['nucleotide']
                else:
                    assert md['nucleotide'] in [0, 1, 2, 3]
        nts3 = pyslim.generate_nucleotides(mts2, keep=False, seed=15)
        self.verify_generate_nucleotides(nts3, check_transitions=False)

    @pytest.mark.parametrize(
            'recipe', ["recipe_long_nonWF.slim"], indirect=True
    )
    def test_generate_and_convert(self, recipe, helper_functions, tmp_path):
        ts = recipe["ts"]
        nts = pyslim.generate_nucleotides(ts, seed=123)
        cts = pyslim.convert_alleles(nts)
        self.verify_converted_nucleotides(nts, cts)
        helper_functions.run_slim_restart(
                nts,
                "restart_nucleotides_nonWF.slim",
                tmp_path,
                WF=False,
        )


class TestDeprecations(tests.PyslimTestCase):
    # test on one arbitrary recipe
    @pytest.mark.skip(reason="TODO")
    @pytest.mark.parametrize('recipe', [next(recipe_eq())], indirect=True)
    def test_slim_tree_sequence(self, recipe):
        ts = recipe["ts"]
        with pytest.warns(FutureWarning):
            _ = pyslim.SlimTreeSequence(ts)

