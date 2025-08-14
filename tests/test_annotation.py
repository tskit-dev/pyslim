"""
Test cases for the metadata reading/writing of pyslim.
"""
import json
import numpy as np
import pytest
import random
import warnings
import contextlib

import msprime
import tskit

import pyslim
import tests
from .recipe_specs import restarted_recipe_eq

def mutcontext(ts):
    if ts.num_mutations > 0:
        handler = pytest.warns(Warning, match="already has.*metadata")
    else:
        handler = contextlib.nullcontext()
    return handler

def verify_slim_restart_equality(in_ts_dict, out_ts_dict, check_prov=True):
    """
    Check for equality, in everything but the last provenance.
    """
    assert in_ts_dict.keys() == out_ts_dict.keys()
    for k in in_ts_dict:
        in_ts = in_ts_dict[k]
        out_ts = out_ts_dict[k]
        if check_prov:
            assert in_ts.num_provenances + 1 == out_ts.num_provenances
        in_tables = in_ts.dump_tables()
        in_tables.sort()
        out_tables = out_ts.dump_tables()
        out_tables.sort()
        in_tables.assert_equals(out_tables, ignore_provenance=True)

class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests for tools to annotate existing msprime-derived tree sequences.
    '''

    def verify_annotated_tables(self, ts1, ts2, check_alleles=True):
        '''
        Verify that the tables returned after annotation are equal, up to the
        expected forgetting of metadata, flags, and other things.
        '''
        assert pyslim.is_current_version(ts2)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        # compare nodes
        assert np.array_equal(tables1.nodes.flags, tables2.nodes.flags)
        assert np.array_equal(tables1.nodes.time, tables2.nodes.time)
        assert np.array_equal(tables1.nodes.population, tables2.nodes.population)
        # compare edges
        assert tables1.edges == tables2.edges
        # compare sites
        assert np.array_equal(tables1.sites.position, tables2.sites.position)
        if check_alleles:
            assert np.array_equal(tables1.sites.ancestral_state, tables2.sites.ancestral_state)
            assert np.array_equal(tables1.sites.ancestral_state_offset,
                                  tables2.sites.ancestral_state_offset)
        # compare mutations
        assert np.array_equal(tables1.mutations.site, tables2.mutations.site)
        assert np.array_equal(tables1.mutations.node, tables2.mutations.node)
        assert np.array_equal(tables1.mutations.time, tables2.mutations.time)
        if check_alleles:
            assert np.array_equal(tables1.mutations.derived_state, tables2.mutations.derived_state)
            assert np.array_equal(tables1.mutations.derived_state_offset,
                                  tables2.mutations.derived_state_offset)

    def verify_annotated_trees(self, ts1, ts2):
        '''
        Verify the *trees* returned before and after annotation are equal.
        '''
        assert ts1.num_trees == ts2.num_trees
        for t1, t2 in zip(ts1.trees(), ts2.trees()):
            assert t1.length == t2.length
            assert t1.get_parent_dict() == t2.get_parent_dict()
            assert t1.total_branch_length == t2.total_branch_length

    def verify_defaults(self, ts):
        '''
        Verify the default values have been entered into metadata.
        '''
        do_pops = [False for _ in ts.populations()]
        for m in ts.mutations():
            md = m.metadata
            assert isinstance(md['mutation_list'], list)
            for mdl in md['mutation_list']:
                assert mdl["mutation_type"] == 0
                assert mdl["selection_coeff"] == 0.0
                assert mdl["subpopulation"] == tskit.NULL
        for n in ts.nodes():
            md = n.metadata
            if not n.is_sample():
                assert md is None
            else:
                assert md["is_vacant"] == [0]
                do_pops[n.population] = True
        for ind in ts.individuals():
            md = ind.metadata
            assert np.array_equal(ind.location, [0, 0, 0])
            assert ind.flags == pyslim.INDIVIDUAL_ALIVE
            assert md["sex"] == pyslim.INDIVIDUAL_TYPE_HERMAPHRODITE
            assert md["flags"] == 0
            assert md["pedigree_p1"] == tskit.NULL
            assert md["pedigree_p2"] == tskit.NULL
        for pop in ts.populations():
            if do_pops[pop.id]:
                md = pop.metadata
                assert md["selfing_fraction"] == 0.0
                assert md["female_cloning_fraction"] == 0.0
                assert md["male_cloning_fraction"] == 0.0
                assert md["sex_ratio"] == 0.0
                assert md["bounds_x0"] == 0.0
                assert md["bounds_x1"] == 1.0
                assert md["bounds_y0"] == 0.0
                assert md["bounds_y1"] == 1.0
                assert md["bounds_z0"] == 0.0
                assert md["bounds_z1"] == 1.0
                assert len(md["migration_records"]) == 0

    def verify_provenance(self, ts):
        for u in ts.provenances():
            tskit.validate_provenance(json.loads(u.record))

    def verify_remapping(self, ts, rts, subpop_map):
        # below we assume the restart has not done additional simulation
        do_remap = (subpop_map is not None)
        if not do_remap:
            subpop_map = { f"p{k}" : k for k in range(ts.num_populations) }
        fwd_map = {-1 : -1}
        rev_map = {-1 : -1}
        for pk in subpop_map:
            # keys are requried to be of the form "pK"
            rts_k = int(pk[1:])
            assert rts_k < rts.num_populations
            ts_k = subpop_map[pk]
            assert ts_k < ts.num_populations
            assert ts_k not in fwd_map
            fwd_map[ts_k] = rts_k
            assert rts_k not in rev_map
            rev_map[rts_k] = ts_k

        # all pops in ts should be present UNLESS they have no metadata
        for pop in ts.populations():
            assert pop.id in fwd_map or pop.metadata is None

        # any extra pops in rts should have empty metadata
        for pop in rts.populations():
            if pop.id not in rev_map:
                assert pop.metadata is None

        for ts_k, rts_k in fwd_map.items():
            if ts_k != -1:
                ts_pop = ts.population(ts_k)
                rts_pop = rts.population(rts_k)
                ts_md = ts_pop.metadata
                rts_md = rts_pop.metadata
                if ts_md is None:
                    assert rts_md is None
                else:
                    if do_remap and (
                            ts_md['name'] == f"p{ts_pop.id}"
                            or ts_md['name'] == f"pop_{ts_pop.id}"
                    ):
                        ts_md['name'] = f"p{rts_pop.id}"
                    assert ts_md['name'] == rts_md['name']
                    if "slim_id" not in ts_md:
                        assert ts_md == rts_md
                    else:
                        assert ts_md["slim_id"] == ts_pop.id
                        assert rts_md["slim_id"] == rts_pop.id
                        k = "migration_records"
                        if k in ts_md:
                            assert k in rts_md or len(ts_md[k]) == 0
                        else:
                            assert k not in rts_md or len(rts_md[k]) == 0
                        if k in ts_md and k in rts_md:
                            assert len(ts_md[k]) == len(rts_md[k])
                            for x, y in zip(ts_md[k], rts_md[k]):
                                assert x['migration_rate'] == y['migration_rate']
                                assert fwd_map[x['source_subpop']] == y['source_subpop']

        # check other tables have been remapped
        # (uses the fact that no additional simulation has been done)
        assert ts.num_nodes == rts.num_nodes
        for n, rn in zip(ts.nodes(), rts.nodes()):
            n.population = fwd_map[n.population]
            # there can be rounding error due to time shift
            assert np.isclose(n.time, rn.time, atol=1e-16)
            n.time = rn.time
            assert n == rn

        assert ts.num_migrations == rts.num_migrations
        for m, rm in zip(ts.migrations(), rts.migrations()):
            m.source = fwd_map[m.source]
            m.dest = fwd_map[m.dest]
            assert m == rm

        assert ts.num_mutations == rts.num_mutations
        for m, rm in zip(ts.mutations(), rts.mutations()):
            md = m.metadata
            rmd = rm.metadata
            assert len(md['mutation_list']) == len(rmd['mutation_list'])
            for x, y in zip(md['mutation_list'], rmd['mutation_list']):
                x['subpopulation'] = fwd_map[x['subpopulation']]
                assert x == y

        assert ts.num_individuals == rts.num_individuals
        for i, ri in zip(ts.individuals(), rts.individuals()):
            md = i.metadata
            rmd = ri.metadata
            md['subpopulation'] = fwd_map[md['subpopulation']]
            assert md == rmd


    def test_annotate_errors(self, helper_functions):
        for ts in helper_functions.get_msprime_examples():
            with pytest.raises(ValueError):
                _ = pyslim.annotate(ts, model_type="WF", tick=0)
            with pytest.raises(ValueError):
                _ = pyslim.annotate(ts, model_type="WF", tick=4.4)
            with pytest.raises(ValueError):
                _ = pyslim.annotate(ts, model_type="foo", tick=4)
            with pytest.raises(ValueError):
                _ = pyslim.annotate(ts, model_type=[], tick=4)
        # no individuals
        ts = msprime.simulate(4)
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate(ts, model_type="WF", tick=1)
        assert "individuals" in str(except_info)
        # non-integer positions
        ts = msprime.sim_mutations(
                msprime.sim_ancestry(4, random_seed=8),
                rate=10,
                discrete_genome=False,
                random_seed=9,
        )
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate(ts, model_type="WF", tick=1)
        assert "not at integer" in str(except_info)

    def test_warns_overwriting_mutations(self, helper_functions):
        ts = msprime.sim_ancestry(
                4,
                population_size=10,
                sequence_length=10,
                recombination_rate=0.01,
                random_seed=100,
        )
        ts = msprime.sim_mutations(
                ts,
                rate=1,
                random_seed=12,
                model=msprime.SLiMMutationModel(type=1)
        )
        assert ts.num_mutations > 0
        with pytest.warns(Warning, match="already has.*metadata"):
            slim_ts = pyslim.annotate(ts, model_type="WF", tick=1) 

    def test_warns_and_overwrites_node_metadata(self, helper_functions):
        ts = msprime.sim_ancestry(
                4,
                population_size=10,
                sequence_length=10,
                recombination_rate=0.01,
                random_seed=100,
        )
        tables = ts.dump_tables()
        tables.nodes.clear()
        for n in ts.nodes():
            tables.nodes.append(n.replace(metadata=b'abc'))
        ts = tables.tree_sequence()
        with pytest.warns(Warning, match="already has.*metadata"):
            slim_ts = pyslim.annotate(ts, model_type="WF", tick=1) 
            for n in slim_ts.nodes():
                assert n.is_sample() or n.metadata is None

    def test_warns_and_overwrites_individual_metadata(self, helper_functions):
        ts = msprime.sim_ancestry(
                4,
                population_size=10,
                sequence_length=10,
                recombination_rate=0.01,
                random_seed=100,
        )
        tables = ts.dump_tables()
        tables.individuals.clear()
        for n in ts.individuals():
            tables.individuals.append(n.replace(metadata=b'abc'))
        ts = tables.tree_sequence()
        with pytest.warns(Warning, match="already has.*metadata"):
            slim_ts = pyslim.annotate(ts, model_type="WF", tick=1) 
            for ind in slim_ts.individuals():
                assert slim_ts.node(ind.nodes[0]).is_sample() or ind.metadata is None

    def test_just_simulate(self, helper_functions, tmp_path):
        ts = msprime.sim_ancestry(
                4,
                population_size=10,
                sequence_length=10,
                recombination_rate=0.01,
                random_seed=100,
        )
        ts = msprime.sim_mutations(ts, rate=0.1, random_seed=11)
        slim_ts = pyslim.annotate(ts, model_type="WF", tick=1) 
        loaded_ts = helper_functions.run_msprime_restart(
                {"default": slim_ts}, tmp_path, multichrom=False, WF=True
        )["default"]
        self.verify_annotated_trees(ts, loaded_ts)

    def test_basic_annotation(self, helper_functions, tmp_path):
        for ts in helper_functions.get_msprime_examples():
            for do_mutations in [False, True]:
                tick = 4
                cycle = 1
                stage = "late"
                if do_mutations:
                    handler = mutcontext(ts)
                else:
                    handler = contextlib.nullcontext()
                with handler:
                    slim_ts = pyslim.annotate(
                            ts, model_type="WF",
                            tick=tick,
                            cycle=cycle,
                            stage=stage,
                            annotate_mutations=do_mutations,
                    )
                assert slim_ts.metadata['SLiM']['model_type'] == 'WF'
                assert slim_ts.metadata['SLiM']['tick'] == tick
                assert slim_ts.metadata['SLiM']['cycle'] == cycle
                assert slim_ts.metadata['SLiM']['stage'] == stage
                assert slim_ts.metadata['SLiM']['file_version'] == pyslim.slim_file_version
                self.verify_annotated_tables(ts, slim_ts, check_alleles=(not do_mutations))
                self.verify_annotated_trees(ts, slim_ts)
                if not do_mutations:
                    self.verify_haplotype_equality(ts, slim_ts)
                self.verify_defaults(slim_ts)
                self.verify_provenance(slim_ts)
                # try loading this into SLiM
                loaded_ts = helper_functions.run_msprime_restart(
                        {"default": slim_ts}, tmp_path, multichrom=False, WF=True
                )["default"]
                self.verify_annotated_tables(loaded_ts, slim_ts)
                self.verify_annotated_trees(loaded_ts, slim_ts)
                self.verify_haplotype_equality(loaded_ts, slim_ts)

    def test_annotate_refseq(self):
        ts = msprime.sim_ancestry(2, sequence_length=10, random_seed=77)
        refseq = "A" * int(ts.sequence_length)
        ats = pyslim.annotate(
                ts, model_type="nonWF", tick=5, reference_sequence=refseq
        )
        assert ats.reference_sequence.data == refseq

    def test_annotate_individuals(self, helper_functions, tmp_path):
        # test the workflow of annotating defaults, then assigning sexes randomly
        for ts in helper_functions.get_msprime_examples():
            with mutcontext(ts):
                slim_ts = pyslim.annotate(ts, model_type="nonWF", tick=1, stage="early")
            tables = slim_ts.dump_tables()
            top_md = tables.metadata
            top_md['SLiM']['separate_sexes'] = True
            tables.metadata = top_md
            metadata = [ind.metadata for ind in tables.individuals]
            sexes = [random.choice([pyslim.INDIVIDUAL_TYPE_FEMALE, pyslim.INDIVIDUAL_TYPE_MALE])
                     for _ in metadata]
            for j in range(len(metadata)):
                metadata[j]["sex"] = sexes[j]
            ims = tables.individuals.metadata_schema
            tables.individuals.packset_metadata(
                    [ims.validate_and_encode_row(r) for r in metadata]
            )
            pop_metadata = [p.metadata for p in tables.populations]
            for j, md in enumerate(pop_metadata):
                # nonWF models always have this
                md['sex_ratio'] = 0.0
            pms = tables.populations.metadata_schema
            tables.populations.packset_metadata(
                    [pms.validate_and_encode_row(r) for r in pop_metadata]
            )
            new_ts = tables.tree_sequence()
            for j, ind in enumerate(new_ts.individuals()):
                md = ind.metadata
                assert md["sex"] == sexes[j]
            self.verify_annotated_tables(new_ts, slim_ts)
            self.verify_annotated_trees(new_ts, slim_ts)
            self.verify_haplotype_equality(new_ts, slim_ts)
            # try loading this into SLiM
            loaded_ts = helper_functions.run_msprime_restart(
                    {"default": new_ts}, tmp_path, multichrom=False, sex="A"
            )["default"]
            self.verify_trees_equal(new_ts, loaded_ts)

    def test_annotate_nodes(self, helper_functions):
        # test workflow of annotating defaults and then editing node metadata
        for ts in helper_functions.get_msprime_examples():
            with mutcontext(ts):
                slim_ts = pyslim.annotate(ts, model_type="nonWF", tick=1)
            tables = slim_ts.dump_tables()
            metadata = [n.metadata for n in tables.nodes]
            sids = list(range(ts.num_nodes))
            random.shuffle(sids)
            for md, sid in zip(metadata, sids):
                if md is not None:
                    md["slim_id"] = sid
            nms = tables.nodes.metadata_schema
            tables.nodes.packset_metadata(
                    [nms.validate_and_encode_row(r) for r in metadata]
            )
            new_ts = tables.tree_sequence()
            for x, sid in zip(new_ts.nodes(), sids):
                if x.metadata is not None:
                    assert x.metadata["slim_id"] == sid
            # not testing SLiM because needs annotation of indivs to make sense

    def test_annotate_mutations(self, helper_functions):
        for ts in helper_functions.get_msprime_examples():
            with mutcontext(ts):
                slim_ts = pyslim.annotate(ts, model_type="nonWF", tick=1)
            tables = slim_ts.dump_tables()
            metadata = [m.metadata for m in tables.mutations]
            selcoefs = [random.uniform(0, 1) for _ in metadata]
            for j in range(len(metadata)):
                metadata[j]['mutation_list'][0]["selection_coeff"] = selcoefs[j]
            ms = tables.mutations.metadata_schema
            tables.mutations.packset_metadata(
                    [ms.validate_and_encode_row(r) for r in metadata]
            )
            new_ts = tables.tree_sequence()
            for j, x in enumerate(new_ts.mutations()):
                md = x.metadata
                assert np.isclose(md['mutation_list'][0]["selection_coeff"], selcoefs[j])

    def test_dont_annotate_mutations(self):
        # Test the option to not overwrite mutation annotations
        ts = msprime.sim_ancestry(10, random_seed=5)
        ts = msprime.sim_mutations(ts, rate=5, random_seed=3)
        assert ts.num_mutations > 0
        tables = ts.dump_tables()
        pre_mutations = tables.mutations.copy()
        pyslim.annotate_tables(
                tables,
                model_type="WF",
                tick=1,
                annotate_mutations=False
        )
        # this is necessary because b'' actually is decoded to
        # an empty mutation_list by the schema
        pre_mutations.metadata_schema = tables.mutations.metadata_schema
        assert tables.mutations.equals(pre_mutations)

    def test_annotate_empty(self):
        t = tskit.TableCollection(sequence_length=10)
        at = pyslim.annotate_tables(t, model_type='nonWF', tick=1)
        assert t.metadata_schema == pyslim.slim_metadata_schemas['tree_sequence']
        for x in ("population", "individual", "node", "site", "mutation"):
            assert getattr(t, x + "s").metadata_schema == pyslim.slim_metadata_schemas[x]

    def test_no_populations_errors(self, helper_functions, tmp_path):
        # test SLiM errors when individuals are in non-SLiM populations
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=456)
        t = ts.dump_tables()
        p = t.populations.add_row(metadata={'name': '', 'description': ''})
        i = t.individuals.add_row()
        t.nodes.add_row(time=0.0, flags=tskit.NODE_IS_SAMPLE, individual=i)
        t.nodes.add_row(time=0.0, flags=tskit.NODE_IS_SAMPLE, individual=i)
        pyslim.annotate_tables(t, model_type='nonWF', tick=1, stage="late")
        pop = t.populations[-1]
        t.populations[-1] = pop.replace(metadata = {'name' : 'not_slim', 'description' : ''})
        ind = t.individuals[-1]
        md = ind.metadata
        md['subpopulation'] = (t.populations.num_rows - 1)
        t.individuals[-1] = ind.replace(metadata=md)
        ts = t.tree_sequence()
        with pytest.raises(RuntimeError, match='individual has a subpopulation id .1.'):
            _ = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    WF=False,
                    multichrom=False,
            )

    def test_empty_populations_errors(self, helper_functions, tmp_path):
        # test SLiM errors on having empty populations in a WF model
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
        ts = pyslim.annotate(ts, model_type='WF', tick=1, stage="late")
        # first make sure does not error without the empty pops
        rts = helper_functions.run_slim_restart(
                {"default": ts},
                "restart_WF.slim",
                tmp_path,
                multichrom=False,
                WF=True,
        )["default"]
        assert rts.num_populations == 1
        # now add empty subpops
        t = ts.dump_tables()
        for k in range(5):
            md = pyslim.default_slim_metadata('population')
            md['slim_id'] = k + 1
            md['name'] = f"new_pop_num_{k}"
            t.populations.add_row(metadata=md)
        ts = t.tree_sequence()
        for p in ts.populations():
            assert "slim_id" in p.metadata
        with pytest.raises(RuntimeError, match='subpopulation.*empty'):
            _ = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_WF.slim",
                    tmp_path,
                    WF=True,
                    multichrom=False,
            )

    def test_empty_populations(self, helper_functions, tmp_path):
        # test SLiM doesn't error on having empty populations in a nonWF model
        # however, it may rewrite the metadata
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=456)
        ts = pyslim.annotate(ts, model_type='nonWF', tick=1)
        t = ts.dump_tables()
        for k in range(5):
            md = pyslim.default_slim_metadata('population')
            md['slim_id'] = k + 1
            md['name'] = f"new_pop_num_{k}"
            md['description'] = f"the {k}-th added pop"
            t.populations.add_row(metadata=md)
        ts = t.tree_sequence()
        for subpop_map in (
                None,
                {"p0": 0, "p1": 1, "p2": 2, "p3": 3, "p4": 4, "p5": 5},
                {"p0": 0, "p2": 1, "p1": 2, "p9": 3, "p4": 4, "p6": 5},
                {"p10": 0, "p2": 1, "p1": 2, "p9": 3, "p4": 4, "p6": 5},
            ):
            sts = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    WF=False,
                    subpop_map=subpop_map,
                    multichrom=False,
            )["default"]
            self.verify_remapping(ts, sts, subpop_map)

    def test_many_populations(self, helper_functions, tmp_path):
        # test we can add more than one population and that names are preserved
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
        t = ts.dump_tables()
        for k in range(1, 6):
            md = pyslim.default_slim_metadata('population')
            md['name'] = f"new_pop_num_{k}"
            md['description'] = f"the {k}-th added pop"
            md['slim_id'] = k
            t.populations.add_row(metadata=md)
            i = t.individuals.add_row()
            for _ in range(2):
                t.nodes.add_row(flags=1, time=0.0, individual=i, population=k)
        ts = t.tree_sequence()
        ts = pyslim.annotate(ts, model_type='WF', tick=1)
        for ind in ts.individuals():
            assert ind.flags == pyslim.INDIVIDUAL_ALIVE
        sts = helper_functions.run_slim_restart(
                {"default": ts},
                "restart_WF.slim",
                tmp_path,
                WF=True,
                multichrom=False,
        )["default"]
        for k in range(1, 6):
            md = sts.population(k).metadata
            assert md['name'] == f"new_pop_num_{k}"

    def test_many_species_errors(self, helper_functions, tmp_path):
        # should error when there are population conflicts between species
        fox_ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
        fox_ts = pyslim.annotate(fox_ts, model_type='WF', tick=1)
        mouse_ts = msprime.sim_ancestry(5, population_size=20, sequence_length=1000, random_seed=234)
        mouse_ts = pyslim.annotate(mouse_ts, model_type='WF', tick=1)
        with pytest.raises(RuntimeError, match='p0 has been used already'):
            _ = helper_functions.run_multispecies_slim_restart(
                    {'fox': fox_ts, 'mouse': mouse_ts},
                    'restart_multispecies_WF.slim',
                    tmp_path,
                    WF=True,
            )
        with pytest.raises(RuntimeError, match='p1 has been used already'):
            _ = helper_functions.run_multispecies_slim_restart(
                    {'fox': fox_ts, 'mouse': mouse_ts},
                    'restart_multispecies_WF.slim',
                    tmp_path,
                    WF=True,
                    subpop_map={'fox': {'p1': 0}, 'mouse': {'p1': 0}},
            )

    def test_many_species(self, helper_functions, tmp_path):
        # verify we can load more than one species
        for wf in (True, False):
            fox_ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
            mouse_ts = msprime.sim_ancestry(5, population_size=20, sequence_length=1000, random_seed=234)
            fox_ts = pyslim.annotate(fox_ts, model_type="WF" if wf else "nonWF", tick=1)
            mouse_ts = pyslim.annotate(mouse_ts, model_type="WF" if wf else "nonWF", tick=1)
            for subpop_map in (
                        {'fox' : None, 'mouse' : {'p1' : 0}},
                        {'fox' : {'p1' : 0}, 'mouse' : None},
                        {'fox' : {'p3' : 0}, 'mouse' : {'p0' : 0}},
                    ):
                tsd = helper_functions.run_multispecies_slim_restart(
                        {'fox': fox_ts, 'mouse': mouse_ts},
                        'restart_multispecies_WF.slim' if wf else 'restart_multispecies_nonWF.slim',
                        tmp_path,
                        WF=wf,
                        subpop_map = subpop_map,
                )
                self.verify_remapping(fox_ts, tsd['fox'], subpop_map['fox'])
                self.verify_remapping(mouse_ts, tsd['mouse'], subpop_map['mouse'])

    def test_keeps_extra_populations(self, helper_functions, tmp_path):
        # test that non-SLiM metadata is not removed
        d = msprime.Demography()
        d.add_population(initial_size=10, name="A", extra_metadata={"pi": 3.1415})
        md = pyslim.default_slim_metadata('population')
        del md['name']
        del md['description']
        md['slim_id'] = 1
        d.add_population(initial_size=10, name="B", extra_metadata=md)
        d.add_population(initial_size=10, name="C", extra_metadata={"abc": 123})
        ts = msprime.sim_ancestry(
                samples={"B": 5},
                demography=d,
                sequence_length=100,
                random_seed=455
        )
        ts = pyslim.annotate(ts, model_type='WF', tick=1)
        sts = helper_functions.run_slim_restart(
                {"default": ts},
                "restart_WF.slim",
                tmp_path,
                multichrom=False,
                WF=True,
        )["default"]
        assert sts.num_populations == 3
        assert np.all(ts.tables.nodes.population == 1)
        p = sts.population(0)
        assert p.metadata['name'] == "A"
        assert p.metadata['pi'] == 3.1415
        p = sts.population(1)
        assert p.metadata['name'] == "B"
        p = sts.population(2)
        assert p.metadata['name'] == "C"
        assert p.metadata['abc'] == 123

    def test_remapping_errors(self, helper_functions, tmp_path):
        d = msprime.Demography()
        d.add_population(initial_size=10, name="p0")
        d.add_population(initial_size=10, name="p1")
        d.add_population(initial_size=10, name="A")
        d.add_population_split(time=1.0, derived=["p0", "p1"], ancestral="A")
        ts = msprime.sim_ancestry(
                {"p0": 5, "p1": 5},
                demography=d,
                sequence_length=100,
                random_seed=455,
        )
        ts = pyslim.annotate(ts, model_type='nonWF', tick=1, stage="late")
        # test must refer to all subpopulations
        with pytest.raises(RuntimeError, match='subpopulation id 2 is used.*but is not remapped'):
            _ = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    multichrom=False,
                    WF=False,
                    subpop_map = {
                        "p0" : 0,
                        "p1" : 1,
                    }
            )
        # test can't map the same subpop to two different pops
        with pytest.raises(RuntimeError, match='subpopMap value .0. is not unique'):
            _ = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    multichrom=False,
                    WF=False,
                    subpop_map = {
                        "p0" : 0,
                        "p1" : 0,
                        "p2" : 1,
                        "p3" : 2,
                    }
            )
        # test errors if it refers to non-existing subpop
        with pytest.raises(RuntimeError):
            _ = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    multichrom=False,
                    WF=False,
                    subpop_map = {
                        "p0" : 0,
                        "p1" : 1,
                        "p2" : 2,
                        "p3" : 3,
                    }
            )

    def test_remapping_empty_pops(self, helper_functions, tmp_path):
        # test that subpop_map doesn't need to explicitly refer to empty pops
        d = msprime.Demography()
        d.add_population(initial_size=10, name="p0")
        d.add_population(initial_size=10, name="p1")
        d.add_population(initial_size=10, name="p2") # unused
        d.add_population(initial_size=10, name="A")
        d.add_population_split(time=1.0, derived=["p0", "p1"], ancestral="A")
        ts = msprime.sim_ancestry(
                {"p0": 5, "p1": 5},
                demography=d,
                sequence_length=100,
                random_seed=456,
        )
        ts = pyslim.annotate(ts, model_type='nonWF', tick=1, stage="late")
        t = ts.dump_tables()
        r = t.populations[2].replace(metadata=None)
        # make both pop 2 and 4 empty
        t.populations[2] = r
        t.populations.append(r)
        ts = t.tree_sequence()
        for subpop_map in (
                {"p0" : 0, "p1" : 1, "p2" : 2, "p3" : 3, "p4" : 4},
                    {"p0" : 0, "p1" : 1, "p3": 3},
                ):
            rts = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_nonWF.slim",
                    tmp_path,
                    multichrom=False,
                    WF=False,
                    subpop_map = subpop_map
            )["default"]
            self.verify_remapping(ts, rts, subpop_map)

    def test_remapping(self, helper_functions, tmp_path):
        d = msprime.Demography()
        # pop 0
        d.add_population(initial_size=10, name="A", extra_metadata={"pi": 3.1415})
        # pop 1
        md = pyslim.default_slim_metadata('population')
        del md['name']
        del md['description']
        md['slim_id'] = 1
        md['migration_records'] = [
            {
                "migration_rate" : 0.1,
                "source_subpop" : 3,
            }
        ]
        d.add_population(initial_size=10, name="pop_1", extra_metadata=md)
        # pop 2
        d.add_population(initial_size=10, name="pop_2", extra_metadata={"abc": 123})
        # pop 3
        md = pyslim.default_slim_metadata('population')
        del md['name']
        del md['description']
        md['slim_id'] = 3
        md['migration_records'] = [
            {
                "migration_rate" : 0.9,
                "source_subpop" : 1,
            }
        ]
        d.add_population(initial_size=10, name="D", extra_metadata=md)
        # pop 4
        d.add_population(initial_size=10, name="XYZ", extra_metadata={"abc": 5123})
        d.add_population_split(time=1.0, derived=["pop_1", "D"], ancestral="A")
        ts = msprime.sim_ancestry(
                samples={"pop_1": 5, "D": 5},
                demography=d,
                sequence_length=100,
                random_seed=455
        )
        ts = pyslim.annotate(ts, model_type='WF', tick=1)
        ts = msprime.sim_mutations(
                ts,
                rate=1e-2,
                random_seed=9,
                model=msprime.SLiMMutationModel(type=1),
        )
        assert ts.num_mutations > 0
        for subpop_map in (
                None,
                {"p0" : 0, "p1" : 1, "p2" : 2, "p3" : 3, "p4" : 4},
                {"p0" : 1, "p1" : 0, "p2" : 2, "p3" : 3, "p4" : 4},
                {"p5" : 0, "p0" : 1, "p7" : 2, "p2" : 3, "p9" : 4},
        ):
            rts = helper_functions.run_slim_restart(
                    {"default": ts},
                    "restart_WF.slim",
                    tmp_path,
                    multichrom=False,
                    WF=True,
                    subpop_map = subpop_map,
            )["default"]
            self.verify_remapping(ts, rts, subpop_map)

    @pytest.mark.skip(reason="migration table not supported by simplify")
    def test_remaps_migration_table(self, helper_functions, tmp_path):
        d = msprime.Demography()
        d.add_population(initial_size=10, name="p0")
        d.add_population(initial_size=10, name="p1")
        d.add_population(initial_size=10, name="A")
        d.add_population_split(time=1.0, derived=["p0", "p1"], ancestral="A")
        d.set_migration_rate(0, 1, 0.1)
        d.set_migration_rate(1, 0, 0.1)
        ts = msprime.sim_ancestry(
                {"p0": 5, "p1": 5},
                demography=d,
                sequence_length=100,
                random_seed=455,
                record_migrations=True
        )
        ts = pyslim.annotate(ts, model_type='WF', tick=1, stage="late")
        assert ts.num_migrations > 0
        rts = helper_functions.run_slim_restart(
                {"default": ts},
                "restart_WF.slim",
                tmp_path,
                multichrom=False,
                WF=True,
                subpop_map = {"p5" : 0, "p2" : 1, "p0" : 2 },
        )["default"]
        self.verify_remapping(ts, rts, subpop_map)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("remove_subpop"), indirect=["recipe"])
    def test_keep_slim_names(self, restart_name, recipe, helper_functions, tmp_path):
        # test that SLiM retains names of SLiM populations
        # even when those populations are removed
        in_ts = recipe["ts"]
        # this recipe moves everyone to subpop 0
        out_ts = helper_functions.run_slim_restart(
                in_ts, restart_name, tmp_path, "multichrom" in recipe
        )
        assert in_ts.keys() == out_ts.keys()
        for k in in_ts:
            it = in_ts[k]
            ot = out_ts[k]
            for n in ot.samples(time=0):
                assert ot.node(n).population == 0
            for j in range(1, it.num_populations):
                in_md = it.population(j).metadata
                out_md = ot.population(j).metadata
                assert in_md['name'] == out_md['name']

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_reload_recapitate(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        # Test the ability of SLiM to load our files after recapitation.
        in_ts = {}
        for chrom, ts in recipe["ts"].items():
            # recapitate, reload
            in_ts[chrom] = self.do_recapitate(
                    ts, recombination_rate=1e-2, ancestral_Ne=10, random_seed=25
            )
        # put it through SLiM (which just reads in and writes out)
        out_ts = helper_functions.run_slim_restart(
                in_ts,
                restart_name,
                tmp_path,
                "multichrom" in recipe,
                WF=False,
        )
        # check for equality, in everything but the last provenance
        verify_slim_restart_equality(in_ts, out_ts)
        # test for restarting with a subpop map
        ex_in_ts = list(in_ts.values())[0]
        for j, k in [(1, 7), (5, 7)]:
            subpop_map = {
                    f"p{i * j % k}" : i
                    for i in range(ex_in_ts.num_populations)
            }
            out_ts = helper_functions.run_slim_restart(
                    in_ts,
                    restart_name,
                    tmp_path,
                    "multichrom" in recipe,
                    subpop_map=subpop_map
            )
            assert in_ts.keys() == out_ts.keys()
            for k in in_ts:
                it = in_ts[k]
                ot = out_ts[k]
                self.verify_remapping(it, ot, subpop_map)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_reload_annotate(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        # Test the ability of SLiM to load our files after annotation.
        in_ts = {}
        for chrom, ts in recipe["ts"].items():
            tables = ts.dump_tables()
            metadata = [m.metadata for m in tables.mutations]
            has_nucleotides = tables.metadata['SLiM']['nucleotide_based']
            if has_nucleotides:
                nucs = [random.choice([0, 1, 2, 3]) for _ in metadata]
                refseq = "".join(
                        random.choices(
                            pyslim.NUCLEOTIDES,
                            k = int(ts.sequence_length),
                        ),
                )
                for n, md in zip(nucs, metadata):
                    for m in md['mutation_list']:
                        m["nucleotide"] = n
                tables.reference_sequence.data = refseq
            for md in metadata:
                for m in md['mutation_list']:
                    m["selection_coeff"] = random.random()
            ms = tables.mutations.metadata_schema
            tables.mutations.packset_metadata(
                    [ms.validate_and_encode_row(r) for r in metadata]
            )
            in_ts[chrom] = tables.tree_sequence()
        # put it through SLiM (which just reads in and writes out)
        out_ts = helper_functions.run_slim_restart(
                in_ts, restart_name, tmp_path, "multichrom" in recipe
        )
        # check for equality, in everything but the last provenance
        verify_slim_restart_equality(in_ts, out_ts)


class TestReload(tests.PyslimTestCase):
    '''
    Tests for basic things related to reloading with SLiM
    '''

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_load_without_provenance(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        cleared_ts = {}
        for chrom, in_ts in recipe["ts"].items():
            in_tables = in_ts.dump_tables()
            in_tables.provenances.clear()
            in_tables.sort()
            cleared_ts[chrom] = in_tables.tree_sequence()
        out_ts = helper_functions.run_slim_restart(
                cleared_ts, restart_name, tmp_path, "multichrom" in recipe
        )
        verify_slim_restart_equality(recipe["ts"], cleared_ts, check_prov=False)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op", "nucleotides"), indirect=["recipe"])
    def test_reload_reference_sequence(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        in_ts = recipe["ts"]
        out_ts = helper_functions.run_slim_restart(
                in_ts, restart_name, tmp_path, "multichrom" in recipe
        )
        assert in_ts.keys() == out_ts.keys()
        for k in in_ts:
            it = in_ts[k]
            ot = out_ts[k]
            assert it.metadata['SLiM']['nucleotide_based'] is True
            assert ot.metadata['SLiM']['nucleotide_based'] is True
            assert it.has_reference_sequence() == ot.has_reference_sequence()
            if it.has_reference_sequence():
                it.reference_sequence.assert_equals(ot.reference_sequence)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq(), indirect=["recipe"])
    def test_restarts_and_runs_simplified(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        in_ts = recipe["ts"]
        py_ts = {}
        n = None
        for chrom, ts in in_ts.items():
            if n is None:
                inds = np.where(ts.individuals_flags & pyslim.INDIVIDUAL_ALIVE > 0)[0][:2]
                n = [u for i in inds for u in ts.individual(i).nodes]
                n.sort()
            py_ts[chrom] = ts.simplify(n, filter_populations=False)
        out_ts = helper_functions.run_slim_restart(
                py_ts,
                restart_name,
                tmp_path,
                "multichrom" in recipe
        )
        assert in_ts.keys() == out_ts.keys()
        for k in in_ts:
            it = in_ts[k]
            ot = out_ts[k]
        assert ot.metadata["SLiM"]["tick"] >= it.metadata["SLiM"]["tick"]
