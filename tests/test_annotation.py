"""
Test cases for the metadata reading/writing of pyslim.
"""
import json
import numpy as np
import pytest
import random
import warnings

import msprime
import tskit

import pyslim
import tests
from .recipe_specs import restarted_recipe_eq

class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests for tools to annotate existing msprime-derived tree sequences.
    '''

    def verify_annotated_tables(self, ts1, ts2, check_alleles=True):
        '''
        Verify that the tables returned after annotation are equal, up to the
        expected forgetting of metadata, flags, and other things.
        '''
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
                assert md["is_null"] is False
                assert md["genome_type"] == pyslim.GENOME_TYPE_AUTOSOME
        for ind in ts.individuals():
            md = ind.metadata
            assert np.array_equal(ind.location, [0, 0, 0])
            assert ind.flags == pyslim.INDIVIDUAL_ALIVE
            assert md["sex"] == pyslim.INDIVIDUAL_TYPE_HERMAPHRODITE
            assert md["flags"] == 0
            assert md["pedigree_p1"] == tskit.NULL
            assert md["pedigree_p2"] == tskit.NULL
        for pop in ts.populations():
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

    def verify_slim_restart_equality(self, in_ts, out_ts):
        """
        Check for equality, in everything but the last provenance.
        """
        assert in_ts.num_provenances + 1 == out_ts.num_provenances
        in_tables = in_ts.dump_tables()
        in_tables.sort()
        out_tables = out_ts.dump_tables()
        out_tables.sort()
        self.assertTableCollectionsEqual(in_tables, out_tables, skip_provenance=-1)

    def test_annotate_errors(self, helper_functions):
        for ts in helper_functions.get_msprime_examples():
            with pytest.raises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="WF",
                                             slim_generation=0)
            with pytest.raises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="WF",
                                             slim_generation=4.4)
            with pytest.raises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="foo",
                                             slim_generation=4)
            with pytest.raises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type=[],
                                             slim_generation=4)
        # odd number of samples
        ts = msprime.simulate(3)
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate_defaults(ts, model_type="WF",
                                         slim_generation=1)
        assert "diploid" in str(except_info)
        # inconsistent populations for diploids
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(sample_size=3),
                msprime.PopulationConfiguration(sample_size=1)],
            migration_matrix=[[0.0, 1.0], [1.0, 0.0]]
            )
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate_defaults(ts, model_type="WF",
                                         slim_generation=1)
        assert "more than one population" in str(except_info)
        # inconsistent times for diploids
        samples = [
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=0, time=1),
        ]
        ts = msprime.simulate(samples=samples)
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate_defaults(ts, model_type="WF",
                                         slim_generation=1)
        assert "more than one time" in str(except_info)
        # non-integer positions
        ts = msprime.simulate(4, mutation_rate=10)
        with pytest.raises(ValueError) as except_info:
            _ = pyslim.annotate_defaults(ts, model_type="WF",
                                         slim_generation=1)
        assert "not at integer" in str(except_info)

    def test_just_simulate(self, helper_functions, tmp_path):
        ts = msprime.simulate(
                sample_size=4, Ne=10, length=10,
                mutation_rate=0.0, recombination_rate=0.01
        )
        ts = msprime.sim_mutations(ts, rate=0.1)
        slim_ts = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1) 
        loaded_ts = helper_functions.run_msprime_restart(slim_ts, tmp_path, WF=True)
        self.verify_annotated_trees(ts, loaded_ts)

    def test_basic_annotation(self, helper_functions, tmp_path):
        for ts in helper_functions.get_msprime_examples():
            for do_mutations in [False, True]:
                slim_gen = 4
                slim_ts = pyslim.annotate_defaults(
                        ts, model_type="WF",
                        slim_generation=slim_gen,
                        annotate_mutations=do_mutations,
                )
                assert slim_ts.metadata['SLiM']['model_type'] == 'WF'
                assert slim_ts.metadata['SLiM']['generation'] == slim_gen
                assert slim_ts.metadata['SLiM']['file_version'] == pyslim.slim_file_version
                self.verify_annotated_tables(ts, slim_ts, check_alleles=(not do_mutations))
                self.verify_annotated_trees(ts, slim_ts)
                if not do_mutations:
                    self.verify_haplotype_equality(ts, slim_ts)
                self.verify_defaults(slim_ts)
                self.verify_provenance(slim_ts)
                # try loading this into SLiM
                loaded_ts = helper_functions.run_msprime_restart(slim_ts, tmp_path, WF=True)
                self.verify_annotated_tables(loaded_ts, slim_ts)
                self.verify_annotated_trees(loaded_ts, slim_ts)
                self.verify_haplotype_equality(loaded_ts, slim_ts)

    def test_annotate_refseq(self):
        ts = msprime.sim_ancestry(2, sequence_length=10, random_seed=77)
        refseq = "A" * int(ts.sequence_length)
        ats = pyslim.annotate_defaults(
                ts, model_type="nonWF", slim_generation=5, reference_sequence=refseq
        )
        assert ats.reference_sequence.data == refseq

    def test_annotate_individuals(self, helper_functions, tmp_path):
        for ts in helper_functions.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
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
                    [ims.validate_and_encode_row(r) for r in metadata])
            pop_metadata = [p.metadata for p in tables.populations]
            for j, md in enumerate(pop_metadata):
                # nonWF models always have this
                md['sex_ratio'] = 0.0
            pms = tables.populations.metadata_schema
            tables.populations.packset_metadata(
                    [pms.validate_and_encode_row(r) for r in pop_metadata])
            new_ts = pyslim.load_tables(tables)
            for j, ind in enumerate(new_ts.individuals()):
                md = ind.metadata
                assert md["sex"] == sexes[j]
            self.verify_annotated_tables(new_ts, slim_ts)
            self.verify_annotated_trees(new_ts, slim_ts)
            self.verify_haplotype_equality(new_ts, slim_ts)
            # try loading this into SLiM
            loaded_ts = helper_functions.run_msprime_restart(new_ts, tmp_path, sex="A")
            self.verify_trees_equal(new_ts, loaded_ts)

    @pytest.mark.skip(reason="not working")
    def test_annotate_XY(self, helper_functions, tmp_path):
        # This is Not Right because as written the null chromosomes have history;
        # we need to *only* simulate the non-null chromosomes.
        random.seed(8)
        for ts in helper_functions.get_msprime_examples():
            for genome_type in ["X", "Y"]:
                slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
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
                        [ims.validate_and_encode_row(r) for r in metadata])
                node_metadata = [n.metadata for n in tables.nodes]
                node_is_null = [False for _ in range(ts.num_nodes)]
                for j in range(slim_ts.num_individuals):
                    nodes = slim_ts.individual(j).nodes
                    node_metadata[nodes[0]]["genome_type"] = pyslim.GENOME_TYPE_X
                    node_metadata[nodes[0]]["is_null"] = (genome_type != "X")
                    if sexes[j] == pyslim.INDIVIDUAL_TYPE_MALE:
                        node_metadata[nodes[1]]["genome_type"] = pyslim.GENOME_TYPE_Y
                        node_metadata[nodes[1]]["is_null"] = (genome_type != "Y")
                    else:
                        node_metadata[nodes[1]]["genome_type"] = pyslim.GENOME_TYPE_X
                        node_metadata[nodes[1]]["is_null"] = (genome_type != "X")
                    for n in nodes:
                        node_is_null[n] = node_metadata[n]["is_null"]
                nms = tables.nodes.metadata_schema
                tables.nodes.packset_metadata(
                        [nms.validate_and_encode_row(r) for r in node_metadata])
                # update populations
                pop_metadata = [p.metadata for p in tables.populations]
                for j, md in enumerate(pop_metadata):
                    # nonWF models always have this
                    md['sex_ratio'] = 0.0
                pms = tables.populations.metadata_schema
                tables.populations.packset_metadata(
                        [pms.validate_and_encode_row(r) for r in pop_metadata])
                new_ts = tables.tree_sequence()
                self.verify_annotated_tables(new_ts, slim_ts)
                self.verify_annotated_trees(new_ts, slim_ts)
                self.verify_haplotype_equality(new_ts, slim_ts)
                # try loading this into SLiM
                loaded_ts = helper_functions.run_msprime_restart(
                    new_ts, tmp_path, sex=genome_type)
                self.verify_trees_equal(new_ts, loaded_ts)
                # these are *not* equal but only due to re-ordering of nodes and individuals
                # ... and for some reason, .subset( ) or .simplify( ) do not produce equality
                # self.assertTableCollectionsEqual(new_ts, loaded_ts,
                #         skip_provenance=-1, reordered_individuals=True)

    def test_annotate_nodes(self, helper_functions):
        for ts in helper_functions.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.dump_tables()
            metadata = [n.metadata for n in tables.nodes]
            gtypes = [random.choice([pyslim.GENOME_TYPE_X, pyslim.GENOME_TYPE_Y])
                      for _ in metadata]
            for md, g in zip(metadata, gtypes):
                if md is not None:
                    md["genome_type"] = g
            nms = tables.nodes.metadata_schema
            tables.nodes.packset_metadata(
                    [nms.validate_and_encode_row(r) for r in metadata])
            new_ts = pyslim.load_tables(tables)
            for x, g in zip(new_ts.nodes(), gtypes):
                if x.metadata is not None:
                    assert x.metadata["genome_type"] == g
            # not testing SLiM because needs annotation of indivs to make sense

    def test_annotate_mutations(self, helper_functions):
        for ts in helper_functions.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.dump_tables()
            metadata = [m.metadata for m in tables.mutations]
            selcoefs = [random.uniform(0, 1) for _ in metadata]
            for j in range(len(metadata)):
                metadata[j]['mutation_list'][0]["selection_coeff"] = selcoefs[j]
            ms = tables.mutations.metadata_schema
            tables.mutations.packset_metadata(
                    [ms.validate_and_encode_row(r) for r in metadata])
            new_ts = pyslim.load_tables(tables)
            for j, x in enumerate(new_ts.mutations()):
                md = x.metadata
                assert np.isclose(md['mutation_list'][0]["selection_coeff"], selcoefs[j])

    def test_dont_annotate_mutations(self, helper_functions):
        # Test the option to not overwrite mutation annotations
        ts = msprime.sim_ancestry(10)
        ts = msprime.sim_mutations(ts, rate=5, random_seed=3)
        assert ts.num_mutations > 0
        tables = ts.dump_tables()
        pre_mutations = tables.mutations.copy()
        pyslim.annotate_defaults_tables(tables, model_type="WF",
                slim_generation=1, annotate_mutations=False)
        # this is necessary because b'' actually is decoded to
        # an empty mutation_list by the schema
        pre_mutations.metadata_schema = tables.mutations.metadata_schema
        assert tables.mutations.equals(pre_mutations)

    def test_empty_populations(self, helper_functions, tmp_path):
        # test SLiM doesn't error on having empty populations
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
        ts = pyslim.annotate_defaults(ts, model_type='WF', slim_generation=1)
        t = ts.dump_tables()
        for k in range(5):
            md = pyslim.default_slim_metadata('population')
            md['name'] = f"new_pop_num_{k}"
            md['description'] = f"the {k}-th added pop"
            t.populations.add_row(metadata=md)
        ts = t.tree_sequence()
        sts = helper_functions.run_slim_restart(
                ts,
                "restart_WF.slim",
                tmp_path,
                WF=True,
        )

    def test_many_populations(self, helper_functions, tmp_path):
        # test we can add more than one population
        ts = msprime.sim_ancestry(5, population_size=10, sequence_length=100, random_seed=455)
        t = ts.dump_tables()
        for k in range(5):
            md = pyslim.default_slim_metadata('population')
            md['name'] = f"new_pop_num_{k}"
            md['description'] = f"the {k}-th added pop"
            t.populations.add_row(metadata=md)
            i = t.individuals.add_row()
            for _ in range(2):
                t.nodes.add_row(flags=1, time=0.0, individual=i, population=k)
        ts = t.tree_sequence()
        ts = pyslim.annotate_defaults(ts, model_type='WF', slim_generation=1)
        for ind in ts.individuals():
            assert ind.flags == pyslim.INDIVIDUAL_ALIVE
        sts = helper_functions.run_slim_restart(
                ts,
                "restart_WF.slim",
                tmp_path,
                WF=True,
        )

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_reload_recapitate(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        # Test the ability of SLiM to load our files after recapitation.
        ts = recipe["ts"]
        # recapitate, reload
        in_ts = self.do_recapitate(ts, recombination_rate=1e-2, ancestral_Ne=10, random_seed=25)
        # put it through SLiM (which just reads in and writes out)
        out_ts = helper_functions.run_slim_restart(in_ts, restart_name, tmp_path)
        # check for equality, in everything but the last provenance
        in_ts.dump("in_ts.trees")
        out_ts.dump("out_ts.trees")
        self.verify_slim_restart_equality(in_ts, out_ts)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_reload_annotate(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        # Test the ability of SLiM to load our files after annotation.
        ts = recipe["ts"]
        tables = ts.dump_tables()
        metadata = [m.metadata for m in tables.mutations]
        has_nucleotides = tables.metadata['SLiM']['nucleotide_based']
        if has_nucleotides:
            nucs = [random.choice([0, 1, 2, 3]) for _ in metadata]
            refseq = "".join(random.choices(pyslim.NUCLEOTIDES,
                                            k = int(ts.sequence_length)))
            for n, md in zip(nucs, metadata):
                for m in md['mutation_list']:
                    m["nucleotide"] = n
            tables.reference_sequence.data = refseq
        for md in metadata:
            for m in md['mutation_list']:
                m["selection_coeff"] = random.random()
        ms = tables.mutations.metadata_schema
        tables.mutations.packset_metadata(
                [ms.validate_and_encode_row(r) for r in metadata])
        in_ts = pyslim.load_tables(tables)
        # put it through SLiM (which just reads in and writes out)
        out_ts = helper_functions.run_slim_restart(in_ts, restart_name, tmp_path)
        # check for equality, in everything but the last provenance
        self.verify_slim_restart_equality(in_ts, out_ts)


class TestReload(tests.PyslimTestCase):
    '''
    Tests for basic things related to reloading with SLiM
    '''

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op"), indirect=["recipe"])
    def test_load_without_provenance(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        in_ts = recipe["ts"]
        in_tables = in_ts.dump_tables()
        in_tables.provenances.clear()
        in_tables.sort()
        cleared_ts = pyslim.SlimTreeSequence(
                in_tables.tree_sequence(),
        )
        out_ts = helper_functions.run_slim_restart(cleared_ts, restart_name, tmp_path)
        out_tables = out_ts.dump_tables()
        out_tables.provenances.clear()
        out_tables.sort()
        in_tables.assert_equals(out_tables)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq("no_op", "nucleotides"), indirect=["recipe"])
    def test_reload_reference_sequence(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        in_ts = recipe["ts"]
        out_ts = helper_functions.run_slim_restart(in_ts, restart_name, tmp_path)
        assert in_ts.metadata['SLiM']['nucleotide_based'] is True
        assert out_ts.metadata['SLiM']['nucleotide_based'] is True
        assert in_ts.has_reference_sequence() == out_ts.has_reference_sequence()
        if in_ts.has_reference_sequence():
            in_ts.reference_sequence.assert_equals(out_ts.reference_sequence)

    @pytest.mark.parametrize(
        'restart_name, recipe', restarted_recipe_eq(), indirect=["recipe"])
    def test_restarts_and_runs_spatial(
        self, restart_name, recipe, helper_functions, tmp_path
    ):
        in_ts = recipe["ts"]
        n = set(in_ts.samples()[:4])
        for k in n:
            n |= set(in_ts.individual(in_ts.node(k).individual).nodes)
        py_ts = in_ts.simplify(list(n))
        out_ts = helper_functions.run_slim_restart(py_ts, restart_name, tmp_path)
        assert out_ts.metadata["SLiM"]["generation"] >= in_ts.metadata["SLiM"]["generation"]

