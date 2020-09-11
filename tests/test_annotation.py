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
import random
import json


class TestAnnotate(tests.PyslimTestCase):
    '''
    Tests for tools to annotate existing msprime-derived tree sequences.
    '''

    def verify_annotated_tables(self, ts1, ts2):
        '''
        Verify that the tables returned after annotation are equal, up to the
        expected forgetting of metadata.
        '''
        tables1 = ts1.tables
        tables2 = ts2.tables
        # compare nodes
        self.assertArrayEqual(tables1.nodes.flags, tables2.nodes.flags)
        self.assertArrayAlmostEqual(tables1.nodes.time, tables2.nodes.time)
        self.assertArrayEqual(tables1.nodes.population, tables2.nodes.population)
        # compare edges
        self.assertEqual(tables1.edges, tables2.edges)
        # compare sites
        self.assertArrayEqual(tables1.sites.position, tables2.sites.position)
        self.assertArrayEqual(tables1.sites.ancestral_state, tables2.sites.ancestral_state)
        self.assertArrayEqual(tables1.sites.ancestral_state_offset,
                              tables2.sites.ancestral_state_offset)
        # compare mutations
        self.assertArrayEqual(tables1.mutations.site, tables2.mutations.site)
        self.assertArrayEqual(tables1.mutations.node, tables2.mutations.node)
        self.assertArrayEqual(tables1.mutations.derived_state, tables2.mutations.derived_state)
        self.assertArrayEqual(tables1.mutations.derived_state_offset,
                              tables2.mutations.derived_state_offset)

    def verify_annotated_trees(self, ts1, ts2):
        '''
        Verify the *trees* returned before and after annotation are equal.
        '''
        self.assertEqual(ts1.num_trees, ts2.num_trees)
        for t1, t2 in zip(ts1.trees(), ts2.trees()):
            self.assertEqual(t1.length, t2.length)
            self.assertEqual(t1.get_parent_dict(), t2.get_parent_dict())
            self.assertAlmostEqual(t1.total_branch_length, t2.total_branch_length)

    def verify_defaults(self, ts):
        '''
        Verify the default values have been entered into metadata.
        '''
        for m in ts.mutations():
            md = m.metadata
            self.assertEqual(md["mutation_type"], 1)
            self.assertEqual(md["selection_coeff"], 0.0)
            self.assertEqual(md["population"], tskit.NULL)
            self.assertEqual(md["slim_time"], 0)
        for n in ts.nodes():
            md = n.metadata
            if not n.is_sample():
                self.assertEqual(md, None)
            else:
                self.assertEqual(md["is_null"], False)
                self.assertEqual(md["genome_type"], pyslim.GENOME_TYPE_AUTOSOME)
        for ind in ts.individuals():
            md = ind.metadata
            self.assertArrayEqual(ind.location, [0, 0, 0])
            self.assertEqual(ind.flags, pyslim.INDIVIDUAL_ALIVE)
            self.assertEqual(md["sex"], pyslim.INDIVIDUAL_TYPE_HERMAPHRODITE)
            self.assertEqual(md["flags"], 0)
        for pop in ts.populations():
            md = pop.metadata
            self.assertEqual(md["selfing_fraction"], 0.0)
            self.assertEqual(md["female_cloning_fraction"], 0.0)
            self.assertEqual(md["male_cloning_fraction"], 0.0)
            self.assertEqual(md["sex_ratio"], 0.5)
            self.assertEqual(md["bounds_x0"], 0.0)
            self.assertEqual(md["bounds_x1"], 0.0)
            self.assertEqual(md["bounds_y0"], 0.0)
            self.assertEqual(md["bounds_y1"], 0.0)
            self.assertEqual(md["bounds_z0"], 0.0)
            self.assertEqual(md["bounds_z1"], 0.0)
            self.assertEqual(len(md["migration_records"]), 0)

    def verify_provenance(self, ts):
        for u in ts.provenances():
            tskit.validate_provenance(json.loads(u.record))

    def verify_slim_restart_equality(self, in_ts, out_ts):
        """
        Check for equality, in everything but the last provenance.
        """
        self.assertEqual(in_ts.num_provenances + 1, out_ts.num_provenances)
        in_tables = in_ts.tables
        in_tables.sort()
        out_tables = out_ts.tables
        out_tables.sort()
        self.assertTableCollectionsEqual(in_tables, out_tables, skip_provenance=-1)

    def test_annotate_errors(self):
        for ts in self.get_msprime_examples():
            with self.assertRaises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="WF",
                                             slim_generation=0)
            with self.assertRaises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="WF",
                                             slim_generation=4.4)
            with self.assertRaises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type="foo",
                                             slim_generation=4)
            with self.assertRaises(ValueError):
                _ = pyslim.annotate_defaults(ts, model_type=[],
                                             slim_generation=4)

    def test_basic_annotation(self):
        for ts in self.get_msprime_examples():
            slim_gen = 4
            slim_ts = pyslim.annotate_defaults(ts, model_type="WF",
                                               slim_generation=slim_gen)
            self.assertEqual(slim_ts.metadata['SLiM']['model_type'], 'WF')
            self.assertEqual(slim_ts.metadata['SLiM']['generation'], slim_gen)
            self.assertEqual(slim_ts.metadata['SLiM']['file_version'], pyslim.slim_file_version)
            self.verify_annotated_tables(ts, slim_ts)
            self.verify_annotated_trees(ts, slim_ts)
            self.verify_haplotype_equality(ts, slim_ts)
            self.verify_defaults(slim_ts)
            self.verify_provenance(slim_ts)
            # try loading this into SLiM
            loaded_ts = self.run_msprime_restart(slim_ts, WF=True)
            self.verify_annotated_tables(loaded_ts, slim_ts)
            self.verify_annotated_trees(loaded_ts, slim_ts)
            self.verify_haplotype_equality(loaded_ts, slim_ts)

    def test_annotate_individuals(self):
        for ts in self.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
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
                self.assertEqual(md["sex"], sexes[j])
            self.verify_annotated_tables(new_ts, slim_ts)
            self.verify_annotated_trees(new_ts, slim_ts)
            self.verify_haplotype_equality(new_ts, slim_ts)
            # try loading this into SLiM
            loaded_ts = self.run_msprime_restart(new_ts, sex="A")
            self.verify_trees_equal(new_ts, loaded_ts)

    def test_annotate_XY(self):
        random.seed(8)
        for ts in self.get_msprime_examples():
            for genome_type in ["X", "Y"]:
                slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
                tables = slim_ts.tables
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
                nms = tables.nodes.metadata_schema
                tables.nodes.packset_metadata(
                        [nms.validate_and_encode_row(r) for r in node_metadata])
                pop_metadata = [p.metadata for p in tables.populations]
                for j, md in enumerate(pop_metadata):
                    # nonWF models always have this
                    md['sex_ratio'] = 0.0
                pms = tables.populations.metadata_schema
                tables.populations.packset_metadata(
                        [pms.validate_and_encode_row(r) for r in pop_metadata])
                new_ts = pyslim.load_tables(tables)
                self.verify_annotated_tables(new_ts, slim_ts)
                self.verify_annotated_trees(new_ts, slim_ts)
                self.verify_haplotype_equality(new_ts, slim_ts)
                # try loading this into SLiM
                loaded_ts = self.run_msprime_restart(new_ts, sex=genome_type)
                self.verify_trees_equal(new_ts, loaded_ts)
                # these are *not* equal but only due to re-ordering of nodes and individuals
                # ... and for some reason, .subset( ) or .simplify( ) do not produce equality
                # self.assertTableCollectionsEqual(new_ts, loaded_ts,
                #         skip_provenance=-1, reordered_individuals=True)

    def test_annotate_nodes(self):
        for ts in self.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
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
                    self.assertEqual(x.metadata["genome_type"], g)
            # not testing SLiM because needs annotation of indivs to make sense

    def test_annotate_mutations(self):
        for ts in self.get_msprime_examples():
            slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
            tables = slim_ts.tables
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
                self.assertEqual(md['mutation_list'][0]["selection_coeff"], selcoefs[j])

    def test_reload_recapitate(self):
        # Test the ability of SLiM to load our files after recapitation.
        for ts, basename in self.get_slim_restarts(no_op=True):
            # recapitate, reload
            in_ts = ts.recapitate(recombination_rate=1e-2, Ne=10, random_seed=25)
            # put it through SLiM (which just reads in and writes out)
            out_ts = self.run_slim_restart(in_ts, basename)
            # check for equality, in everything but the last provenance
            self.verify_slim_restart_equality(in_ts, out_ts)

    def test_reload_annotate(self):
        # Test the ability of SLiM to load our files after annotation.
        for ts, basename in self.get_slim_restarts(no_op=True):
            tables = ts.tables
            metadata = [m.metadata for m in tables.mutations]
            has_nucleotides = tables.metadata['SLiM']['nucleotide_based']
            if has_nucleotides:
                nucs = [random.choice([0, 1, 2, 3]) for _ in metadata]
                refseq = "".join(random.choices(pyslim.NUCLEOTIDES,
                                                k = int(ts.sequence_length)))
                for n, md in zip(nucs, metadata):
                    for m in md['mutation_list']:
                        m["nucleotide"] = n
            else:
                refseq = None
            for md in metadata:
                for m in md['mutation_list']:
                    m["selection_coeff"] = random.random()
            ms = tables.mutations.metadata_schema
            tables.mutations.packset_metadata(
                    [ms.validate_and_encode_row(r) for r in metadata])
            in_ts = pyslim.load_tables(tables, reference_sequence=refseq)
            # put it through SLiM (which just reads in and writes out)
            out_ts = self.run_slim_restart(in_ts, basename)
            # check for equality, in everything but the last provenance
            self.verify_slim_restart_equality(in_ts, out_ts)


class TestReload(tests.PyslimTestCase):
    '''
    Tests for basic things related to reloading with SLiM
    '''

    def test_load_without_provenance(self):
        # with 0.5, SLiM should read info from metadata, not provenances
        for in_ts, basename in self.get_slim_restarts(no_op=True):
            in_tables = in_ts.tables
            in_tables.provenances.clear()
            cleared_ts = pyslim.SlimTreeSequence(
                    in_tables.tree_sequence(),
                    reference_sequence=in_ts.reference_sequence
                    )
            out_ts = self.run_slim_restart(cleared_ts, basename)
            out_tables = out_ts.tables
            out_tables.provenances.clear()
            self.assertEqual(in_tables, out_tables)

    def test_reload_reference_sequence(self):
        for in_ts, basename in self.get_slim_restarts(no_op=True, nucleotides=True):
            out_ts = self.run_slim_restart(in_ts, basename)
            self.assertEqual(in_ts.metadata['SLiM']['nucleotide_based'], True)
            self.assertEqual(out_ts.metadata['SLiM']['nucleotide_based'], True)
            self.assertEqual(in_ts.reference_sequence, out_ts.reference_sequence)
