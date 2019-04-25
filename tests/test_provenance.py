"""
Test cases for provenance handling.
"""
from __future__ import print_function
from __future__ import division

import pyslim
import msprime
import tests
import unittest
import random
import json
import tskit


class TestProvenance(tests.PyslimTestCase):

    def get_0_1_slim_examples(self):
        for filename in ['tests/examples/recipe_WF.v3.0.trees',
                         'tests/examples/recipe_nonWF.v3.0.trees']:
            with self.assertWarns(Warning):
                yield tskit.load(filename)

    def get_0_2_slim_examples(self):
        for filename in ['tests/examples/recipe_WF.v3.2.trees',
                         'tests/examples/recipe_nonWF.v3.2.trees']:
            with self.assertWarns(Warning):
                yield tskit.load(filename)

    def test_provenance_creation(self):
        record = pyslim.make_pyslim_provenance_dict()
        tskit.provenance.validate_provenance(record)

        record = pyslim.make_slim_provenance_dict("nonWF", 100)
        tskit.provenance.validate_provenance(record)

    def test_upgrade_provenance(self):
        ts = msprime.simulate(10)
        # test bad input
        with self.assertRaises(ValueError):
            pyslim.upgrade_slim_provenance(ts.tables)
        # test good input
        record = {"program": "SLiM", "version": "3.0", "file_version": "0.1",
                  "model_type": "WF", "generation": 10}
        prov = msprime.Provenance(timestamp='2018-08-25T14:59:13', record=json.dumps(record))
        is_slim, version = pyslim.provenance._slim_provenance_version(json.loads(prov.record))
        self.assertTrue(is_slim)
        self.assertEqual(version, "0.1")
        tables = ts.dump_tables()
        tables.provenances.add_row(json.dumps(record))
        pyslim.upgrade_slim_provenance(tables) # modifies the tables
        new_ts = tables.tree_sequence()
        self.assertEqual(new_ts.num_provenances, 3)
        new_record = json.loads(new_ts.provenance(2).record)
        is_slim, version = pyslim.provenance._slim_provenance_version(new_record)
        self.assertTrue(is_slim)
        self.assertEqual(version, "0.3")
        self.assertEqual(record['model_type'], new_record['parameters']['model_type'])
        self.assertEqual(record['generation'], new_record['slim']["generation"])

    def test_convert_0_1_files(self):
        for ts in self.get_0_1_slim_examples():
            pts = pyslim.SlimTreeSequence(ts)
            self.assertEqual(ts.num_provenances, 1)
            self.assertEqual(pts.num_provenances, 2)
            self.assertEqual(ts.provenance(0).record, pts.provenance(0).record)
            record = json.loads(ts.provenance(0).record)
            new_record = json.loads(pts.provenance(1).record)
            self.assertEqual(record['model_type'], new_record['parameters']['model_type'])
            self.assertEqual(record['generation'], new_record['slim']["generation"])
            self.assertListEqual(list(ts.samples()), list(pts.samples()))
            self.assertArrayEqual(ts.tables.nodes.flags, pts.tables.nodes.flags)
            samples = list(ts.samples())
            t = ts.first()
            pt = pts.first()
            for _ in range(20):
                u = random.sample(samples, 1)[0]
                self.assertEqual(t.parent(u), pt.parent(u))
                if t.parent(u) != msprime.NULL_NODE:
                    self.assertEqual(t.branch_length(u), pt.branch_length(u))

    def test_convert_0_2_files(self):
        for ts in self.get_0_2_slim_examples():
            pts = pyslim.SlimTreeSequence(ts)
            self.assertEqual(ts.num_provenances, 1)
            self.assertEqual(pts.num_provenances, 2)
            self.assertEqual(ts.provenance(0).record, pts.provenance(0).record)
            record = json.loads(ts.provenance(0).record)
            new_record = json.loads(pts.provenance(1).record)
            self.assertEqual(record['parameters']['model_type'],
                             new_record['parameters']['model_type'])
            self.assertEqual(record['slim']['generation'],
                             new_record['slim']["generation"])
            self.assertListEqual(list(ts.samples()), list(pts.samples()))
            self.assertArrayEqual(ts.tables.nodes.flags, pts.tables.nodes.flags)
            samples = list(ts.samples())
            t = ts.first()
            pt = pts.first()
            for _ in range(20):
                u = random.sample(samples, 1)[0]
                self.assertEqual(t.parent(u), pt.parent(u))
                if t.parent(u) != msprime.NULL_NODE:
                    self.assertEqual(t.branch_length(u), pt.branch_length(u))


