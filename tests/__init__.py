"""
Common code for the pyslim test cases.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import random
import unittest
import base64


def setUp():
    # Make random tests reproducible.
    random.seed(210)


class PyslimTestCase(unittest.TestCase):
    '''
    Base class for test cases in pyslim.
    '''

    def assertArrayEqual(self, x, y):
        self.assertListEqual(list(x), list(y))

    def assertArrayAlmostEqual(self, x, y):
        self.assertEqual(len(x), len(y))
        for a, b in zip(x, y):
            self.assertAlmostEqual(a, b)

    def assertTablesAlmostEqual(self, t1, t2):
        # nodes
        self.assertEqual(t1.nodes.num_rows, t2.nodes.num_rows)
        self.assertArrayEqual(t1.nodes.flags, t2.nodes.flags)
        self.assertArrayEqual(t1.nodes.individual, t2.nodes.individual)
        self.assertArrayEqual(t1.nodes.metadata, t2.nodes.metadata)
        self.assertArrayEqual(t1.nodes.metadata_offset, t2.nodes.metadata_offset)
        self.assertArrayEqual(t1.nodes.population, t2.nodes.population)
        self.assertArrayAlmostEqual(t1.nodes.time, t2.nodes.time)
        # other tables
        self.assertEqual(t1.edges, t2.edges)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)
        self.assertEqual(t1.migrations, t2.migrations)
        self.assertEqual(t1.individuals, t2.individuals)
        self.assertEqual(t1.provenances, t2.provenances)
