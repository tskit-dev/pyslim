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


