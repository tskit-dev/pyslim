"""
Common code for the pyslim test cases.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import pyslim
import random
import unittest
import base64
import os

_example_files = ["tests/examples/recipe_WF",
                 "tests/examples/recipe_nonWF"]

def setUp():
    # Make random tests reproducible.
    random.seed(210)

    # run SLiM
    for filename in _example_files:
        treefile = filename + ".trees"
        print(treefile)
        try:
            os.remove(treefile)
        except FileNotFoundError:
            pass
        outdir = os.path.dirname(filename)
        slimfile = os.path.basename(filename) + ".slim"
        print("running " + "cd " + outdir + " && slim -s 23 " + slimfile)
        out = os.system("cd " + outdir + " && slim -s 23 " + slimfile + ">/dev/null")
        assert out == 0


def tearDown():
    for filename in _example_files:
        treefile = filename + ".trees"
        try:
            os.remove(treefile)
            pass
        except FileNotFoundError:
            pass


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

    def verify_haplotype_equality(self, ts, slim_ts):
        self.assertEqual(ts.num_sites, slim_ts.num_sites)
        for j, v1, v2 in zip(range(ts.num_sites), ts.variants(),
                             slim_ts.variants()):
            g1 = [v1.alleles[x] for x in v1.genotypes]
            g2 = [v2.alleles[x] for x in v2.genotypes]
            self.assertArrayEqual(g1, g2)

    def get_slim_example_files(self):
        for filename in _example_files:
            yield filename + ".trees"

    def get_slim_examples(self):
        for filename in self.get_slim_example_files():
            print("---->", filename)
            yield pyslim.load(filename)

