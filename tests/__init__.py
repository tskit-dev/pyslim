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

# recipes that record everyone ever
_everyone_example_files = ["tests/examples/recipe_record_everyone"]

_example_files = ["tests/examples/recipe_{}".format(x)
                  for x in ['WF', 'nonWF', 'nucleotides']] + \
                 _everyone_example_files

# this is of the form (input, basename)
# TODO: test restarting of nucleotides after reference sequence dumping is enabled
_restart_files = [("tests/examples/recipe_{}.trees".format(x),
                   "tests/examples/restart_{}".format(x))
                  for x in ['WF', 'nonWF']] # , 'nucleotides']]

def run_slim_script(slimfile, args=''):
    outdir = os.path.dirname(slimfile)
    script = os.path.basename(slimfile)
    print("running " + "cd " + outdir + " && slim -s 23 " + args + " " + script)
    out = os.system("cd " + outdir + " && slim -s 23 " + args + " " + script + ">/dev/null")
    return out

def setUp():
    # Make random tests reproducible.
    random.seed(210)

    # run SLiM
    for basename in _example_files:
        treefile = basename + ".trees"
        print(treefile)
        try:
            os.remove(treefile)
        except FileNotFoundError:
            pass
        slimfile = basename + ".slim"
        out = run_slim_script(slimfile)
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
        for treefile in self.get_slim_example_files():
            print("---->", treefile)
            self.assertTrue(os.path.isfile(treefile))
            yield pyslim.load(treefile)

    def get_slim_everyone_examples(self):
        for filename in _everyone_example_files:
            treefile = filename + ".trees"
            print("---->", treefile)
            self.assertTrue(os.path.isfile(treefile))
            yield pyslim.load(treefile)

    def get_slim_restarts(self):
        for treefile, basename in _restart_files:
            self.assertTrue(os.path.isfile(treefile))
            ts = pyslim.load(treefile)
            yield ts, basename

    def run_slim_restart(self, in_ts, basename, args=''):
        infile = basename + ".init.trees"
        outfile = basename + ".trees"
        slimfile = basename + ".slim"
        for treefile in infile, outfile:
            try:
                os.remove(treefile)
            except FileNotFoundError:
                pass
        in_ts.dump(infile)
        out = run_slim_script(slimfile, args=args)
        print("out:", out)
        try:
            os.remove(infile)
        except FileNotFoundError:
            pass
        assert out == 0
        self.assertTrue(os.path.isfile(outfile))
        out_ts = pyslim.load(outfile)
        try:
            os.remove(outfile)
        except FileNotFoundError:
            pass
        return out_ts

    def run_msprime_restart(self, in_ts, sex=None):
        basename = "tests/examples/restart_msprime"
        args = " -d L={}".format(int(in_ts.sequence_length))
        if sex is not None:
            args = args + " -d \"SEX='{}'\"".format(sex)
        out_ts = self.run_slim_restart(in_ts, basename, args=args)
        return out_ts

