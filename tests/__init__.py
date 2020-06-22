"""
Common code for the pyslim test cases.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import pyslim
import msprime
import random
import unittest
import base64
import os
import attr

# possible attributes to simulation scripts are
#  WF, nonWF
#  nucleotides
#  everyone: records everyone ever
#  pedigree: writes out accompanying info file containing the pedigree
#  remembered_early: remembering and saving the ts happens during early
# All files are of the form `tests/examples/{key}.slim`
example_files = {}
example_files['recipe_nonWF'] = {"nonWF": True, "pedigree": True}
example_files['recipe_WF'] = {"WF": True, "pedigree": True}
example_files['recipe_nucleotides'] = {"WF": True, "pedigree": True, "nucleotides": True}
example_files['recipe_long_nucleotides'] = {"WF": True, "nucleotides": True}
example_files['recipe_roots'] = {"WF": True}
for t in ("WF", "nonWF"):
    for s in ("early", "late"):
        value = {t: True, "everyone": True, "pedigree": True}
        if s == 'early':
            value['remembered_early'] = True
        example_files[f'recipe_record_everyone_{t}_{s}'] = value


for f in example_files:
    print(f, example_files[f])
    example_files[f]['basename'] = os.path.join("tests", "examples", f)


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
    for f in example_files:
        basename = example_files[f]['basename']
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
    for f in example_files:
        basename = example_files[f]['basename']
        treefile = basename + ".trees"
        try:
            os.remove(treefile)
            pass
        except FileNotFoundError:
            pass
        infofile = treefile + ".pedigree"
        try:
            os.remove(infofile)
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

    def get_slim_examples(self, *args, return_info=False):
        for ex in example_files.values():
            basename = ex['basename']
            use = True
            for a in args:
                if a not in ex:
                    use = False
            if use:
                treefile = basename + ".trees"
                print("---->", treefile)
                self.assertTrue(os.path.isfile(treefile))
                ts = pyslim.load(treefile)
                if return_info:
                    infofile = treefile + ".pedigree"
                    if os.path.isfile(infofile):
                        ex['info'] = self.get_slim_info(infofile)
                    else:
                        ex['info'] = None
                    yield (ts, ex)
                else:
                    yield ts

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

    def get_msprime_examples(self):
        demographic_events = [
            msprime.MassMigration(
            time=5, source=1, destination=0, proportion=1.0)
        ]
        for n in [2, 10, 20]:
            for mutrate in [0.0]:
                for recrate in [0.0, 0.01]:
                    yield msprime.simulate(n, mutation_rate=mutrate,
                                           recombination_rate=recrate,
                                           length=200)
                    population_configurations =[
                        msprime.PopulationConfiguration(
                        sample_size=n, initial_size=100),
                        msprime.PopulationConfiguration(
                        sample_size=n, initial_size=100)
                    ]
                    yield msprime.simulate(
                        population_configurations=population_configurations,
                        demographic_events=demographic_events,
                        recombination_rate=recrate,
                        mutation_rate=mutrate,
                        length=250)

    def get_slim_info(self, fname):
        # returns a dictionary whose keys are SLiM individual IDs, and whose values
        # are dictionaries with two entries:
        # - 'parents' is the SLiM IDs of the parents
        # - 'age' is a dictionary whose keys are tuples (SLiM generation, stage)
        #   and whose values are ages (keys not present are ones the indivdiual was
        #   not alive for)
        self.assertTrue(os.path.isfile(fname))
        out = {}
        with open(fname, 'r') as f:
            header = f.readline().split()
            self.assertSequenceEqual(
                    header,
                    ['generation', 'stage', 'individual', 'age', 'parent1', 'parent2'])
            for line in f:
                gen, stage, ind, age, p1, p2 = line.split()
                gen = int(gen)
                ind = int(ind)
                age = int(age)
                parents = tuple([int(p) for p in (p1, p2) if p != "-1"])
                if ind not in out:
                    out[ind] = {
                            "parents" : parents,
                            "age" : {}
                            }
                else:
                    for p in parents:
                        self.assertIn(p, out[ind]['parents'])
                out[ind]['age'][(gen, stage)] = age
        return out
