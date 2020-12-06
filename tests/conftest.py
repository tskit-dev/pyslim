"""
To make pytest run this setup/teardown code it needs to go here, I guess.
"""

import os
import random
import pytest


from . import example_files, run_slim_script


@pytest.fixture(scope="session", autouse=True)
def setup_slim_examples():
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

    # the yield makes pytest wait and run tests before doing the rest of this
    yield

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
