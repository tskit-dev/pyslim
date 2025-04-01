from collections import namedtuple
import os
import pytest
import random
import subprocess

from filelock import FileLock
import json
import tskit, msprime

from .recipe_specs import recipe_specs

class Outfiles:
    """
    Simple wrapper for a list of Outfile items: these specify a file path that SLiM can
    use to output a file, the SLiM variable name containing that path, for use in the
    recipe, a function to apply to the path to create an object for the unit test to
    access and a "name" (key) in the dictionary to access that object in the unit test.
    """
    Outfile = namedtuple("Outfile", "path, slim_name, post_process, key")

    @staticmethod
    def read_dictionary(fname):
        # reads in the json-written dictionary
        assert os.path.isfile(fname)
        with open(fname, "r") as f:
            out = json.load(f)
        return out

    @staticmethod
    def parse_pedigree_info(fname):
        # returns a dictionary whose keys are SLiM individual IDs, and whose values
        # are dictionaries with two entries:
        # - 'parents' is the SLiM IDs of the parents
        # - 'age' is a dictionary whose keys are tuples (SLiM generation, stage)
        #   and whose values are ages (keys not present are ones the indivdiual was
        #   not alive for)
        assert os.path.isfile(fname)
        out = {}
        with open(fname, 'r') as f:
            header = f.readline().split()
            assert header == ['generation', 'stage', 'individual', 'age', 'parent1', 'parent2']
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
                        assert p in out[ind]['parents']
                out[ind]['age'][(gen, stage)] = age
        return out

    def load_ts(self, path):
        out = {}
        if os.path.isfile(path):
            out['default'] = tskit.load(path)
        elif os.path.isdir(path):
            chroms = os.listdir(path)
            for cfile in os.listdir(path):
                c, e = os.path.splitext(cfile)
                if e == ".trees":
                    out[c] = tskit.load(os.path.join(path, cfile))
        return out

    def __init__(self, out_dir):
        # Note: the 'key' below cannot match a 'key' in recipe_specs
        self._outfiles = [
            self.Outfile(
                path=os.path.join(out_dir, "out.trees"),
                slim_name="TREES_FILE",  # The var containing the path name for SLiM output
                post_process=self.load_ts,  # function applied on path to create returned obj
                key="ts",  # The key to use for the post-processed item in the returned dict
            ),
            self.Outfile(
                path=os.path.join(out_dir, "out.pedigree"),
                slim_name="PEDIGREE_FILE",
                post_process=self.parse_pedigree_info,
                key="info",
            ),
            self.Outfile(
                path=os.path.join(out_dir, "out_mutations.json"),
                slim_name="MUTATIONS_FILE",
                post_process=self.read_dictionary,
                key="mutation_info",
            ),
        ]

    def __getitem__(self, index):
        return self._outfiles[index]

    def __len__(self, index):
        return len(self._outfiles)
    
    def results(self):
        res = {"path": {}}
        for o in self._outfiles:
            res["path"][o.key] = o.path
            if os.path.isfile(o.path) or os.path.isdir(o.path):
                assert o.key not in res
                res[o.key] = o.post_process(o.path) 
        return res


class MultiOutfiles(Outfiles):
    
    def __init__(self, out_dir, species):
        # Note: the 'key' below cannot match a 'key' in recipe_specs
        self._outfiles = []
        for sp in species:
            self._outfiles.extend([
                self.Outfile(
                    path=os.path.join(out_dir, f"{sp}_out.trees"),
                    slim_name=f"TREES_FILE_{sp}",  # The var containing the path name for SLiM output
                    post_process=tskit.load,  # function applied on path to create returned obj
                    key=f"ts_{sp}",  # The key to use for the post-processed item in the returned dict
                ),
                self.Outfile(
                    path=os.path.join(out_dir, f"{sp}_out.pedigree"),
                    slim_name=f"PEDIGREE_FILE_{sp}",
                    post_process=self.parse_pedigree_info,
                    key=f"info_{sp}",
                ),
                self.Outfile(
                    path=os.path.join(out_dir, f"{sp}_out_mutations.json"),
                    slim_name=f"MUTATIONS_FILE_{sp}",
                    post_process=self.read_dictionary,
                    key=f"mutation_info_{sp}",
                ),
            ])


def run_slim(recipe, out_dir, recipe_dir="test_recipes", species=None, **kwargs):
    """
    Run the recipe, present in the recipe_dir, outputting to files in out_dir
    kwargs are passed as variables to SLiM (in addition to the outfiles)
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    full_recipe = os.path.abspath(os.path.join(script_dir, recipe_dir, recipe))
    if not os.path.isfile(full_recipe):
        raise RuntimeError(f"{full_recipe} cannot be found")
    assert os.path.isdir(out_dir)  # should have been generated by the calling function
    if species is None:
        outfiles = Outfiles(out_dir)
    else:
        outfiles = MultiOutfiles(out_dir, species)
    slim_vars = []
    for o in outfiles:
        if o.slim_name != "":
            # for happy windows filepaths
            tmp_str = o.path.replace('\\', '\\\\')
            slim_vars += ["-d", f"{o.slim_name}=\"{tmp_str}\""]
    for k in kwargs:
        x = kwargs[k]
        if x is not None:
            if isinstance(x, str) and x[:10] != 'Dictionary':
                # for happy windows filepaths
                x = x.replace('\\', '\\\\')
                x = f"'{x}'"
            if isinstance(x, bool):
                x = 'T' if x else 'F'
            slim_vars += ["-d", f"{k}={x}"]
    cmd = ["slim", "-s", "23"] + slim_vars + [full_recipe]
    print(f"Running {' '.join(cmd)}, outputting errors etc. to '{out_dir}/SLiM_run_output.log'")
    with open(os.path.join(out_dir, "SLiM_run_output.log"), "w") as out:
        retval = subprocess.call(cmd, stderr=subprocess.STDOUT, stdout=out)
    if retval != 0:
        errmsg = f"Could not run {' '.join(cmd)}"
        errline = False
        with open(os.path.join(out_dir, "SLiM_run_output.log"), "r") as err:
            for line in err:
                if line[:5] == "ERROR" or line[:5] == "Error":
                    errline = True
                    errmsg += "\n"
                if errline and len(line.strip()) > 0:
                    errmsg += line
        raise RuntimeError(errmsg)
    return outfiles.results()


def _obj_to_eidos(x):
    if isinstance(x, str):
        out = '"' + x + '"'
    elif isinstance(x, dict):
        out = _dict_to_eidos(x)
    elif x is None:
        out = 'NULL'
    else:
        out = x
    return out


def _dict_to_eidos(d):
    txt = "Dictionary(" + ", ".join([f"\"{a}\", {_obj_to_eidos(b)}" for a, b in d.items()]) + ")"
    return txt


class HelperFunctions:
    @classmethod
    def run_slim_restart(cls, in_ts, recipe, out_dir, multichrom, subpop_map=None, **kwargs):
        # Saves out the tree sequence to a trees file and run the SLiM recipe
        # on it, saving to files in out_dir.
        assert isinstance(in_ts, dict)
        infile = os.path.join(out_dir, "in.trees")
        if not multichrom:
            assert len(in_ts) == 1
            ts = list(in_ts.values())[0]
            ts.dump(infile)
        else:
            os.mkdir(infile)
            for chrom, ts in in_ts.items():
                ts.dump(os.path.join(infile, f"{chrom}.trees"))
        # for happy windows filepaths       
        infile_str = infile.replace('\\', '\\\\')
        kwargs['TREES_IN'] = infile
        if 'STAGE' not in kwargs:
            kwargs['STAGE'] = ts.metadata['SLiM']['stage']
        if subpop_map is not None:
            # subpopMap argument has entries like { "p1" : 2 },
            # meaning that when loaded in, subpop p1 will be the second row
            # in the population table
            spm = _dict_to_eidos(subpop_map)
            kwargs['SUBPOP_MAP'] = spm
        results = run_slim(recipe, out_dir, **kwargs)
        return results["ts"]

    @classmethod
    def run_multispecies_slim_restart(cls, in_ts, recipe, out_dir, subpop_map=None, **kwargs):
        # in_ts should now be a dictionary with keys corresponding to species names
        infiles = {}
        for sp in in_ts:
            infiles[sp] = os.path.join(out_dir, f"{sp}_in.trees")
            in_ts[sp].dump(infiles[sp])
        infiles_str_dict = {}
        for k, v in infiles.items():
            # for happy windows filepaths
            infiles_str_dict[k] = v.replace('\\', '\\\\')        
        kwargs['TREES_IN'] = _dict_to_eidos(infiles_str_dict)
        if 'STAGE' not in kwargs:
            kwargs['STAGE'] = in_ts[sp].metadata['SLiM']['stage']
        if subpop_map is not None:
            # subpopMap argument has entries like { "p1" : 2 },
            # meaning that when loaded in, subpop p1 will be the second row
            # in the population table
            kwargs['SUBPOP_MAP'] = _dict_to_eidos(subpop_map)
        results = run_slim(recipe, out_dir, species=list(in_ts.keys()), **kwargs)
        out = { sp: results[f"ts_{sp}"] for sp in in_ts }
        return out

    @classmethod
    def run_msprime_restart(cls, in_ts, out_dir, multichrom=False, sex=None, WF=False):
        assert not multichrom
        L = int(list(in_ts.values())[0].sequence_length)
        out_ts = cls.run_slim_restart(
            in_ts,
            "restart_msprime.slim",  # This is standard script defined in test_recipes
            out_dir,
            multichrom=False,
            WF=WF,
            SEX=sex,
            L=L,
        )
        return out_ts

    @staticmethod
    def get_msprime_examples():
        # NOTE: we use DTWF below to avoid rounding of floating-point times
        # that occur with a continuous-time simulator
        seed = 6
        mutrate = 0.01
        for n in [2, 20]:
            for recrate in [0.0, 0.01]:
                ts = msprime.sim_ancestry(
                        n,
                        recombination_rate=recrate,
                        population_size=10,
                        sequence_length=200,
                        random_seed=seed,
                        model="dtwf"
                )
                yield ts
                mts = msprime.sim_mutations(
                        ts,
                        model=msprime.SLiMMutationModel(type=0),
                        rate=mutrate
                )
                yield mts
                seed += 1
                demography = msprime.Demography()
                demography.add_population(name="A", initial_size=20)
                demography.add_population(name="B", initial_size=30)
                demography.add_population(name="C", initial_size=10)
                demography.add_population_split(
                        time=5,
                        derived=["A", "B"],
                        ancestral="C"
                )
                ts = msprime.sim_ancestry(
                        {"A": n, "B": n},
                        demography=demography,
                        recombination_rate=recrate,
                        sequence_length=250,
                        random_seed=seed,
                        model="dtwf"
                )
                yield ts
                seed += 1


@pytest.fixture
def helper_functions():
    return HelperFunctions


@pytest.fixture(scope="session", params=recipe_specs.keys())
def recipe(request, tmp_path_factory, worker_id):
    """
    Most logic in this fixture is to check whether we are running as a single proc, or as a worker.
    If a worker, and this is the first run, we need to lock the process until the simulation has finished.
    """
    recipe_name = request.param
    if worker_id == "master":
        out_dir = tmp_path_factory.getbasetemp() / recipe_name
    else:
        # get the temp directory shared by all workers, so that simulation files
        # are shared between workers
        out_dir = tmp_path_factory.getbasetemp().parent / recipe_name
    with FileLock(str(out_dir) + ".lock"):
        if out_dir.is_dir():
            ret = Outfiles(out_dir).results()
        else:
            os.makedirs(out_dir, exist_ok=True)
            ret = run_slim(recipe=recipe_name, out_dir=out_dir)
    
    ret.update(recipe_specs[recipe_name])
    return ret
