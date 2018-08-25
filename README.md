# pySLiM

SLiM can now read, and write tree sequences, which store genealogical data of populations.
SLiM can use a tree sequence produced by the coalescent simulator `msprime` to initialize
a simulation, but to do so we need to add the relevant metadata. SLiM can also write out
history as a tree sequence, and in so doing it records extra information in metadata.
This package makes it easy to add the relevant metadata to a tree sequence so that SLiM
can use it, and to read the metadata in a SLiM-produced tree sequence.

The SLiM manual documents how the extra metadata is stored in the tree sequence files,
and provides additional examples of how to use this package.

## Installation

To install `pyslim`, do
```
git clone https://github.com/tskit-dev/pyslim.git
cd pyslim
python setup.py install --user
```
You should also be able to install it with `pip install pyslim`.
You'll also need an up-to-date [msprime](https://github.com/tskit-dev/msprime) and [SLiM](https://messerlab.org/slim/).

To run the tests to make sure everything is working, do:
```
python -m nose tests
```

*Note:* if you use `python3` you may need to replace `python` with `python3` above.


## Quickstart: coalescent simulation for SLiM

Coalescent simulation is more biologically limited, but is much faster than
forwards simluation, so it can be helpful to run a "hybrid" simulation, by
endowing a SLiM simulation with a history derived from a msprime simulation.
The `pyslim.annotate()` command helps make this easy, by adding default
information to a tree sequence, allowing it to be read in by SLiM. This will
simulate a tree sequence with msprime, add SLiM information, and write it out
to a `.trees` file:
```
import msprime
import pyslim

# simulate a tree sequence of 12 nodes
ts = msprime.simulate(12, mutation_rate=1.0, recombination_rate=1.0)
new_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
new_ts.dump("initialize_nonWF.trees")
```

The resulting file `slim_ts.trees` can be read into SLiM to be used as a starting state,
as illustrated in this minimal example:
```
initialize()
{
    setSeed(23);
    initializeSLiMModelType("nonWF");
    initializeTreeSeq();
    initializeMutationRate(1e-2);
    initializeMutationType("m1", 0.5, "f", -0.1);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, 99);
    initializeRecombinationRate(1e-2);
}

1 early() { 
    sim.readFromPopulationFile("initialize_nonWF.trees");
}

10 {
    sim.treeSeqOutput("nonWF_restart.trees");
    catn("Done.");
    sim.simulationFinished();
}
```
See the SLiM manual for more about this operation.


## Quickstart: finishing an uncoalesced SLiM simulation

Although we can initialize a SLiM simulation with the results of a coalescent simulation,
if during the simulation we don't actually use the genotypes for anything, it can be much
more efficient to only coalesce the portions of the first-generation ancestors that have
not yet coalesced. (See the SLiM manual for more explanation, and the upcoming paper.)
Doing this is as simple as:
```
ts = pyslim.load("unfinished.trees")
ts.recapitate(recombination_rate = 1e-6, Ne=1000)
```

## Quickstart: adding neutral mutations to a SLiM simulation

If you have recorded a tree sequence in SLiM, likely you have not included any neutral mutations.
To add these (in a completely equivalent way to having included them during the simulation),
you can use the `msprime.mutate( )` function, which returns a new tree sequence with additional mutations.
This works as follows:
```
ts = pyslim.load("unmutated.trees")
mut_ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=1e-8, keep=True))
```
This will add infinite-sites mutations at a rate of 1e-8 per site, and will
`keep` any existing mutations (and not add any new ones to the sites where they
exist already). We have wrapped the call to `msprime.mutate` in a call to
`pyslim.SlimTreeSequence`, because `msprime.mutate` returns an *msprime* tree sequence,
and by converting it back into a `pyslim` tree sequence we can still use the methods
defined by `pyslim`. (The conversion does not modify the tree sequence at all,
it only adds the `.slim_generation` attribute.) The output of other `msprime`
functions that return tree sequences may be converted back to
`pyslim.SlimTreeSequence` in the same way.


## Quickstart: reading SLiM metadata

To retrieve the extra information that SLiM stores in a tree sequence, use the `extract_X_metadata()`
functions, where `X` is one of `mutation`, `population`, `node`, or `individual`.
For instance, to see the age of each Individual produced by annotation
in the previous example:
```
for ind in pyslim.extract_individual_metadata(new_ts.tables):
    print(ind.age)
```
In this example, all the ages are 0 (the default).

## Quickstart: modifying SLiM metadata

To *modify* the metadata that `pyslim` has introduced into a coalescent simulation,
or the metadata in a SLiM-produced tree sequence, use the `annotate_X_metadata()` functions.
For instance, to set the ages of the individuals in the tree sequence to random numbers between 1 and 4,
and write out the resulting tree sequence:
```
import random

tables = new_ts.dump_tables()
ind_md = list(pyslim.extract_individual_metadata(tables))
for ind in ind_md:
    ind.age = random.choice([1,2,3,4])

pyslim.annotate_individual_metadata(tables, ind_md)
mod_ts = pyslim.load_tables(tables, slim_format=True)

for ind in pyslim.extract_individual_metadata(mod_ts.tables):
    print(ind.age)

mod_ts.dump("modified_ts.trees")
```

# Documentation

Here we describe the technical details.
Currently, the python package `msprime` simulates tree sequences
and provides tools to work with them.
In the future, tools for working with tree sequences will be separated
into a package called `tskit`.
`pyslim` provides a thin interface between `msprime/tskit`.

## Metadata entries

SLiM records additional information in the metadata columns of Population, Individual, Node, and Mutation tables.
The information is recorded in a binary format, and is extracted and written by `pyslim` using the python `struct` module.
Nothing besides this binary information can be stored in the metadata of these tables for the tree sequence to be used by SLiM,
and so when `pyslim` annotates an existing tree sequence, anything in those columns is overwritten.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.

## Other important notes:

1. `tskit` "nodes" correspond to SLiM "genomes".  Individuals in SLiM are diploid, so each has two nodes.

2. Since in SLiM, all individual are diploid, every individual will be associated with two nodes.
    The Individual table contains entries for 

    a. the currently alive individuals, 
    b. any individuals that have been remembered with `treeSeqRememberIndividuals()`, and
    c. the *first* generation of the SLiM simulation.

    This last category is here because they are necessary for recapitation (described above);
    but they are *not* marked as samples, so will most likely be removed if you `simplify` the tree sequence.

3. The "remembered nodes" will be the *first* nodes.
