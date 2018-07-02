# pySLiM

SLiM can now read, and write tree sequences, which store genealogical data of populations.
SLiM can use a tree sequence produced by the coalescent simulator `msprime` to initialize
a simulation, but to do so we need to add the relevant metadata. SLiM can also write out
history as a tree sequence, and in so doing it records extra information in metadata.
This package makes it easy to add the relevant metadata to a tree sequence so that SLiM
can use it, and to read the metadata in a SLiM-produced tree sequence.

The SLiM manual documents how the extra metadata is stored in the tree sequence files,
and provides examples of how to use this package.

## Installation

To install `pyslim`, do
```
git clone https://github.com/tskit-dev/pyslim.git
cd pyslim
python setup.py install --user
```
You'll also need an up-to-date [msprime](https://github.com/tskit-dev/msprime) and [SLiM](https://messerlab.org/slim/), of course.

To run the tests to make sure everything is working, do:
```
cd tests/examples
for x in *.slim; do slim $x; done
cd -
python -m nose tests
```

*Note:* if you use `python3` you may need to replace `python` with `python3` above.


## Quickstart: coalescent simulation for SLiM

The `pyslim.annotate()` command will add default information to a tree sequence, allowing
it to be read in by SLiM. This will simulate a tree sequence with msprime, add SLiM information,
and write it out to a `.trees` file:
```
import msprime
import pyslim

# simulate a tree sequence of 12 nodes
ts = msprime.simulate(12, mutation_rate=1.0, recombination_rate=1.0)
new_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
new_ts.dump("slim_ts.trees")
```

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

## Metadata entries

SLiM records additional information in the metadata columns of Population, Individual, Node, and Mutation tables.
The information is recorded in a binary format, and is extracted and written by `pyslim` using the python `struct` module.
This binary information can be the *only* thing in the metadata of these tables for the tree sequence to be used by SLiM,
and so when `pyslim` annotates an existing tree sequence, anything in those columns is overwritten.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.

## Time and SLiM Tree Sequences

The "time" in SLiM records the number of generations *since* the beginning of the simulation.
However, since tree sequences naturally deal with history in retrospectively,
properties related to "time" of tree sequences are measured in generations *ago*,
i.e., generation *prior* to a given point.
To distinguish these two notions of time, we'll talk about "SLiM time" and "tskit time".
When SLiM records a tree sequence,
it records tskit time in units of generations before the start of the simulation;
so that all tskit times it records in the tree sequence are *negative*
(because it measures how long *before* the start, but happened *after* the start).
This is terribly counterintuitive. The fact that there are two notions of time 
- one moving forwards, the other backwards -
is unavoidable, but `pyslim` does one thing to make this easier to work with:
when `pyslim` loads in a tree sequence file, it checks to see what the current SLiM time was
at the end of the simulation, and shifts all times in the tree sequence so that tskit time
is measured in generations before the *end* of the simulation.

The upshot is that:

1. The ``time`` attribute of tree sequence nodes gives the number of generations
    before then end of the simulation at which those nodes were born.

2. These numbers will *not* match the values in the `.trees` file, but you should not
    need to worry about that, as long as you always `load()` and `dump()` using `pyslim`.

3. The conversion factor, the "SLiM time" when the tree sequence file was written,
    is stored in the `slim_generation` attribute of a `SlimTreeSequence`.

An example should help clarify things.
Suppose that `my.trees` is a file that was saved by SLiM at the end of a simulation run
for 100 generations, and that
we want to find the list of nodes in the tree sequence that were born during the first
20 generations of the simulation.
Since the birth time of a node is recorded in the `.time` attribute of a node
in tskit time, a node that was born in the first 20 generations of the simulation,
i.e., more than 80 generations before the end of the simulation,
will have a `.time` attribute of at least 80.
```
ts = pyslim.load("my.trees", slim_format=True)
old_nodes = []
for n in ts.nodes():
   if n.time > ts.slim_generation - 20:
        old_nodes.append(n)
```

In the future, we may change this behavior,
but if so will provide an upgrade path for old files.

## SLiM Tree Sequences

Because SLiM adds additional information to tree sequences,
`pyslim` defines a subclass of `msprime.TreeSequence` to make it easy to access this information,
and to make the time shift described above seamless.
When you run `pyslim.load('my.trees', slim_format=True)`,
you get a `SlimTreeSequence` object.
This has all the same properties and methods as a plain `TreeSequence`,
with the following differences:

1. It has a `slim_generation` attribute.
2. It's `.dump()` method shifts times by `slim_generation` before writing them out,
    so that `pyslim.load("my.trees", slim_format=True).dump("my2.trees")`
    writes out a tree sequence identical to the one that was read in (up to floating point error).


## Other important notes:

1. `tskit` "nodes" correspond to SLiM "genomes".  Individuals are diploid, so each individual has two nodes.

2. The currently alive individuals will be those in the Individual table; all other nodes will *not*
    have a corresponding individual.

3. The "remembered nodes" will be the *first* nodes.
