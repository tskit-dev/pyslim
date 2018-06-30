# pySLiM

SLiM can now read, and write tree sequences, which store genealogical data of populations.
SLiM can use a tree sequence produced by the coalescent simulator `msprime` to initialize
a simulation, but to do so we need to add the relevant metadata. SLiM can also write out
history as a tree sequence, and in so doing it records extra information in metadata.
This package makes it easy to add the relevant metadata to a tree sequence so that SLiM
can use it, and to read the metadata in a SLiM-produced tree sequence.

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
new_ts = pyslim.annotate_defaults(ts)
new_ts.dump("slim_ts.trees")
```

## Quickstart: reading SLiM metadata

To retrieve the extra information that SLiM stores in a tree sequence, use the `extract_X_metadata()`
functions. For instance, to see the age of each Individual produced by annotation
in the previous example:
```
for ind in extract_mutation_metadata(new_ts.tables):
    print(ind.age)
```


## Important notes:

1. `tskit` "nodes" correspond to SLiM "genomes".  Individuals are diploid, so each individual has two nodes.

2. The currently alive individuals will be those in the Individual table; all other nodes will *not*
    have a corresponding individual.

3. The "remembered nodes" will be the *first* nodes.
