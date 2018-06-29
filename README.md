# pySLiM

SLiM can now read, and write tree sequences, which store genealogical data of populations.
SLiM can use a tree sequence produced by the coalescent simulator `msprime` to initialize
a simulation, but to do so we need to add the relevant metadata. SLiM can also write out
history as a tree sequence, and in so doing it records extra information in metadata.
This package makes it easy to add the relevant metadata to a tree sequence so that SLiM
can use it, and to read the metadata in a SLiM-produced tree sequence.


## Usage

The `pyslim.annotate()` command will add default information to a tree sequence, allowing
it to be read in by SLiM. To further modify this information, you must extract the underlying
tables and modify them, like so:
```
import msprime
import pyslim

# simulate a tree sequence of 12 nodes
ts = msprime.simulate(12, mutation_rate=1.0, recombination_rate=1.0)
tables = pyslim.annotate(ts.tables)

# assign 3 females, 3 males
orig_individuals = tables.individuals.copy()
tables.individuals.clear()
for sex, ind in zip([1, 1, 1, 2, 2, 2], orig_individuals):
    tables.individuals.add_row(flags=sex, location=ind.location, metadata=ind.metadata)

new_slim_ts = pyslim.load_tables(tables)
```


## Important notes:

1. `tskit` "nodes" correspond to SLiM "genomes".  Individuals are diploid, so each individual has two nodes.

2. The currently alive individuals will be those in the Individual table; all other nodes will *not*
    have a corresponding individual.

3. The "remembered nodes" will be the *first* nodes.


## TODO

3. Save population metadata as ASCII somewhere.

5. Properly set (comma-separated) strings of mutation IDs in derived states, and set ancestral states to be ''.

6. 'time' in pyslim should be "time ago"; but all arguments should say "time_since_start" or "time_before_end".
    This is relevant for the "reference_time" argument and the "time" information in mutation metadata.

7. Check that remembered individuals/samples are dealt with correctly.

8. Assign mutation IDs in annotate().
