---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell}
:tags: [remove-cell]

import pyslim, tskit, msprime
from IPython.display import SVG
import numpy as np
import util

np.random.seed(1234)
```


(sec_vignette_continuing)=


# Vignette: Following up with more coalescent simulation

Previously, we saw how to use recapitation
to simulate the period *before* a SLiM simulation
with the coalescent simulator, msprime.
We can do the same thing *after* a period of SLiM simulation.
To demonstrate this,
below we'll run a simulation in which

1. A population evolves neutrally for a long time, but then
2. it experiences strong positive selection on new mutations for 100 generations, and
3. evolves neutrally for another 1000 generations.

To do this, we'll simulate step (2) first, with SLiM,
then recapitate to add step (1), and then "continue" the simulation using msprime to add in (3).


## Positive selection

Here's a SLiM script that has rapid, strong selection acting genome-wide for 20 generations.
It is perhaps not very realistic, but it's dramatic.

```{literalinclude} rapid_adaptation.slim
```
```{code-cell}
%%bash
slim -s 5 rapid_adaptation.slim
```

We can see what happened in the GUI,
but let's pull some more statistics out of the tree sequence:
```{code-cell}
ts = pyslim.load("rapid_adaptation.trees")

# allele frequencies
p = ts.sample_count_stat(
                [ts.samples()], lambda x: x/20000, 1, windows='sites',
                span_normalise=False, polarised=True, strict=False)
print(f"There are {ts.num_sites} segregating sites, of which {np.sum(p > 0.25)}")
print(f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")
```

The selection was, indeed, strong.

## Recapitation

Ok, now let's do phase (1), recapitating and mutating the result.
We'll add SLiM mutations with "mutation type" 0
(so in SLiM these would be called type `m0`),
so first we check that all the existing mutations are of a different type.

```{code-cell}
rts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-8, random_seed=6)

# check type m0 is not used:
mut_types = set([md['mutation_type']
                for mut in ts.mutations()
                for md in mut.metadata['mutation_list']])
print(f"Keeping {rts.num_mutations} existing mutations of type(s) {mut_types}.")
assert 0 not in mut_types

# add type m0 mutations
rts = pyslim.SlimTreeSequence(
        msprime.sim_mutations(
            rts, rate=1e-8, random_seed=7, keep=True,
            model=msprime.SLiMMutationModel(type=0))
        )

p = rts.sample_count_stat(
                [rts.samples()], lambda x: x/20000, 1, windows='sites',
                span_normalise=False, polarised=True, strict=False)
print(f"After mutation, there are {rts.num_sites} segregating sites, of which {np.sum(p > 0.25)}")
print(f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")
```

Now, there are more segregating sites - neutral ones.


## Continuing the simulation

To "continue" the simulation neutrally, we'll 

1. simulate the desired period of time in msprime
2. randomly match the initial ancestors in the msprime simulation
   with the final individuals of the SLim simulation, and
3. merge the two together, using the {meth}`tskit.TreeSequence.union` method.


**(1)** Simulating for a given period of time in msprime requires the ``end_time`` argument
(remembering that this is *time ago*);
we'll do this to simulate an additional 1000 generations.

This is almost what we need, but there is one more detail:
if complete coalescence occurs on any region of the genome,
msprime will stop simulating the history of that region.
This is a problem, since we need all lineages to extend back to ``end_time``.
To make sure all lineages trace back to ``end_time``,
we'll add one "fake" sample from a separate population, that *can't* coalesce with the rest,
then remove it before the next step, using the ``keep_input_roots=True`` argument to ``simplify()``.


```{code-cell}
new_time = 1000
demog_model = msprime.Demography()
demog_model.add_population(initial_size=10000, name='real')
demog_model.add_population(initial_size=10000, name='fake')
new_ts = msprime.sim_ancestry(
              samples={'real' : 10000, 'fake' : 1},
              demography=demog_model,
              end_time=new_time,
              sequence_length=rts.sequence_length,
              recombination_rate=1e-8,
              random_seed=9)
new_ts = msprime.sim_mutations(
                 new_ts, rate=1e-8, random_seed=10, keep=True,
                 model=msprime.SLiMMutationModel(type=0)
        )
new_tables = new_ts.tables
# check that the spurious samples are 20000 and 20001
for n in (20000, 20001):
   assert n in new_ts.samples()
   assert new_ts.node(n).population == 1
new_tables.simplify(samples=np.arange(20000), keep_input_roots=True)
print(f"Remaining number of populations: {new_tables.populations.num_rows}")
```

**(2)** Now we'll pull out the IDs of the nodes from 1000 generations ago,
shift the times in the SLiM tree sequence back 1000 generations,
randomly assign each to a node at the end of the SLiM simulation,
and merge them.

```{code-cell}

new_nodes = np.where(new_tables.nodes.time == new_time)[0]
print(f"There are {len(new_nodes)} nodes from the start of the new simulation.")
# There are 4425 nodes from the start of the new simulation.

slim_nodes = rts.samples(time=0)
assert(len(slim_nodes) == 20000)

# randomly give new_nodes IDs in rts
node_map = np.repeat(tskit.NULL, new_tables.nodes.num_rows)
node_map[new_nodes] = np.random.choice(slim_nodes, len(new_nodes), replace=False)

# shift times: in nodes and mutations
# since tree sequences are not mutable, we do this in the tables directly
# also, unmark the nodes at the end of the SLiM simulation as samples
tables = rts.tables
tables.nodes.flags = tables.nodes.flags & ~np.uint32(tskit.NODE_IS_SAMPLE)
tables.nodes.time = tables.nodes.time + new_time
tables.mutations.time = tables.mutations.time + new_time

# merge the two sets of tables
tables.union(new_tables, node_map,
            add_populations=False,
            check_shared_equality=False)

# get back the tree sequence
full_ts = pyslim.SlimTreeSequence(tables.tree_sequence())

p = full_ts.sample_count_stat(
                [full_ts.samples()], lambda x: x/20000, 1,
                windows='sites', span_normalise=False,
                polarised=True, strict=False)
print(f"There are {full_ts.num_sites} segregating sites, of which {np.sum(p > 0.25)}")
print(f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")
```

Well, allele frequencies have drifted.
Don't worry, we'll explain what happened there in a minute.

Let's do a consistency check.
First, here's the root of the first tree in the recapitated SLiM simulation:
```{code-cell}
:tags: ["remove-output"]
t = rts.first()
assert(t.num_roots == 1)
r = rts.node(t.root)
print(r)
```
```{code-cell}
:tags: ["remove-input"]
util.pp(r)
```
Now, here's the root of the first tree *after* continuing
for 1000 generations, which should be the same:
```{code-cell}
:tags: ["remove-output"]
ft = full_ts.first()
assert(ft.num_roots == 1)
fr = full_ts.node(ft.root)
print(fr)
```
```{code-cell}
:tags: ["remove-input"]
util.pp(fr)
```
That matches up - the time of what should be the same node in the "continued" tree sequence
is 1000 generations earlier.

So, what happened with ``union`` back there?
Well, the basic usage is ``tables.union(other, node_map)``,
where ``node_map`` is an array of length equal to the number of nodes in ``other``,
whose entries are either ``tskit.NULL`` or the ID of a node in ``tables``.
The entries that *aren't* NULL indicate that
``union`` should glue together ``tables`` and ``other`` by saying that that pair of nodes are the same.
(So, e.g., if ``node_map[3]`` is equal to ``25``, then it says that node 25 in ``tables``
is actually the same, really, as node 3 in ``other``.)
We then asked ``union`` to please not create new populations,
since otherwise it would have assigned all the new nodes to a new population.
We also asked it to not "check for overlap equality":
sometimes, when unioning together two tree sequences,
we really expect everything having to do with the set of nodes we're saying are identical
to be identical in the two tree sequences, so ``union`` by default throws an error if it's not.
We don't expect that in this case, because, for instance,
there could be a mutation above one of the terminal nodes in the SLiM tree sequence;
this would clearly not be present in the new tree sequence.
