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
import pandas as pd 
import os

np.random.seed(1234)
```


(sec_vignette_parallel)=


# Vignette: Parallelizing SLiM simulations in a phylogenetic tree

Imagine you want to simulate the evolutionary history of a group. If there is
no migration between any of the branches in your tree, any branches stemming
from the same node can be simulated in parallel (see {numref}`phylo`).

```{figure} _static/phylo.png
:height: 200px
:name: phylo

Example of phylogeny we might want to simulate. Note how branches with the same color can be simulated in parallel when there is no migration.
```

To do this, we'll need to do two things:
(1) be able to *simulate* branches in parallel, and
(2) glue the resulting simulations (one per tip) back together.


## Simulating the branches

First, we need to write a SLiM script that will be used for simulating the
history of each branch in our phylogeny.
We will perform a simple simulation, in which each branch can have a different
(but fixed) population size and length (number of generations).
Also, we will allow deleterious mutations to happen across the entire chromosome
at a fixed rate.

Here is a SLiM script that would do this:

```{literalinclude} phylo_bgs.slim
```

For each branch, the presence or absence of ``infile`` tells SLiM
whether we want to start it from a previous branch or not.
If so, SLiM will read the previous tree sequence and change the
population size accordingly.
Note that when you read a tree sequence into SLiM, the generation counter will
be updated with the time encoded in the tree sequence, so we need to set the end
of the simulation as the length of the branch (`num_gens`) plus the current
"time" at the end of the loaded tree sequence.
At the end of the simulation, we call `sim.treeSeqRememberIndividuals` right
before saving the resulting tree sequence. This is necessary because we need to
ensure the individuals in the final generation are never dropped from the tree
sequence in future runs of SLiM which are started from the output of the
simulation, as they will later be used to glue the tree sequences together.

I encoded the phylogeny we will simulate in a simple table,
which we'll use as ``df`` in the code below:

```{code-cell}
:tags: ["hide-input"]
df = pd.read_csv("_static/phylo.tsv", sep="\t")
df = df.fillna('')
df["infile"] = df.parent + ".trees"
df["outfile"] = df.child + ".trees"
df.loc[df["infile"]==".trees", "infile"] = ""
df["is_leaf"] = ~df.child.isin(df.parent)
df
```

With our phylogeny and the simulation parameters, we are ready to run our
simulations.
One way to parallelize the simulation of sister branches is to use `make`.
You do not need to know much about this tool (though it is totally worthwhile
to check it out).
The main idea here is that you can specify dependency between files and `make`
works its magic to run the simulations in the right order.
Here is python code that will write out a makefile from the information in ``df``:

```{code-cell}
f = open("sims.make", "w")
print(f"all: {' '.join(df.outfile.to_list())}\n", file=f)
for i, row in df.iterrows():
    print(f"{row.outfile}: {row.infile}", file=f)
    print(f"\tslim -d \"infile='{row.infile}'\" -d popsize={row.popsize} -d num_gens={row.edgelen} -d \"outfile='{row.child}.trees'\" phylo_bgs.slim\n", file=f)
f.close()
```

Here's the result. Again, don't worry about the details,
but you can see that the file encodes the phylogeny
through a bunch of ``child : parent`` "rules":
```{code-cell}
:tags: ["hide-input"]
%%bash
cat sims.make
```

With the makefile in hand,
we can now run make, specifying the maximum number of simulations
to be run simultaneously the ``-j``.

```python
%%bash
make -f sims.make -j 3
```

```{dropdown} Alternative to using make

You would have to write a recursion over the branches in your tree (starting 
from the root) and then parallelize the runs of sister branches somehow.

```python
def phylo_recursion(parent, df):
    print(parent)
    childs = df[df.parent==parent]
    print(childs)
    if len(childs) == 0:
        return
    # you could parallelize this loop over childs with same parent
    for i, row in childs.iterrows():
        if not os.path.exists(row.outfile):
            os.system(f"slim -d \"infile='{row.infile}'\" -d popsize={row.popsize} -d num_gens={row.edgelen} -d \"outfile='{row.child}.trees'\" phylo_bgs.slim")
            phylo_recursion(row.child, df)

phylo_recursion("", df)
```

## Putting it all together: unioning the tree sequences

With the tree sequences in hand, we now need to glue them together.
This can be done using
[**union**](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.union)
from tskit.
For two tree sequences which share some of its past history is shared, **union**
works by copying the non-shared parts of one of the tree sequence onto the other.
The trickiest part of this operation is defining the parts that are equivalent
in the two tree sequences. For that, you will have to create an array that serves
as a map of node IDs between the two tree sequences.

Here is a function that will construct a map of the node IDs of two SLiM tree sequences 
that correspond to the same chromosomes in SLiM
at any time older than the given time ago at which the two populations split.
Given two tree sequences ``other`` and ``ts``,
the goal here is to find,
for each node born before ``split_time`` ago in ``other``,
the matching node in ``ts``, where we can identify matching using the SLiM ID in metadata.
The code could be made easier to read by iterating over nodes,
but the following numpy-based version is much faster:

```{code-cell}
def match_nodes(other, ts, split_time):
    """
    Given SLiM tree sequences `other` and `ts`, builds a numpy array with length
    `other.num_nodes` in which the indexes represent the node id in `other` and the
    entries represent the equivalent node id in `ts`. If a node in `other` has no
    equivalent in `ts`, then the entry takes the value `tskit.NULL` (-1). The
    matching is done by comparing the IDs assigned by SLiM which are kept in
    node metadata. This matching of SLiM IDs is *only* done for nodes with time
    older than the specified `split_time`.
    """
    node_mapping = np.full(other.num_nodes, tskit.NULL)
    sids0 = np.array([n.metadata["slim_id"] for n in ts.nodes()])
    sids1 = np.array([n.metadata["slim_id"] for n in other.nodes()])
    alive_before_split1 = (other.tables.nodes.time >= split_time)
    is_1in0 = np.isin(sids1, sids0)
    both = np.logical_and(alive_before_split1, is_1in0)
    sorted_ids0 = np.argsort(sids0)
    matches = np.searchsorted(
        sids0,
        sids1[both],
        side='left',
        sorter=sorted_ids0
    )
    node_mapping[both] = sorted_ids0[matches]
    return node_mapping
```

Now we are finally ready to **union** our tree sequences. For that, I wrote a
recursive function that goes through our data frame with the phylogeny and
returns a dictionary with the merged tree sequences from the tip to the root.

```{code-cell}
def union_recursion(df, focal="root", merged=None):
    if merged is None:
        merged = {
            row.child : (
                tskit.load(row.outfile),
                row.edgelen,
                [row.child]
            ) for i, row in df[df.is_leaf].iterrows()
        }
    print(f"Going in: {focal}")
    childs = df[df.parent == focal]
    assert (len(childs) == 2) or (len(childs) == 0)
    if len(childs) == 2:
        for i, row in childs.iterrows():
            merged = union_recursion(df, row.child, merged)
        cname1 = childs.iloc[0,]["child"]
        cname2 = childs.iloc[1,]["child"]
        split_time = merged[cname1][1]
        assert split_time == merged[cname2][1] # ultrametric
        print(f'Unioning: {childs["child"].to_list()}, Split time: {split_time}')
        ts1 = merged[cname1][0]
        ts2 = merged[cname2][0]
        node_map = match_nodes(ts2, ts1, split_time)
        tsu = ts1.union(ts2, node_map, check_shared_equality=True)
        # the time from tip to start of simulation is split_time plus the
        # length of the edge
        merged[focal] = (
            tsu,
            split_time + df[df.child==focal].edgelen.item(),
            merged[cname1][2] + merged[cname2][2]
        )
    return merged

merged = union_recursion(df)
# union of all three species tree sequences is in the root.
tsu, _, pops = merged["root"]
```

A slightly tricky thing we had to do there was to make sure we kept track of
which population in the union'ed tree sequence corresponds to
which population in our phylogeny.
That's the role of the ``pops`` variable, that
contains the names of the populations
in the order in which they are indexed in the resulting tree sequence ``tsu``.
(Since the population in SLiM is ``p1``, the zero-th population will be unused.)

Let's make sure we have the right number of samples in each of the populations
we specified.

```{code-cell}
for i, name in enumerate(pops):
    n_samples = len(tsu.samples(i + 1)) // 2
    print(f"Union-ed tree sequence has {n_samples} samples in population {name},\n"
          f"\tand we specified {df[df.child==name].popsize.item()} individuals in our simulations.")
```

Finally, we should recapitate the result,
in case some trees on the root branch haven't coalesced,
and write out the result:
```{code-cell}
tsu = pyslim.SlimTreeSequence(tsu)
tsu = tsu.recapitate(recombination_rate=1e-8)
tsu.dump("final.trees")
```

Now we're done, and can analyse the final tree sequence!
Just for fun, I'll look at the trees produced by the simulation.
For instance, we might be curious how often there are
disagreements between the species tree and the simulated gene trees
(also called incomplete lineage sorting, or ILS).

To make it possible to look at the trees,
I will first simplify the union-ed tree sequence to keep only two diploid
samples per population. Then, I'll recapitate the tree sequence so that we have
a single root in all the trees.

```{code-cell}
rng = np.random.default_rng(seed=123)
ind_alive = tsu.individuals_alive_at(0)
ind_pops = tsu.individual_populations[ind_alive]
subsample_indivs = [
    rng.choice(ind_alive[ind_pops == i + 1], 2)
    for i, _ in enumerate(pops)
]
subsample_nodes = [
    np.concatenate([tsu.individual(i).nodes for i in x])
    for x in subsample_indivs
]
tsus = pyslim.SlimTreeSequence(
    tsu.simplify(
        np.concatenate(subsample_nodes),
        filter_populations=False,
    )
)
SVG(tsus.draw_svg(
    node_labels={
        node.id: pops[node.population - 1]
        for node in tsus.nodes()
        if not node.time > 0.0
    },
    x_lim=[0,2200],
    size=(800, 300),
))
```

Again, we have to index ``pops`` by ``node.population - 1``
because the zero-th population is unused.
In an upcoming version of SLiM, populations will have a ``name`` property
in metadata, so we will not have to keep track of things
in this awkward way.
