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
import random
random.seed(23)

ts = tskit.load("example_sim.trees")
tables = ts.tables
```

```{eval-rst}
.. currentmodule:: pyslim
```


(sec_metadata)=

# Metadata

(sec_metdata_overview)=

## Overview

SLiM puts SLiM-specific information into the *metadata* for the tree sequence,
as well as for each populations, individuals, nodes and mutations.
Here is a quick reference to what information is available:
see the SLiM manual for the more technical writeup.
A good way to get a generic metadata example is with {func}`.default_slim_metadata`.

**Top-level:**
If `ts` is your tree sequence, then `ts.metadata` is a dict,
and `ts.metadata["SLiM"]` contains information about the simulation:

- `file_version`: the version of the SLiM tree sequence file format
- `tick`: the value of `community.tick` within SLiM when the file was written out
- `cycle`: the value of `sim.cycle` within SLiM when the file was written out
- `model_type`: either `"WF"` or `"nonWF"`
- `nucleotide_based`: whether this is a nucleotide-based simulation
- `separate_sexes`: whether the simulation has separate sexes or not
- `spatial_dimensionality`: for instance, `""` or `"x"` or `"xy"` (etcetera)
- `spatial_periodicity`: whether space wraps around in some directions (same format as dimensionality)
- `stage`: the *stage* of the life cycle at which the file was written out (either `"first"`, `"early"`, or `"late"`)

**Populations:**
Information about each SLiM-produced population is written to metatadata.
The format uses JSON and is extensible, so other keys may be present
and some keys may be missing (for instance, there are no spatial bounds
in a nonspatial simulation). The metadata may be `None` for populations
that SLiM did not use. The keys that SLiM uses are:

- `slim_id`: the ID of this population in SLiM
- `name`: the name of the population (by default, `p0`, `p1`, etcetera)
- `description`: a string describing the population
- `selfing_fraction`, `female_cloning_fraction`, `male_cloning_fraction`, and `sex_ratio`: only present when applicable (e.g., in WF simulations)
- `bounds_x0`, `bounds_x1`, `bounds_y0`, `bounds_y1`, `bounds_z0`, and `bounds_z1`: the spatial bounds, when applicable
- `migration_records`: A *list* of entries decribing migration between populations in a WF model.

**Individuals:**
Each individual produced by SLiM contains the following metadata:

- `pedigree_id`: the "pedigree ID", unique within the SLiM simulation
- `pedigree_p1`, `pedigree_p2`: the pedigree IDs of the individuals' two
  parents (they may be equal in the case of selfing, or `-1` to indicate no
  parent, in the case of the initial generation or for cloning)
- `age`: the `.age` property within SLiM at the time the file was written out
- `subpopulation`: the subpopulation within SLiM the individual was in at the time the file was written out
- `sex`: the sex of the individual (either {data}`.INDIVIDUAL_TYPE_FEMALE`, {data}`.INDIVIDUAL_TYPE_MALE`, or {data}`.INDIVIDUAL_TYPE_HERMAPHRODITE`)
- `flags`: additional information; currently only recording whether the individual was a "migrant" or not (see the SLiM manual)

**Nodes:**
Each "node" produced by SLiM (i.e., "genome" within SLiM) has:

- 'slim_id': the unique ID associated with the genome by SLiM
- 'is_null': whether the genome is a "null" genome (in which case it isn't
  really there, so shouldn't have any mutations or relationships in the tree
  sequence!)
- 'genome_type': the 'type' of this genome (0 for autosome, 1 for X, 2 for Y)

**Mutations:**
Each mutation's metadata is a dictionary with a single key, `"mutation_list"`,
whose entry is a *list* of metadata dictionaries corresponding to the mutations that are "stacked",
i.e., all present, in all genomes inheriting from this (tskit) mutation.
So, `ts.mutation(12).metadata["mutation_list"]` is a list, each of whose entries contains:

- `mutation_type`: the numeric ID of the `MutationType` within SLiM
- `selection_coeff`: the selection coefficient
- `subpopulation`: the numeric ID of the subpopulation the mutation occurred in
- `slim_time`: the value of `community.tick` when the mutation occurred
- `nucleotide`: either `-1` if there is no associated nucleotide, or the numeric code for the nucleotide (see {data}`.NUCLEOTIDES`)


(sec_metadata_tools)=

## Metadata tools

The dictionaries describing the schema for these metadata entries
are available in {data}`.slim_metadata_schemas`.
Furthermore, this method may be useful in working with metadata:

```{eval-rst}
.. autofunction:: default_slim_metadata
```


## Modifying SLiM metadata
For more on working with metadata,
see {ref}`tskit's metadata documentation <tskit:sec_metadata>`.


### Top-level metadata

The entries of the top-level metadata dict are *read-only*.
So, you might think that
`tables.metadata["SLiM"]["model_type"] = "nonWF"`
would switch the model type,
but this in fact (silently) does nothing. To modify the top-level metadata,
we must (a) work with tables (as tree sequences are immutable, and (b)
extract the metadata dict, modify the dict, and copy it back in.
Instead, you should do
```{code-cell}
md = tables.metadata
md["SLiM"]["model_type"] = "nonWF"
tables.metadata = md
```
Modifying the top-level metadata
could be used to set spatial bounds on an annotated msprime simulation, for instance.
(This is recorded in the population metadata.)


### Modifying SLiM metadata in tables


To modify the metadata that ``pyslim`` has introduced into
the tree sequence produced by a coalescent simulation,
or the metadata in a SLiM-produced tree sequence,
we need to edit the TableCollection that forms the editable data behind the tree sequence.
For instance, to set the ages of the individuals in the tree sequence to random numbers between 1 and 4,
we will extract a copy of the underlying tables, clear it,
and then iterate over the individuals in the tree sequence,
as we go re-inserting them into the tables
after replacing their metadata with a modified version:

```{code-cell}
tables = ts.dump_tables()
tables.individuals.clear()
for ind in ts.individuals():
    md = ind.metadata
    md["age"] = random.choice([1,2,3,4])
    _ = tables.individuals.append(
        ind.replace(metadata=md)
    )

mod_ts = tables.tree_sequence()

# check that it worked:
print("First ten ages:", [mod_ts.individual(i).metadata["age"] for i in range(10)])
for ind in mod_ts.individuals():
    assert ind.metadata['age'] in [1, 2, 3, 4]

# save out the tree sequence
mod_ts.dump("modified_ts.trees")
```

## Technical details

### Metadata entries

SLiM records additional information in the metadata columns of Individual, Node, and Mutation tables,
in a binary format using the python ``struct`` module.
See {ref}`tskit's metadata documentation <tskit:sec_metadata>`
for details on how this works.
Nothing besides this binary information can be stored in the metadata of these tables if the tree sequence is to be used by SLiM,
and so when ``pyslim`` annotates an existing tree sequence, anything in those columns is overwritten.
Population metadata is stored as JSON, however, which is more flexible.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.

Of particular note is that *nodes* and *populations* may have empty metadata.
SLiM will not use the metadata of nodes that are not associated with alive individuals,
so this can safely be omitted (and makes recapitation easier).
And, populations not used by SLiM will have empty metadata.
All remaining metadata are required (besides edges and sites, whose metadata is not used at all).
