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

ts = pyslim.load("example_sim.trees")
tables = ts.tables
```


(sec_metadata)=

# Metadata

(sec_metadata_converting_times)=

## Converting from SLiM time to tskit time

:::{note}
This is a nitpicky, document-the-details section.
Hopefully, you don't have to deal with the specifics of converting between tskit and SLiM time,
but this page is here for you if you do.
:::

SLiM is a forwards simulator, while the tree sequence format thinks about things
*retrospectively*, and so works with times in units of *time ago*.
Mostly, you don't have to convert between the two,
unless you want to match up information in a tree sequence
with information written out by SLiM itself.
In other words, SLiM's time counter measures the number of time steps
("generations") since the start of the simulation,
and times in the tree sequence record how long before the end of the simulation.
However, there are Some Details, and off-by-one errors are easy to make,
so we'll spell it out in detail.

SLiM's time counter is called the "generation"
(although a "year" or "life cycle" would be a more appropriate name for a nonWF model).
The SLiM generation starts at 1, and records which round of the life cycle the simulation is in.
However, the order of the life cycle differs between WF and nonWF models:
in a WF model, it is "*early* {math}`\to` *birth* {math}`\to` *late*",
while in a nonWF model, it is "*birth* {math}`\to` *early* {math}`\to` *late*".
Usually, the first set of individuals are created in the *early()* phase of generation 1,
and so in a WF model reproduce immediately, in the same generation they were "born".
Parents and offspring cannot have the same birth time in the tree sequence,
and so some clever bookkeeping was required.
You'll want to refer to the tables below to see what's going on.
"Time" in a tree sequence is actually *time ago*,
or *time before the tree sequence was recorded*.
To obtain this number, and ensure that offspring cannot have the same birth time-ago
in the tree sequence as their parents,
SLiM also keeps track of "how many birth phases of the life cycle have happened so far"
(the column "# births" in the tables).
As the simulation goes along,
tskit time ago is recorded as minus one times the number of birth phases so far.
When the tree sequence is output, the current cumulative number of birth phases
is added to this,
so "tskit time ago" is, equivalently, "how many birth phases happened since this time".
In a nonWF model, the two counters ("generation" and "number of birth phases")
are always in sync; but in a WF they are not (during *early*).
The extra wrinkle this introduces is that the correspondence between "tskit time ago"
and "SLiM time" depends on *which phase the tree sequence was recorded in*,
but only for WF models.

To help keep all this straight, here are schematics for WF and nonWF models.
(To see the nonWF model, click on the tab.)

```{tabbed} WF model

For a WF model, the SLiM generation (first column) can be obtained by subtracting the
tskit time ago from the SLiM generation at time of output only during the same stage that output occured in.

|    generation      |       stage         |  # births          |                                |  tskit time ago, early output |   tskit time ago, late output |
|--------------------|---------------------|--------------------|--------------------------------|-------------------------------|-------------------------------|
|       1            |       early         |       0            | {math}`\leftarrow` add subpops |        n-1                    |         n                     |
|       1            |       birth         |       1            |                                |        n-2                    |         n-1                   |
|       1            |       late          |       1            |                                |        n-2                    |         n-1                   |
|       2            |       early         |       1            |                                |        n-2                    |         n-1                   |
|       2            |       birth         |       2            |                                |        n-3                    |         n-2                   |
|       2            |       late          |       2            |                                |        n-3                    |         n-2                   |
|       3            |       early         |       2            |                                |        n-3                    |         n-2                   |
|       3            |       birth         |       3            |                                |        n-4                    |         n-3                   |
|       3            |       late          |       3            |                                |        n-4                    |         n-3                   |
| {math}`\downarrow` | {math}`\cdots`      | {math}`\downarrow` |                                | {math}`\uparrow`              | {math}`\uparrow`              |
|       n-2          |       early         |       n-3          |                                |        2                      |         2                     |
|       n-2          |       birth         |       n-2          |                                |        1                      |         2                     |
|       n-2          |       late          |       n-2          |                                |        1                      |         2                     |
|       n-1          |       early         |       n-2          |                                |        1                      |         2                     |
|       n-1          |       birth         |       n-1          |                                |        0                      |         1                     |
|       n-1          |       late          |       n-1          |                                |        0                      |         1                     |
|       n            |       early         |       n-1          |  treeSeqOutput {math}`\to`     |        0                      |         1                     |
|       n            |       birth         |       n            |                                |                               |         0                     |
|       n            |       late          |       n            |                                | treeSeqOutput {math}`\to`     |         0                     |

```

```{tabbed} nonWF model

Note that for nonWF models the SLiM generation (first column) can always be obtained by subtracting the
tskit time ago from the SLiM generation at time of output.

|    generation      |       stage         |  # births          |                                |  tskit time ago, early output |   tskit time ago, late output |
|--------------------|---------------------|--------------------|--------------------------------|-------------------------------|-------------------------------|
|       1            |       birth         |       1            |                                |        n-1                    |         n-1                   |
|       1            |       early         |       1            | {math}`\leftarrow` add subpops |        n-1                    |         n-1                   |
|       1            |       late          |       1            |                                |        n-1                    |         n-1                   |
|       2            |       birth         |       2            |                                |        n-2                    |         n-2                   |
|       2            |       early         |       2            |                                |        n-2                    |         n-2                   |
|       2            |       late          |       2            |                                |        n-2                    |         n-2                   |
|       3            |       birth         |       3            |                                |        n-3                    |         n-3                   |
|       3            |       early         |       3            |                                |        n-3                    |         n-3                   |
|       3            |       late          |       3            |                                |        n-3                    |         n-3                   |
| {math}`\downarrow` | {math}`\cdots`      | {math}`\downarrow` |                                | {math}`\uparrow`              | {math}`\uparrow`              |
|       n-2          |       birth         |       n-2          |                                |        2                      |         2                     |
|       n-2          |       early         |       n-2          |                                |        2                      |         2                     |
|       n-2          |       late          |       n-2          |                                |        2                      |         2                     |
|       n-1          |       birth         |       n-1          |                                |        1                      |         1                     |
|       n-1          |       early         |       n-1          |                                |        1                      |         1                     |
|       n-1          |       late          |       n-1          |                                |        1                      |         1                     |
|       n            |       birth         |       n            |                                |        0                      |         0                     |
|       n            |       early         |       n            |  treeSeqOutput {math}`\to`     |        0                      |         0                     |
|       n            |       late          |       n            |                                | treeSeqOutput {math}`\to`     |         0                     |

```
When the tree sequence is written out, SLiM records the value of its current generation,
which can be found in the metadata: ``ts.metadata['SLiM']['generation']``
(or, the ``ts.slim_generation`` attribute).
In most cases, the "SLiM time" referred to by a ``time`` in the tree sequence
(i.e., the value that would be reported by ``sim.generation``
within SLiM at the point in time thus referenced)
can be obtained by subtracting ``time`` from ``ts.slim_generation``.
**However,** in WF models, birth happens between the "early()" and "late()" stages,
so if the tree sequence was written out using ``sim.treeSeqOutput()`` during "early()" in a WF model,
the tree sequence's times measure time before the last set of individuals are born,
i.e., before SLiM time step ``ts.slim_generation - 1``.
The stage that the tree sequence was saved is recorded in the metadata of the tree sequence,
as ``ts.metadata['SLiM']['stage']``.
Using this, we can convert from the times of a tree sequence ``ts``
to SLiM time as follows:

```{code-cell}
def slim_time(ts, time, stage):
  slim_time = ts.slim_generation - time
  if ts.metadata['SLiM']['model_type'] == "WF":
    if (ts.metadata['SLiM']['stage'] == "early"
        and stage == "late"):
        slim_time -= 1
    if (ts.metadata['SLiM']['stage'] == "late"
        and stage == "early"):
        slim_time += 1
  return slim_time
```

This is what is computed by the {meth}`.SlimTreeSequence.slim_time` method
(which also has a ``stage`` argument).

Some of the other methods in pyslim -- those that depend on {meth}`.SlimTreeSequence.individuals_alive_at`
-- need you to tell them during which stage the tree sequence was saved with ``sim.treeSeqOutput``,
and need this to be the same as the stage that any individuals were saved with ``sim.treeSeqRememberIndividuals``.
This argument, ``remembered_stage``, defaults to "late()";
we recommend that you also default to always Remembering individuals, and saving out the tree sequence,
during "late()" as well, unless you have good reason not to.
(This means you *must specify* the stage of the block in your SLiM script,
since the stage defaults to "early()"!)


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
what we do is (a) extract the metadata (as a list of dicts),
(b) modify them, and then (c) write them back into the tables.
For instance, to set the ages of the individuals in the tree sequence to random numbers between 1 and 4,
and write out the resulting tree sequence:

```{code-cell}
tables = ts.tables
ind_md = [ind.metadata for ind in tables.individuals]
for md in ind_md:
   md["age"] = random.choice([1,2,3,4])

ims = tables.individuals.metadata_schema
tables.individuals.packset_metadata(
  [ims.validate_and_encode_row(md) for md in ind_md])
mod_ts = pyslim.load_tables(tables)

# check that it worked:
print("First ten ages:", [mod_ts.individual(i).metadata["age"] for i in range(10)])
for ind in mod_ts.individuals():
    assert ind.metadata['age'] in [1, 2, 3, 4]

# save out the tree sequence
mod_ts.dump("modified_ts.trees")
```

## Technical details

### Metadata entries

SLiM records additional information in the metadata columns of Population, Individual, Node, and Mutation tables,
in a binary format using the python ``struct`` module.
See {ref}`tskit's metadata documentation <tskit:sec_metadata>`
for details on how this works.
Nothing besides this binary information can be stored in the metadata of these tables if the tree sequence is to be used by SLiM,
and so when ``pyslim`` annotates an existing tree sequence, anything in those columns is overwritten.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.

Of particular note is that *nodes* and *populations* may have empty metadata.
SLiM will not use the metadata of nodes that are not associated with alive individuals,
so this can safely be omitted (and makes recapitation easier).
And, populations not used by SLiM will have empty metadata.
All remaining metadata are required (besides edges and sites, whose metadata is not used at all).


(sec_legacy_metadata)=

## Legacy metadata

In previous versions of pyslim,
SLiM-specific metadata was provided as customized objects:
for instance, for a node ``n`` provided by a ``SlimTreeSequence``,
we'd have ``n.metadata`` as a ``NodeMetadata`` object,
with attributes ``n.metadata.slim_id`` and ``n.metadata.is_null`` and ``n.metadata.genome_type``.
However, with tskit 0.3,
the capacity to deal with structured metadata
was implemented in {ref}`tskit itself <tskit:sec_metadata>`,
and so pyslim shifted to using the tskit-native metadata tools.
As a result, parsed metadata is provided as a dictionary instead of an object,
so that now ``n.metadata`` would be a dict,
with entries ``n.metadata["slim_id"]`` and ``n.metadata["is_null"]`` and ``n.metadata["genome_type"]``.
Annotation should be done with tskit methods (e.g., ``packset_metadata``).

For now, the old-style metadata is still available:
passing the argument ``legacy_metadata=True`` to {meth}`load`
will produce a tree sequence whose metadata is just as before,
and so all previously-written scripts that depend on metadata processing should work, unchanged.
Restating this:

:::{note}
   To make an script that relied on previous metadata parsing work,
   it should suffice to add `legacy_metadata=True` to calls producing
   SlimTreeSequences, e.g., replacing ``pyslim.load("file.trees")`` with
   ``pyslim.load("file.trees", legacy_metadata=True)``, and
   ``ts.simplify(nodes)`` with
   ``pyslim.SlimTreeSequence(ts.simplify(nodes), legacy_metadata=True)``.
   If this fails, please file an issue on github.
:::

Here are more detailed notes on how to migrate a script from the legacy
metadata handling.

**1.** Use top-level metadata instead of ``slim_provenance``:
previously, information about the model type and the time counter (generation)
in SLiM was provided in the Provenances table, made available through
the ``ts.slim_provenance`` object.  This is still available but deprecated,
and should be obtained from the *top-level* metadata object, ``ts.metadata["SLiM"]``.
So, in your scripts ``ts.slim_provenance.model_type`` should be replaced with
``ts.metadata["SLiM"]["model_type"]``,
and (although it's not deprecated), probably ``ts.slim_generation`` should
probably be replaced with
``ts.metadata["SLiM"]["generation"]``.

**2.** Switch metadata objects to dicts:
if ``md`` is the ``metadata`` property of a population, individual, or node,
this means replacing ``md.X`` with ``md["X"]``.
The ``migration_records`` property of population metadata is similarly
a list of dicts rather than a list of objects, so instead of
``ts.population(1).metadata.migration_records[0].source_subpop``
we would write
``ts.population(1).metadata["migration_records"][0]["source_subpop"]``.

Mutations were previously a bit different - if ``mut`` is a mutation
(e.g., ``mut = ts.mutation(0)``)
then ``mut.metadata`` was previously a list of MutationMetadata objects.
Now, ``mut.metadata`` is a dict, with a single entry:
``mut.metadata["mutation_list"]`` is a list of dicts, each containing the information
that was previously in the MutationMetadata objects.
So, for instance, instead of ``mut.metadata[0].selection_coeff``
we would write ``mut.metadata["mutation_list"][0]["selection_coeff"]``.

**3.** The ``decode_X`` and ``encode_X`` methods are now deprecated,
as this is handled by tskit itself.
For instance, ``encode_node`` would take a NodeMetadata object
and produce the raw bytes necessary to encode it in a Node table,
and ``decode_node`` would do the inverse operation.
This is now handled by the relevant MetadataSchema object:
for nodes one can obtain this as ``nms = ts.tables.nodes.metadata_schema``,
which has the methods ``nms.validate_and_encode_row`` and ``nms.decode_row``.
Decoding is for the most part not necessary,
since the metadata is automatically decoded,
but ``pyslim.decode_node(raw_md)`` could be replaced by ``nms.decode_row(raw_md)``.
Encoding is necessary to modify tables,
and ``pyslim.encode_node(md)`` can be replaced by ``nms.validate_and_encode_row(md)``
(where furthermore ``md`` should now be a dict rather than a NodeMetadata object).

**4.** The ``annotate_X_metadata`` methods are deprecated,
as again tskit has tools to do this.
These methods would set the metadata column of a table -
for instance, if ``metadata`` is a list of NodeMetadata objects, then
``annotate_node_metadata(tables, metadata)`` would modify ``tables.nodes`` in place
to contain the (encoded) metadata in the list ``metadata``.
Now, this would be done as follows (where now ``metadata`` is a list of metadata dicts):

```{code-cell}
metadata = [ {'slim_id': k, 'is_null': False, 'genome_type': 0}
            for k in range(tables.nodes.num_rows) ]
nms = tables.nodes.metadata_schema
tables.nodes.packset_metadata(
  [nms.validate_and_encode_row(r) for r in metadata])
```

If speed is an issue, then ``encode_row`` can be substituted for ``validate_and_encode_row``,
but at the risk of missing errors in metadata.

**5.** the ``extract_X_metadata`` methods are not necessary,
since the metadata in the tables of a TableCollection are automatically decoded.
For instance, ``[ind.metadata["sex"] for ind in tables.individuals]`` will obtain
a list of sexes of the individuals in the IndividualTable.

:::{warning}
   It is our intention to remain backwards-compatible for a time.
   However, the legacy code will disappear at some point in the future,
   so please migrate over scripts you intend to rely on.
:::
