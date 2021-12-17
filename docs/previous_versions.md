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

ts = tskit.load("example_sim.trees")
tables = ts.tables
```


(sec_previous_versions)=


# Migrating from previous versions of pyslim

A number of features that were first introduced in pyslim have been made part of core
tskit functionality. For instance, reference sequence support was provided (although
loosely) inpyslim to support SLiM's nucleotide models, but is now part of a standard
tskit {class}`tskit.TreeSequence`. Similarly, metadata processing in tskit made
code to do this within pyslim obsolete; this "legacy metadata" code has been removed
and instructions for how to migrate your code are {ref}`below <sec_legacy_metadata>`.

In fact, we are now at the (very good) place where we don't really need
the {class}`pyslim.SlimTreeSequence` class any longer,
and it will soon be deprecated.
So, pyslim is migrating to be purely functional: instead of providing the SlimTreeSequence
class with specialized methods, all methods will be functions of TreeSequences,
that take in a tree sequence and return something
(a modified tree sequence or some summary of it).
Backwards compatibility will be maintained for some time, but we request that you
switch over sooner, as your code will be cleaner and faster.

To migrate, you should:


1. Replace `ts.slim_generation` with `ts.metadata['SLiM']['generation']`,
    and `ts.model_type` with `ts.metadata['SLiM']['model_type']`.
2. Replace `ts.reference_sequence` with `ts.reference_sequence.data`.
3. Replace calls to `ts.recapitate(...)` with `pyslim.recapitate(ts, ...)`,
    and similarly with other SlimTreeSequence methods.

If you encounter difficulties, please post an
[issue](https://github.com/tskit-dev/pyslim/issues)
or [discussion](https://github.com/tskit-dev/pyslim/discussions) on github.


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

.. note::

    Until pyslim version 0.600, the old-style metadata was still available,
    but this functionality has been removed.

Here are more detailed notes on how to migrate a script from the legacy
metadata handling. If you run into issues, please ask (open a discussion on github).

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
Now, this could be done as follows (where now ``metadata`` is a list of metadata dicts):

```{code-cell}
metadata = [ {'slim_id': k, 'is_null': False, 'genome_type': 0}
            for k in range(tables.nodes.num_rows) ]
nms = tables.nodes.metadata_schema
tables.nodes.packset_metadata(
  [nms.validate_and_encode_row(r) for r in metadata]
)
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
=======
>>>>>>> 483184a (deprecation start)
