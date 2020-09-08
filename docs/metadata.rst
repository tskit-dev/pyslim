.. _sec_metadata:

========
Metadata
========

***************************************
Converting from SLiM time to tskit time
***************************************

SLiM is a forwards simulator, while the tree sequence format thinks about things
*retrospectively*, and so works with times in units of *time ago*.
Mostly, you don't have to convert between the two,
unless you want to match up information in a tree sequence
with information written out by SLiM itself.
In other words, SLiM's time counter measures the number of time steps
("generations") since the start of the simulation,
and times in the tree sequence record how long before the end of the simulation.
However, off-by-one errors are easy to make, so we'll spell it out in detail.

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

.. code-block:: python

   def slim_time(ts, time):
      slim_time = ts.slim_generation - time
      time_adjust =  (self.metadata['SLiM']['model_type'] == "WF"
                      and self.metadata['SLiM']['stage'] == "early")
      return slim_time - time_adjust

Some of the other methods in pyslim -- those that depend on :meth:`.SlimTreeSequence.individuals_alive_at`
-- need you to tell them during which stage the tree sequence was saved with ``sim.treeSeqOutput``,
and need this to be the same as the stage that any individuals were saved with ``sim.treeSeqRememberIndividuals``.
This argument, ``remembered_stage``, defaults to "late()";
we recommend that you also default to always Remembering individuals, and saving out the tree sequence,
during "late()" as well, unless you have good reason not to.
(This means you *must specify* the stage of the block in your SLiM script,
since the stage defaults to "early()"!)

***********************
Modifying SLiM metadata
***********************

For more on working with metadata,
see `tskit's metadata documentation <https://tskit.readthedocs.io/en/latest/metadata.html#sec-metadata>`_.

++++++++++++++++++
Top-level metadata
++++++++++++++++++

The entries of the top-level metadata dict are *read-only*: so,
you might think that
``tables.metadata["SLiM"]["model_type"] = "nonWF"`` would switch the model type,
but this in fact (silently) does nothing. To modify the top-level metadata,
we must (a) work with tables (as tree sequences are immutable, and (b)
extract the metadata dict, modify the dict, and copy it back in.
Instead, you should do

.. code-block:: python

   md = tables.metadata
   md["SLiM"]["model_type"] = "nonWF"
   tables.metadata = md

Modifying the top-level metadata
could be used to set spatial bounds on an annotated msprime simulation, for instance.


+++++++++++++++++++++++++++++++++
Modifying SLiM metadata in tables
+++++++++++++++++++++++++++++++++


To modify the metadata that ``pyslim`` has introduced into 
the tree sequence produced by a coalescent simulation,
or the metadata in a SLiM-produced tree sequence,
what we do is (a) extract the metadata (as a list of dicts),
(b) modify them, and then (c) write them back into the tables.
For instance, to set the ages of the individuals in the tree sequence to random numbers between 1 and 4,
and write out the resulting tree sequence:

.. code-block:: python

   import random

   tables = ts.tables
   ind_md = [ind.metadata for ind in tables.individuals]
   for md in ind_md:
       md["age"] = random.choice([1,2,3,4])

   ims = tables.individuals.metadata_schema
   tables.individuals.packset_metadata(
      [ims.validate_and_encode_row(md) for md in ind_md])
   mod_ts = pyslim.load_tables(tables, slim_format=True)

   # check that it worked:
   for ind in mod_ts.individuals():
       print(ind.metadata["age"])

   # save out the tree sequence
   mod_ts.dump("modified_ts.trees")


*****************
Technical details
*****************

++++++++++++++++
Metadata entries
++++++++++++++++

SLiM records additional information in the metadata columns of Population, Individual, Node, and Mutation tables,
in a binary format using the python ``struct`` module.
See `tskit's metadata documentation <https://tskit.readthedocs.io/en/latest/metadata.html#sec-metadata>`_
for details on how this works.
Nothing besides this binary information can be stored in the metadata of these tables if the tree sequence is to be used by SLiM,
and so when ``pyslim`` annotates an existing tree sequence, anything in those columns is overwritten.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.

Of particular note is that *nodes* and *populations* may have empty metadata.
SLiM will not use the metadata of nodes that are not associated with alive individuals,
so this can safely be omitted (and makes recapitation easier).
And, populations not used by SLiM will have empty metadata.
All remaining metadata are required (besides edges and sites, whose metadata is not used at all).


.. _sec_legacy_metadata:

===============
Legacy metadata
===============

In previous versions of pyslim,
SLiM-specific metadata was provided as customized objects:
for instance, for a node ``n`` provided by a ``SlimTreeSequence``,
we'd have ``n.metadata`` as a ``NodeMetadata`` object,
with attributes ``n.metadata.slim_id`` and ``n.metadata.is_null`` and ``n.metadata.genome_type``.
However, with tskit 0.3, 
the capacity to deal with structured metadata
was implemented in `tskit itself <https://tskit.readthedocs.io/en/latest/metadata.html#sec-metadata>`_,
and so pyslim shifted to using the tskit-native metadata tools.
As a result, parsed metadata is provided as a dictionary instead of an object,
so that now ``n.metadata`` would be a dict,
with entries ``n.metadata["slim_id"]`` and ``n.metadata["is_null"]`` and ``n.metadata["genome_type"]``.
Annotation should be done with tskit methods (e.g., ``packset_metadata``).

For now, the old-style metadata is still available:
passing the argument ``legacy_metadata=True`` to :meth:`load`
will produce a tree sequence whose metadata is just as before,
and so all previously-written scripts that depend on metadata processing should work, unchanged.
Restating this:

.. note::

   To make an script that relied on previous metadata parsing work,
   it should suffice to replace ``pyslim.load("file.trees")`` with
   ``pyslim.load("file.trees", legacy_metadata=True)``.
   If this fails, please file an issue on github.

Here are more detailed notes on how to migrate a script from the legacy
metadata handling.

First, switch metadata objects to dicts:
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

Second, the ``decode_X`` and ``encode_X`` methods are now deprecated,
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

Third, the ``annotate_X_metadata`` methods are deprecated,
as again tskit has tools to do this.
These methods would set the metadata column of a table -
for instance, if ``metadata`` is a list of NodeMetadata objects, then
``annotate_node_metadata(tables, metadata)`` would modify ``tables.nodes`` in place
to contain the (encoded) metadata in the list ``metadata``.
Now, this would be done as follows (where now ``metadata`` is a list of metadata dicts):

.. code-block:: python

   nms = tables.nodes.metadata_schema
   tables.nodes.packset_metadata(
      [nms.validate_and_encode_row(r) for r in metadata])

If speed is an issue, then ``encode_row`` can be substituted for ``validate_and_encode_row``,
but at the risk of missing errors in metadata.

Fourth, the ``extract_X_metadata`` methods are not necessary,
since the metadata in the tables of a TableCollection are automatically decoded.
For instance, ``[ind.metadata["sex"] for ind in tables.individuals]`` will obtain
a list of sexes of the individuals in the IndividualTable.

.. warning::

   It is our intention to remain backwards-compatible for a time.
   However, the legacy code will disappear at some point in the future,
   so please migrate over scripts you intend to rely on.
