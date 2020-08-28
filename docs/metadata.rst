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
which can be found in the metadata: `ts.metadata['SLiM']['generation']`
(or, the `ts.slim_generation` attribute).
In most cases, the "SLiM time" referred to by a ``time`` in the tree sequence
(i.e., the value that would be reported by ``sim.generation``
within SLiM at the point in time thus referenced)
can be obtained by subtracting ``time`` from ``ts.slim_generation``.
**However,** in WF models, birth happens between the "early()" and "late()" stages,
so if the tree sequence was written out using ``sim.treeSeqOutput()`` during "early()" in a WF model,
the tree sequence's times measure time before the last set of individuals are born,
i.e., before SLiM time step ``ts.slim_generation - 1``.
The stage that the tree sequence was saved is recorded in the metadata of the tree sequence,
as `ts.metadata['SLiM']['stage']`.
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

To *modify* the metadata that ``pyslim`` has introduced into a coalescent simulation,
or the metadata in a SLiM-produced tree sequence, use the ``annotate_X_metadata()`` functions.
For instance, to set the ages of the individuals in the tree sequence to random numbers between 1 and 4,
and write out the resulting tree sequence:

.. code-block:: python

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
