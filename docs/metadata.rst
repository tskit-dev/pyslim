.. _sec_metadata:

========
Metadata
========


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

Tools for working with tree sequences are provided by a package called `tskit`.
`pyslim` provides a thin interface between `msprime/tskit`.

++++++++++++++++
Metadata entries
++++++++++++++++

SLiM records additional information in the metadata columns of Population, Individual, Node, and Mutation tables.
The information is recorded in a binary format, and is extracted and written by ``pyslim`` using the python ``struct`` module.
Nothing besides this binary information can be stored in the metadata of these tables for the tree sequence to be used by SLiM,
and so when ``pyslim`` annotates an existing tree sequence, anything in those columns is overwritten.
For more detailed documentation on the contents and format of the metadata, see the SLiM manual.
