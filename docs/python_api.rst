.. _sec_python_api:

==========
Python API
==========

This page provides detailed documentation for the ``pyslim`` Python API.

******************************
Additions to the tree sequence
******************************

The :class:`tskit.TreeSequence` class represents a sequence of correlated trees
that describes the genealogical history of a population.
Here, we describe pyslim's additions to this class.


++++++++++++++++++++++
The SLiM Tree Sequence
++++++++++++++++++++++


.. autoclass:: pyslim.SlimTreeSequence()
    :members:


*********
Constants
*********

.. data:: slim_tree_sequence.py:NUCLEOTIDES == ['A', 'C', 'G', 'T']

   Nucleotide states in nucleotide models are encoded as integers (0, 1, 2, 3);
   this gives the mapping.


********
Metadata
********

SLiM-specific metadata is made visible to the user by ``.metadata`` properties.
For instance::

   >>> ts.node(3).metadata
   NodeMetadata(slim_id=3, is_null=0, genome_type=0)

shows that the fourth node in the tree sequence was given pedigree ID ``3`` by SLiM,
is *not* a null genome, and has ``genome_type`` zero, which corresponds to an autosome 
(see below).


+++++++++
Mutations
+++++++++

.. autoclass:: pyslim.MutationMetadata
   :members:

.. autofunction:: pyslim.decode_mutation

.. autofunction:: pyslim.encode_mutation


+++++
Nodes
+++++

.. autoclass:: pyslim.NodeMetadata
   :members:

.. autofunction:: pyslim.decode_node

.. autofunction:: pyslim.encode_node

These constants are used to encode the "genome type" of a Node.

.. data:: GENOME_TYPE_AUTOSOME == 0

.. data:: GENOME_TYPE_X == 1

.. data:: GENOME_TYPE_Y == 2


+++++++++++
Individuals
+++++++++++

.. autoclass:: pyslim.IndividualMetadata
   :members:

.. autofunction:: pyslim.decode_individual

.. autofunction:: pyslim.encode_individual

These constants are used to encode other information about Individuals.


.. data:: INDIVIDUAL_TYPE_HERMAPHRODITE == -1

.. data:: INDIVIDUAL_TYPE_FEMALE == 0

.. data:: INDIVIDUAL_TYPE_MALE == 1

.. data:: INDIVIDUAL_FLAG_MIGRATED == 0x01

And, these are used in the `Individual.flags`:

.. data:: INDIVIDUAL_ALIVE == 2**16

   This flag is used by SLiM to record information in the :class:`tskit.Individual` metadata.

.. data:: INDIVIDUAL_REMEMBERED == 2**17

   This flag is used by SLiM to record information in the :class:`tskit.Individual` metadata.

.. data:: INDIVIDUAL_FIRST_GEN == 2**18

   This flag is used by SLiM to record information in the :class:`tskit.Individual` metadata.


+++++++++++
Populations
+++++++++++

The :class:`Population` metadata contains quite a bit of information about a SLiM population.

.. autoclass:: pyslim.PopulationMigrationMetadata
   :members:

.. autoclass:: pyslim.PopulationMetadata
   :members:

.. autofunction:: pyslim.decode_population

.. autofunction:: pyslim.encode_population


++++++++++
Annotation
++++++++++

These two functions will add default SLiM metadata to a tree sequence (or the
underlying tables), which can then be modified and loaded into SLiM.

.. autofunction:: pyslim.annotate_defaults

.. autofunction:: pyslim.annotate_defaults_tables


+++++++++++
Provenances
+++++++++++

.. autoclass:: pyslim.ProvenanceMetadata
   :members:

.. autofunction:: pyslim.get_provenance

