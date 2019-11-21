.. _sec_introduction:

============
Introduction
============

This is the documentation for pyslim, a Python API
for reading and modifying `tskit <https://tskit.readthedocs.io/>`_ tree sequence files
produced by `SLiM <https://messerlab.org/slim/>`_, 
or modifying files produced by other programs (e.g.,
`msprime <https://msprime.readthedocs.io/en/stable/>`_,
`fwdpy11
<https://fwdpy11.readthedocs.io/en/stable/pages/tsoverview.html>`_
and `tsinfer <https://tsinfer.readthedocs.io/>`_) for use in SLiM. 

SLiM can read and write *tree sequences*, which store genealogical data of entire populations.
These can be used to efficiently store both the state of the population at various points
during a simulation *as well as* its genealogical history. Furthermore, SLiM can "load" a saved tree sequence
file to recreate the exact state of the population at the time it was saved.
To do this, SLiM has added several additional types of information to the basic tree sequence file
(in "metadata"); this package makes it easy to read and write this information.

********
Overview
********

A tree sequence is a way of storing both genealogies and genotypes
of a bunch of genomes.
See `the tskit documentation <https://tskit.readthedocs.io/en/latest/>`_
for more description of the tree sequence and underlying data structure.
Each (haploid) genome is associated with a *node*,
and the "focal" nodes are called *samples*.
Many operations act by default on the samples;
and the tree sequence always describes the entire genome of a sample
(unlike other, ancestral nodes, about which we might have only partial information).
SLiM is diploid, so each *individual* has two nodes;
many operations you might want to do involve first finding the individuals you want,
and then looking at their nodes.

****************************
Who is in the tree sequence?
****************************

Which individuals are present in a SLiM-produced tree sequence?
There are *three* types:

1. Everyone who was alive at the end of the simulation.
   Their nodes are samples.

2. Everyone who you asked SLiM to *remember*,
   using the ``treeSeqRememberIndividuals()`` method.
   Their nodes are also samples.

3. Everyone who was alive at the *start* of the simulation.
   Their nodes are *not* samples.

This last category is a common source of confusion:
why are these individuals there?
This is to allow *recapitation*, described in the :ref:`sec_tutorial`:
we need the first generation to be able to trace ancestry back from.
You can tell which individuals are in which categories
by looking at their *flags*;
these three categories are marked with the ``INDIVIDUAL_ALIVE``,
``INDIVIDUAL_REMEMBERED``, and ``INDIVIDUAL_FIRST_GEN`` flags, respectively.
You can also pull out individuals alive at a particular time
with the :meth:`.SlimTreeSequence.individuals_alive_at()` method.
See examples below or in the Vignette for how to do these things.

************************************************
What else can I find out from the tree sequence?
************************************************

Enough information is stored in the tree sequence
to completely reconstruct the state of the SLiM simulation
(except for user-defined data, like a `tag`.
Most of this is stored as *metadata*, which pyslim makes accessible:
see :ref:`sec_metadata`.

