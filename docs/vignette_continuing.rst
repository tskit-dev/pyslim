.. _sec_vignette_continuing:


======================================================
Vignette: Following up with more coalescent simulation
======================================================

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


******************
Positive selection
******************

Here's a SLiM script that has rapid, strong selection acting genome-wide for 20 generations.
It is perhaps not very realistic, but it's dramatic.

.. literalinclude:: rapid_adaptation.slim

We can see what happened in the GUI,
but let's pull some more statistics out of the tree sequence:

.. code-block:: python

   import pyslim, tskit, msprime
   import numpy as np

   ts = pyslim.load("rapid_adaptation.trees")

   # allele frequencies
   p = ts.sample_count_stat([ts.samples()], lambda x: x/20000, 1, windows='sites', strict=False)
   print(f"There are {ts.num_sites} segregating sites, of which {np.sum(p > 0.25)} "
         f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")

This tells us that:

.. code-block:: none

   There are 60142 segregating sites, of which 95 are at frequency above 25%, and 695 are above 5%.

The selection was, indeed, strong.

************
Recapitation
************

Ok, now let's do phase (1).
This just means recapitating and mutating the result:

.. code-block:: python

   rts = ts.recapitate(Ne=1000, recombination_rate=1e-8, random_seed=6)
   rts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=1e-8, random_seed=7))

   p = rts.sample_count_stat([rts.samples()], lambda x: x/20000, 1, windows='sites', strict=False)
   print(f"There are {rts.num_sites} segregating sites, of which {np.sum(p > 0.25)} "
         f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")

Now, there are more segregating sites - neutral ones.

.. code-block:: none

   There are 81031 segregating sites, of which 243 are at frequency above 25%, and 1311 are above 5%.


*************************
Continuing the simulation
*************************

To "continue" the simulation neutrally, we'll 

a. simulate the desired period of time in msprime
b. randomly match the initial ancestors in the msprime simulation
   with the final individuals of the SLim simulation, and
c. merge the two together, using the :meth:`tskit.TreeSequence.union` method.


*(a)* Simulating for a given period of time in msprime requires the ``end_time`` argument
(remembering that this is *time ago*);
we'll do this to simulate an additional 1000 generations.

This is almost what we need, but there is one more detail:
if complete coalescence occurs on any region of the genome,
msprime will stop simulating the history of that region.
This is a problem, since we need all lineages to extend back to ``end_time``
To make sure all lineages trace back to ``end_time``,
we'll add one "fake" sample from a separate population, that *can't* coalesce with the rest,
then remove it before the next step, using the ``keep_input_roots=True`` argument to ``simplify()``.


.. code-block:: python

   new_time = 1000
   popcons = [msprime.PopulationConfiguration(
                        sample_size=20000, initial_size=10000),
              msprime.PopulationConfiguration(
                        sample_size=1, initial_size=1)]
   new_ts = msprime.simulate(population_configurations=popcons,
                  end_time=new_time,
                  length=rts.sequence_length,
                  mutation_rate=1e-8, recombination_rate=1e-8,
                  random_seed=9)
   new_tables = new_ts.tables
   # check that the spurious sample is number 20000
   assert 20000 in new_ts.samples()
   assert new_ts.node(20000).population == 1
   new_tables.simplify(samples=np.arange(20000), keep_input_roots=True)

*(b)* Now we'll pull out the IDs of the nodes from 1000 generations ago,
shift the times in the SLiM tree sequence back 1000 generations,
randomly assign each to a node at the end of the SLiM simulation,
and merge them.

.. code-block:: python

   new_nodes = np.where(new_tables.nodes.time == new_time)[0]
   print(f"There are {len(new_nodes)} nodes from the start of the new simulation.")
   # There are 4425 nodes from the start of the new simulation.


   slim_indivs = rts.individuals_alive_at(0)
   slim_nodes = []
   for ind in slim_indivs:
      slim_nodes.extend(ts.individual(ind).nodes)

   slim_nodes = np.array(slim_nodes)
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

   p = full_ts.sample_count_stat([full_ts.samples()], lambda x: x/20000, 1,
                                 windows='sites', strict=False)
   print(f"There are {full_ts.num_sites} segregating sites, of which {np.sum(p > 0.25)} "
         f"are at frequency above 25%, and {np.sum(p > 0.05)} are above 5%.")

Don't worry, we'll explain what happened there in a minute.
Now, things have drifted:

.. code-block:: none

   There are 328814 segregating sites, of which 4312 are at frequency above 25%,
   and 20912 are above 5%.

Let's do a consistency check:

.. code-block:: python

   t = rts.first()
   assert(t.num_roots == 1)
   r = rts.node(t.root)
   # Root of the first tree in the recapitated SLiM simulation:
   print(r)
   # {'id': 79305, 'time': 5920.012447565916, 'population': 1, 'individual': -1,
   #  'flags': 0, 'metadata': None}

   ft = full_ts.first()
   assert(ft.num_roots == 1)
   fr = full_ts.node(ft.root)
   # Root of the first tree after continuing for 1000 generations:
   print(fr)
   # {'id': 79305, 'time': 6920.012447565916, 'population': 1, 'individual': -1,
   #  'flags': 0, 'metadata': None}

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
