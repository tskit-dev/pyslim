.. _sec_tutorial:

========
Tutorial
========


Coalescent simulation is more limited in the degree of biological realism
it can attain, but is much faster than
forward simulation, so it can be helpful to run a "hybrid" simulation, by
endowing a SLiM simulation with a history derived from an msprime coalescent simulation.
For instance, suppose we have a SLiM simulation of a population of 100,000 individuals
that we have run for 10,000 generations without neutral mutations.
Now, we wish to extract whole-genome genotype data for only 1,000 individuals.
A typical way to do this would be to


1. :meth:`.SlimTreeSequence.recapitate` :
   The simulation has likely not reached demographic equilibrium - it has not
   *coalesced* entirely; recapitation uses coalescent simulation to provide
   a "prior history" for the initial generation of the simulation.

2. :meth:`.SlimTreeSequence.simplify` : For efficiency, subset the tree
   sequence to only the information relevant for those 1,000 individuals
   we wish to sample. This needs to come *after* recapitation (see below).

3. :meth:`msprime.mutate` : This adds neutral mutations on top of the tree sequence.

These steps are described below.
First, to get something to work with,
you can run this simple SLiM script
of a single population of sexual organisms,
fluctuating around 1000 individuals,
for 1000 generations:

.. literalinclude:: example_sim.slim

(Note: by setting the random seed in the simulation,
you should get exactly the same results as me,
when you run the code below.)


************
Recapitation
************

.. figure:: _static/pedigree_recapitate.png
   :scale: 42%
   :align: right

Although we can initialize a SLiM simulation with the results of a coalescent simulation,
if during the simulation we don't actually use the genotypes for anything, it can be much
more efficient to only coalesce the portions of the first-generation ancestors that have
not yet coalesced. (See the SLiM manual for more explanation.)
This is depicted in the figure at the right:
imagine that the common ancestors of all samples did not exist at all sites within the
SLiM simulation. Recapitation starts at the *top* of the genealogies,
and runs a coalescent simulation back through time
to fill out the rest of genealogical history relevant to the samples.
The purple nodes are new ancestral genomes that have been added to the tree sequence.
This is important - if we did not do this,
then effectively the initial population would be genetically homogeneous,
and so our simulation would have less genetic variation than it should have
(since the component of variation from the initial population would be omitted).

Doing this is as simple as:

.. code-block:: python

   orig_ts = pyslim.load("example_sim.trees")
   rts = orig_ts.recapitate(recombination_rate = 1e-8, Ne=200)

We can check that this worked as expected, by verifying that after recapitation
all trees have only one root:

.. code-block:: python

   orig_max_roots = max(t.num_roots for t in orig_ts.trees())
   recap_max_roots = max(t.num_roots for t in rts.trees())
   print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
         f"and after recapitation, it was {recap_max_roots}.")
   # Before recapitation, the max number of roots was 15, and after recapitation, it was 1.

Note that demography needs to be set up explicitly - if you have more than one population,
you must set migration rates or else coalescence will never happen
(see below for an example, and :meth:`.SlimTreeSequence.recapitate` for more).


**************
Simplification
**************

.. figure:: _static/pedigree_simplify.png
   :scale: 42%
   :align: right

In this example, we imagine that we have many more individuals in the simulation
than we actually want to analyze.
We can get rid of the extra information using an operation called *simplification*.
This is depicted in the figure at the right:
we have only retained information relevant to the genealogies of the remaining samples,
substantially simplifying the tree sequence.
(Precisely, simplification retains only nodes of the tree sequence that are
branching points of some marginal genealogy -- see
`Kelleher et al 2018 <https://doi.org/10.1371/journal.pcbi.1006581>`_
for details.)
While simplification sounds very appealing - it makes things simpler after all -
it is often not necessary in practice, because tree sequences are very compact,
and many operations with them are quite fast.
So, you should probably not make simplification a standard step in your workflow,
only using it if necessary.

Simplification to history of 100 individuals alive today
can be done with the :meth:`.SlimTreeSequence.simplify` method:

.. code-block:: python

   keep_indivs = np.random.choice(rts.individuals_alive_at(0), 100, replace=False)
   keep_nodes = []
   for i in keep_indivs:
      keep_nodes.extend(rts.individual(i).nodes)
   sts = rts.simplify(keep_nodes)

   print(f"Before, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals) "
          f"in the tree sequence, and now there are {sts.num_samples} sample nodes "
          f"(and {sts.num_individuals} individuals).")
   # Before, there were 1930 sample nodes (and 1965 individuals) in the tree sequence,
   # and now there are 200 sample nodes (and 129 individuals).

**Note** that you must pass simplify a list of *node IDs*, not individual IDs.
Here, we used the :meth:`.SlimTreeSequence.individuals_alive_at` method to obtain the list
of individuals alive today.
Also note that there are *still* more than 100 individuals remaining - 29 non-sample individuals
have not been simplified away, from either the initial population or the final population,
because they have nodes that are required to describe the genealogies of the samples.


*********************************************
Adding neutral mutations to a SLiM simulation
*********************************************

.. figure:: _static/pedigree_mutate.png
   :scale: 42%
   :align: right

If you have recorded a tree sequence in SLiM, likely you have not included any neutral mutations,
since it is much more efficient to simply add these on afterwards.
To add these (in a completely equivalent way to having included them during the simulation),
you can use the :meth:`msprime.mutate` function, which returns a new tree sequence with additional mutations.
Continuing with the cartoons from above, these are added to each branch of the tree sequence
at the rate per unit time that you request.
This works as follows:

.. code-block:: python

   ts = pyslim.SlimTreeSequence(msprime.mutate(sts, rate=1e-8, keep=True))

   print(f"The tree sequence now has {ts.num_mutations} mutations, "
         f"and mean pairwise nucleotide diversity is {ts.diversity()}.")
   # The tree sequence now has 28430 mutations, and mean pairwise nucleotide diversity is 2.3319e-05.


This adds infinite-sites mutations at a rate of 1e-8 per site, making sure to
``keep`` any existing mutations.
We have wrapped the call to :meth:`msprime.mutate` in a call to
:class:`pyslim.SlimTreeSequence`, because :meth:`msprime.mutate` returns an *msprime* tree sequence,
and by converting it back into a ``pyslim`` tree sequence we can still use the methods
defined by ``pyslim``. (The conversion does not modify the tree sequence at all,
it only adds the ``.slim_generation`` attribute.) The output of other ``msprime``
functions that return tree sequences may be converted back to
:class:`pyslim.SlimTreeSequence` in the same way.


.. _sec_extracting_individuals:

**************************************
Extracting particular SLiM individuals
**************************************

To get another example with discrete subpopulations,
let's run another SLiM simulation, similar to the above
but with two populations exchanging migrants:

.. literalinclude:: migrants.slim

The first, most common method to extract individuals is simply to
get all those that were alive at a particular time, using
:meth:`.SlimTreeSequence.individuals_alive_at()`.
For instance, to get the list of individual IDs of
all those alive at the end of the simulation
(i.e., zero time units ago), we could do:

.. code-block:: python

   orig_ts = pyslim.load("migrants.trees")
   alive = orig_ts.individuals_alive_at(0)
   
   print(f"There are {len(alive)} individuals alive from the final generation.")
   # There are 2012 individuals alive from the final generation.

These are individual IDs, and we can use ``ts.individual( )`` to get information
about each of these individuals from their ID.
For instance,
to then count up how many of these individuals are in each population,
we could do:

.. code-block:: python

   num_alive = [0 for _ in range(orig_ts.num_populations)]
   for i in alive:
      ind = orig_ts.individual(i)
      num_alive[ind.population] += 1

   for pop, num in enumerate(num_alive):
      print(f"Number of individuals in population {pop}: {num}")

   # Number of individuals in population 0: 0
   # Number of individuals in population 1: 1014
   # Number of individuals in population 2: 998

Our SLiM script started numbering populations at 1, while tskit starts counting at 0,
so there is an empty "population 0" in a SLiM-produced tree sequence.

Now, let's recapitate and mutate the tree sequence.
Recapitation takes a bit more thought, because we have to specify a migration matrix
(or else it will run forever, unable to coalesce).

.. code-block:: python

   pop_configs = [msprime.PopulationConfiguration(initial_size=1000)
                  for _ in range(orig_ts.num_populations)]
   rts = orig_ts.recapitate(population_configurations=pop_configs,
                            migration_matrix=[[0.0, 0.0, 0.0],
                                              [0.0, 0.0, 0.1],
                                              [0.0, 0.1, 0.0]])
   ts = pyslim.SlimTreeSequence(
            msprime.mutate(rts, rate=1e-8))

Again, there are *three* populations because SLiM starts counting at 1;
the first population is unused (no migrants can go to it).
Let's compute genetic diversity within and between each of the two populations
(we compute the mean density of pairwise nucleotide differences,
often denoted :math:`\pi` and :math:`d_{xy}`).
To do this, we need to extract the node IDs from the individuals of the two populations
that are alive at the end of the simulation.

.. code-block:: python

   pop_indivs = [[], [], []]
   pop_nodes = [[], [], []]
   for i in ts.individuals_alive_at(0):
      ind = ts.individual(i)
      pop_indivs[ind.population].append(i)
      pop_nodes[ind.population].extend(ind.nodes)

   diversity = ts.diversity(pop_nodes[1:])
   divergence = ts.divergence(pop_nodes[1:], indexes=[(0,1)])

   print(f"There are {ts.num_mutations} mutations across {ts.num_trees} distinct "
         f"genealogical trees describing relationships among {ts.num_samples} "
         f"sampled genomes, with a mean genetic diversity of {diversity[0]} and "
         f"{diversity[1]} within the two populations, and a mean divergence of "
         f"{divergence[0]} between them.")
   # There are 117453 mutations across 28972 distinct genealogical trees describing relationships
   # among 4024 sampled genomes, with a mean genetic diversity of 9.697e-05 and 9.661e-05 within
   # the two populations, and a mean divergence of 9.738e-05 between them.


Other information about individuals is available,
for instance, spatial location, SLiM's internal pedigree ID,
their age at time of death (or the end of the simulation, if they are still alive),
their sex, and their birth time.
For instance, we can create an age distribution by sex:

.. code-block:: python

   max_age = max([ind.metadata.age for ind in ts.individuals()])
   age_table = np.zeros((max_age + 1, 2))
   age_labels = {pyslim.INDIVIDUAL_TYPE_FEMALE: 'females',
                 pyslim.INDIVIDUAL_TYPE_MALE: 'males'}
   for i in ts.individuals_alive_at(0):
      ind = ts.individual(i)
      age_table[ind.metadata.age, ind.metadata.sex] += 1

   print(f"number\t{age_labels[0]}\t{age_labels[1]}")
   for age, x in enumerate(age_table):
      print(f"{age}\t{x[0]}\t{x[1]}")

   # number	females	males
   # 0	        337.0	350.0
   # 1	        228.0	209.0
   # 2	        149.0	132.0
   # 3	        109.0	96.0
   # 4	        62.0	70.0
   # 5	        37.0	47.0
   # 6	        23.0	34.0
   # 7	        23.0	24.0
   # 8	        11.0	16.0
   # 9	        9.0	10.0
   # 10	        5.0	5.0
   # 11	        3.0	4.0
   # 12	        2.0	1.0
   # 13	        1.0	3.0
   # 14	        2.0	0.0
   # 15	        2.0	3.0
   # 16	        2.0	0.0
   # 17	        0.0	0.0
   # 18	        0.0	2.0
   # 19	        1.0	0.0

We have looked up how to interpret the `sex` attribute
by using the values of :data:`.INDIVIDUAL_TYPE_FEMALE` (which is 0)
and :data:`.INDIVIDUAL_TYPE_MALE` (which is 1).
In a simulation without separate sexes,
all individuals would have sex equal to :data:`.INDIVIDUAL_TYPE_HERMAPHRODITE`
(which is -1).

The "flags" object also tells us the reasons that each individual
was retained in the tree sequence: whether they were there in the initial population
(the "first generation"), whether they were Remembered, or whether they were alive at the
end of the simulation. (Note these are not mutually exclusive.)
To count these up, we could do this:

.. code-block:: python

   indiv_types = {"first_gen" : 0,
                  "remembered" : 0,
                  "alive" : 0}
   for ind in ts.individuals():
      if ind.flags & pyslim.INDIVIDUAL_FIRST_GEN:
         indiv_types['first_gen'] += 1
      if ind.flags & pyslim.INDIVIDUAL_REMEMBERED:
         indiv_types['remembered'] += 1
      if ind.flags & pyslim.INDIVIDUAL_ALIVE:
         indiv_types['alive'] += 1

   for k in indiv_types:
      print(f"Number of individuals that are {k}: {indiv_types[k]}")

   # Number of individuals that are first_gen: 2000
   # Number of individuals that are remembered: 0
   # Number of individuals that are alive: 2012


As a final example,
suppose that we want to randomly sample 10 individuals alive and older than 2 time steps
from each of the populations at the end of the simulation,
and simplify the tree sequence to retain only those individuals.
We could do this by iterating over individuals as above,
but the numpy arrays :attr:`.SlimTreeSequence.individual_ages`
and :attr:`.SlimTreeSequence.individual_populations` will make this easier:

.. code-block:: python

   alive = ts.individuals_alive_at(0)
   adults = alive[ts.individual_ages[alive] > 2]
   pops = [np.where(
              ts.individual_populations[adults] == k)[0] for k in [1, 2]]
   sample_inds = [np.random.choice(pop, 10, replace=False) for pop in pops]
   sample_nodes = []
   for samp in sample_inds:
      for i in samp:
         sample_nodes.extend(ts.individual(i).nodes)
   sub_ts = ts.simplify(sample_nodes)

The resulting tree sequence does indeed have fewer individuals and fewer trees:

.. code-block:: python

   print(f"There are {sub_ts.num_mutations} mutations across {sub_ts.num_trees} distinct "
         f"genealogical trees describing relationships among {sub_ts.num_samples} "
         f"sampled genomes, with a mean overall genetic diversity of {sub_ts.diversity()}.")
   # There are 46163 mutations across 6376 distinct genealogical trees describing relationships
   # among 40 sampled genomes, with a mean overall genetic diversity of 9.789e-05.



*********************
Reading SLiM metadata
*********************

Each ``Mutation``, ``Population``, ``Node``, and ``Individual`` carries additional information
stored by SLiM in its ``metadata`` property. The precise metadata stored in each is detailed in the SLiM manual.
For instance, here is the information available about an individual
in the previous example:

.. code-block:: python

    print(ts.individual(0))

   # {'id': 0, 'flags': 65536,
   #  'location': array([0., 0., 0.]),
   #  'metadata': IndividualMetadata(
   #                pedigree_id=985545,
   #                age=16,
   #                population=1,
   #                sex=0,
   #                flags=0),
   #  'nodes': array([4000, 4001], dtype=int32),
   #  'population': 1,
   #  'time': 16.0}

For a precise description of the metadata fields, see the SLiM manual.
Briefly,
the ``id`` is the ID internal to the tree sequence;
``location`` is the (x,y,z) coordinates of the individual,
``pedigree_id`` is SLiM's internal ID for the individual,
``age`` and ``population`` are their age and population at death, or at the time the simulation stopped if they were still alive,
``sex`` is their sex
(as an integer, one of :data:`.INDIVIDUAL_TYPE_FEMALE`,
:data:`.INDIVIDUAL_TYPE_MALE`, or :data:`.INDIVIDUAL_TYPE_HERMAPHRODITE`),
``nodes`` is an array of the node IDs that represent the genomes of this individual,
and ``time`` is the time (in units of "time ago") that the individual was born.
Several of these are available as numpy arrays,
across all individuals at once:
:attr:`.SlimTreeSequence.individual_locations`,
:attr:`.SlimTreeSequence.individual_populations`,
:attr:`.SlimTreeSequence.individual_ages`,
and :attr:`.SlimTreeSequence.individual_times`.
Also see :meth:`.SlimTreeSequence.individual_ages_at`.


******************************
Coalescent simulation for SLiM
******************************

The :func:`.annotate` command helps make this easy, by adding default
information to a tree sequence, allowing it to be read in by SLiM. This will
simulate a tree sequence with msprime, add SLiM information, and write it out
to a ``.trees`` file:

.. code-block:: python

   import msprime
   import pyslim

   # simulate a tree sequence of 12 sample genomes
   ts = msprime.simulate(12, mutation_rate=0.0, recombination_rate=1e-8, length=1e6)
   new_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
   new_ts.dump("initialize_nonWF.trees")


Note that we have set the mutation rate to ``0.0``:
this is because any mutations that are produced will be read in by SLiM...
which *could* be a very useful thing, if you want to generate mutations with msprime
that provide standing variation for selection within SLiM...
**but**, currently msprime only produces mutations
with an infinite-sites model, while SLiM requires mutation positions to be at integer positions.
We `plan to fix this <https://github.com/tskit-dev/msprime/issues/553>`_,
but in the meantime you'll have to generate any pre-existing mutations by hand.
*However*, if you intend the pre-existing mutations to be *neutral*,
then there is no need to add them at this point;
you can add them after the fact, as discussed above.
Also note that we have set ``slim_generation`` to 1;
this means that as soon as we load the tree sequence into SLiM,
SLiM will set the current time counter to 1.
(If we set ``slim_generation`` to 100, then any script blocks scheduled to happen before 100
would not execute after loading the tree sequence.)

The resulting file ``slim_ts.trees`` can be read into SLiM to be used as a starting state,
as illustrated in this minimal example::

   initialize()
   {
       initializeSLiMModelType("nonWF");
       initializeTreeSeq();
       initializeMutationRate(1e-2);
       initializeMutationType("m1", 0.5, "f", -0.1);
       initializeGenomicElementType("g1", m1, 1.0);
       initializeGenomicElement(g1, 0, 1e6-1);
       initializeRecombinationRate(1e-8);
   }

   1 early() { 
       sim.readFromPopulationFile("initialize_nonWF.trees");
   }

   10 {
       sim.treeSeqOutput("nonWF_restart.trees");
       catn("Done.");
       sim.simulationFinished();
   }

See the SLiM manual for more about this operation.



**********************************
Possibly important technical notes
**********************************

Also known as "gotchas".

1. If you use msprime to simulate a tree sequence, and then use that to initialize a SLiM simulation,
    you have to specify the same sequence length in both: as in the examples above,
    the ``length`` argument to :py:meth:`msprime.simulate` should be equal to the SLiM sequence length *plus 1.0* (e.g., if the base positions in SLiM are 0 to 99, then there are 100 bases in all,
    so the sequence length should be 100).

2. Make sure to distinguish *individuals* and *nodes*!
   ``tskit`` "nodes" correspond to SLiM "genomes".
   Individuals in SLiM are diploid, so each has two nodes.

3. Since in SLiM, all individual are diploid, every individual will be associated with two nodes.
   As described above, the Individual table contains entries for 

   a. the currently alive individuals, 
   b. any individuals that have been remembered with ``treeSeqRememberIndividuals()``, and
   c. the *first* generation of the SLiM simulation.

   This last category is here because they are necessary for recapitation (described above);
   but they are *not* marked as samples, so will most likely be removed if you `simplify` the tree sequence.


