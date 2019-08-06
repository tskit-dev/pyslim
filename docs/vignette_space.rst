.. _sec_vignette_space:


==============================
Vignette: A spatial simulation
==============================

Here we'll talk through a typical workflow with pyslim,
which will:

1. Simulate data with SLiM, remembering some ancestral individuals.
2. Recapitate and mutate.
3. Take a subsample of the modern and ancestral individuals.
4. Get these individual locations and make a map.
5. Compute divergences between individuals, and plot against geographic distance.
6. Write out a VCF file of these individuals' genotypes and other data for use by other programs.


**********
Simulation
**********

Here is a simple spatial SLiM recipe that simulates 1000 individuals on a spatial landscape.
The focus of this vignette is not on SLiM, so we won't go into detail here.
Here are notes:

1. It does not have *any* mutations: we'll add these on afterwards.
2. There is local fecundity regulation of population density: individuals with more neighbors
   have fewer offspring.
3. We run the simulation for 1000 time steps, and "remember" everyone who is alive at time step 500.

.. code-block:: none

   initialize() {
       setSeed(23);
       initializeSLiMModelType("nonWF");
       initializeSLiMOptions(dimensionality="xy");
       initializeTreeSeq();
       initializeMutationRate(0.0);
       initializeMutationType("m1", 0.5, "f", 0.0);
       initializeGenomicElementType("g1", m1, 1.0);
       initializeGenomicElement(g1, 0, 1e8-1);
       initializeRecombinationRate(1e-8);    

       defineConstant("LAMBDA", 2.0); // birth rate
       defineConstant("K", 1);      // carrying capacity per unit area
       defineConstant("W", 20);      // width and height of the area
       defineConstant("MU", 0.5);     // death rate
       defineConstant("SIGMA", 0.5);  // interaction distance
       
       // spatial interaction for local competition
       initializeInteractionType("i1", "xy", reciprocal=T,
                                 maxDistance = 3 * SIGMA);
       i1.setInteractionFunction("n", 1.0/(2*PI*SIGMA^2), SIGMA);
   }

   reproduction() {
       neighbor_density = i1.totalOfNeighborStrengths(individual);
       num_offspring = rpois(1, LAMBDA / (1 + neighbor_density / K));
       mate = i1.drawByStrength(individual, 1);  // single mating
       if (size(mate) > 0) {
           for (k in seqLen(num_offspring)) {
               offspring = p1.addCrossed(individual, mate);
               pos = individual.spatialPosition + rnorm(2, 0, SIGMA);
               offspring.setSpatialPosition(p1.pointReflected(pos));
           }
       }
   }

   1 early() {
       sim.addSubpop("p1", K * W * W);
       p1.setSpatialBounds(c(0.0, 0.0, W, W));
       for (ind in p1.individuals) {
           ind.setSpatialPosition(p1.pointUniform());
       }
   }

   early() { // survival probabilities
       p1.fitnessScaling = 1 - MU;
   }

   late() {
       i1.evaluate();
   }

   500 early() {
       sim.treeSeqRememberIndividuals(p1.individuals);
   }

   1000 {
       sim.treeSeqOutput("spatial_sim.trees");
       catn("Done.");
       sim.simulationFinished();
   }


Ok, now let's have a quick look at the output:

.. code-block:: python

   import pyslim, tskit
   import numpy as np

   slim_ts = pyslim.load("spatial_sim.trees")
   print(f"Tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
         f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
         f" and {slim_ts.num_mutations} mutations.")

Running this code, we get

.. code-block:: none

   Tree sequence has 7139 trees along the genome of length 100000000.0, 1605 individuals,
   2410 'sample' genomes, and 0 mutations.

It makes sense we have no mutations: we haven't added any yet.
The tree sequence is recording the relationship between 2,410 genomes (the "samples"),
which requires 7,139 distinct trees along the genome.
Individuals are diploid, so why is the number of individuals not equal to twice the number of samples? 
Recall that for the *next* step, recapitation, it is necessary to keep the genomes from the first
generation around in the tree sequence so we can trace back lineages from them.
These extra individuals are from the first generation.
Let's count up when the individuals in the tree sequence were born:

.. code-block:: python

   for t in np.unique(slim_ts.individual_times):
     print(f"There are {np.sum(slim_ts.individual_times == t)} individuals from time {t}.")
   
This gets us

.. code-block:: none

   There are 327 individuals from time 0.0.
   There are 162 individuals from time 1.0.
   There are 82 individuals from time 2.0.
   There are 25 individuals from time 3.0.
   There are 20 individuals from time 4.0.
   There are 13 individuals from time 5.0.
   There are 5 individuals from time 6.0.
   There are 4 individuals from time 7.0.
   There are 1 individuals from time 8.0.
   There are 283 individuals from time 500.0.
   There are 143 individuals from time 501.0.
   There are 66 individuals from time 502.0.
   There are 28 individuals from time 503.0.
   There are 27 individuals from time 504.0.
   There are 9 individuals from time 505.0.
   There are 6 individuals from time 506.0.
   There are 2 individuals from time 507.0.
   There are 2 individuals from time 508.0.
   There are 400 individuals from time 999.0.

These are now *tskit* times, which, confusingly,
is measured in units of "time since the end of the simulation".
So, this tells us that 400 of the individuals are those that we initialized the simulation with
(we ran it for 1000 time steps and it began at SLiM time step 1, which is 999 time steps ago),
and there's some more individuals around 500 time steps ago and some more in the past few time steps.
This is a non-Wright-Fisher simulation,
and so individuals may live for more than one time step.
Let's check that all these individuals are alive at either (a) today, (b) 500 time steps ago,
or (c) the start of the simulation.

.. code-block:: python

   for t in [0, 500, 999]:
      alive = slim_ts.individuals_alive_at(t)
      print(f"There were {len(alive)} individuals alive {t} time steps in the past.")


This tells us that

.. code-block:: none

   There were 639 individuals alive 0 time steps in the past.
   There were 566 individuals alive 500 time steps in the past.
   There were 400 individuals alive 999 time steps in the past.

And, 639 + 566 + 400 is 1605, the total number of individuals.
So, this all checks out.


*************************
Recapitation and mutation
*************************

Next, we want to (a) simulate some ancestral diversity and (b) add in neutral mutations.
Please see `Haller et al (2019) <https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12968>`_
for the why and how of these steps.
But, first let's see if recapitation is necessary:
on how much of the genome is the tree sequence not coalesced?
In other words, recapitation adds diversity present in the initial generation;
will it make a difference? 
In fact, there is *no* segment of the genome that has coalesced;
*every* tree has more than one root:

.. code-block:: none

   >> sum([t.num_roots == 1 for t in slim_ts.trees()])
   0

Next, we will:

1. Recapitate, running a coalescent simulation to build ancestral trees.
2. Mutate, adding neutral variation.

We *won't* simplify at this point, although it would not hurt,
and if we did it would have to come *after* these steps.

.. code-block:: python

   import msprime

   recap_ts = slim_ts.recapitate(recombination_rate=1e-8, Ne=2000)
   ts = pyslim.SlimTreeSequence(
         msprime.mutate(recap_ts, rate=1e-8, keep=True))

   print(f"The tree sequence now has {ts.num_trees} trees,"
         f" and {ts.num_mutations} mutations.")

Note since :meth:`mutate <msprime.mutate>` is an msprime method, it does not return a pyslim
tree sequence, so we need to convert it back.
This has added mutations according to an infinite-sites model of mutation,
resulting in

.. code-block:: none

   The tree sequence now has 16571 trees,and 31125 mutations.

Two notes:
We will have no further use for `slim_ts` or for `recap_ts`;
we've just given them separate names for tidiness.
And, since the original SLiM mutation had no mutations, we didn't need to specify `keep=True`
in :meth:`mutate <msprime.mutate>`, but if we *had* put down selected mutations with SLiM
we'd probably want to keep them around.


****************************
Take a sample of individuals
****************************

Ok, it's time to compute some things.
In real life we don't get to work with *everyone* usually,
so we'll take a subset of individuals.
The range we have simulated has width and height 20 units,
with a population density of around 1 per unit area.
We'll get genomes to work with by pulling out

1. All the modern individuals in the squares of width 4 in the corners of the range, and
2. Five individuals sampled randomly from everyone alive 500 time steps ago.

.. code-block:: python


   np.random.seed(23)

   alive = ts.individuals_alive_at(0)
   locs = ts.individual_locations[alive, :]

   cwidth = 3
   corners = {
      'topleft' : alive[np.logical_and(locs[:, 0] < cwidth, locs[:, 1] < 5)],
      'topright' : alive[np.logical_and(locs[:, 0] < cwidth, locs[:, 1] > 15)],
      'bottomleft' : alive[np.logical_and(locs[:, 0] > 20 - cwidth, locs[:, 1] < 5)],
      'bottomright' : alive[np.logical_and(locs[:, 0] > 20 - cwidth, locs[:, 1] > 15)],
      'center' : alive[np.logical_and(np.abs(locs[:, 0] - 10) < cwidth,
                                      np.abs(locs[:, 1] - 10) < cwidth)]
      }

   old_ones = ts.individuals_alive_at(500)
   corners['ancient'] = np.random.choice(old_ones, size=5)

   for k in corners:
      print(f"We have {len(corners[k])} individuals in the {k} group.")


.. code-block:: none

   We have 51 individuals in the topleft group.
   We have 16 individuals in the topright group.
   We have 22 individuals in the bottomleft group.
   We have 35 individuals in the bottomright group.
   We have 46 individuals in the center group.
   We have 5 individuals in the ancient group.


To keep names associated with each subset of individuals,
we've kept the individuals in a dict, so that for instance
`corners["topleft"]` is an array of all the individual IDs that are in the top left corner.
The IDs of the ancient individuals we will work with are kept in the array `ancient`.


******************
Plotting locations
******************

We should check this: plot where these individuals lie
relative to everyone else.
The individuals locations are available in individual metadata,
but to make things easier, it's also present in a `num_individuals x 2`
numpy array as `ts.individual_locations`.
Since `corners["topleft"]` is an array of individual IDs,
we can pull out the locations of the "topleft" individuals
by indexing the rows of the individual location array:

.. code-block:: none

   >>> ts.individual_locations
   array([[ 7.89486192, 19.10248491,  0.        ],
          [ 6.70972958,  4.79636467,  0.        ],
          [ 5.70948388,  3.93691031,  0.        ],
          ...,
          [19.58402983,  1.70062556,  0.        ],
          [ 8.74254249,  7.87962527,  0.        ],
          [ 8.48533111,  6.63358757,  0.        ]])
   >>> ts.individual_locations.shape
   (1605, 3)
   >>> ts.individual_locations[corners["topleft"], :].shape
   (51, 3)

Note that the locations arrays have *three* columns because SLiM allows for
`(x, y, z)` coordinates. We'll just use the first two columns, for `x` and `y`.

Using this, we can easily plot the locations of all the individuals from today
(on the left) and 500 time steps ago (on the right).
We have to do a bit of mucking around to set the colors so that they reflect
which group each individual is in.

.. code-block:: python

   import matplotlib
   matplotlib.use('Agg')
   import matplotlib.pyplot as plt

   group_order = ['topleft', 'topright', 'bottomleft', 'bottomright', 'center', 'ancient']
   ind_colors = np.repeat(0, ts.num_individuals)
   for j, k in enumerate(group_order):
      ind_colors[corners[k]] = 1 + j

   old_locs = ts.individual_locations[old_ones, :]

   fig = plt.figure(figsize=(12, 6))
   ax = fig.add_subplot(121)
   ax.set_title("today")
   ax.scatter(locs[:,0], locs[:,1], s=10, c=ind_colors[alive])
   ax = fig.add_subplot(122)
   ax.set_title("long ago")
   ax.scatter(old_locs[:, 0], old_locs[:, 1], s=10, c=ind_colors[old_ones])
   fig.savefig("spatial_sim_locations.png")


.. image:: _static/spatial_sim_locations.png
   :width: 1200px
   :alt: Spatial location of all individuals and the genotyped ones.



*********************
Isolation by distance
*********************

Now, let's look at *isolation by distance*, i.e.,
let's compare geographic and genetic distances.
Here, "genetic distance" will be mean pairwise sequence divergence.
First, we'll compute mean genetic distance between each of our five groups.

The first thing we need to do is some bookkeeping.
So far, we've just worked with *individuals*,
but tree sequence tools, in particular the statistics computation methods from tskit,
are designed to work with *genomes*, also known as "nodes".
So, first we need to pull out the *node IDs* corresponding to the individuals we want.
The things that make up a tree sequence - individuals, nodes, mutations, etcetera -
can generally be examined individually. 
For instance, here's what we have for the five "ancient" individuals:

.. code-block:: none

   >>> for i in ancient:
   ...   print(ts.individual(i))
   ... 
   {'id': 1079, 'flags': 131072, 'location': array([17.80637122, 10.01626   ,  0.        ]),
      'metadata': IndividualMetadata(pedigree_id=136192, age=4, population=1, sex=-1, flags=0),
      'nodes': array([880, 881], dtype=int32), 'population': 1, 'time': 504.0}
   {'id': 1527, 'flags': 131072, 'location': array([7.36765899, 1.98245938, 0.        ]),
      'metadata': IndividualMetadata(pedigree_id=137290, age=0, population=1, sex=-1, flags=0),
      'nodes': array([1776, 1777], dtype=int32), 'population': 1, 'time': 500.0}
   {'id': 1070, 'flags': 131072, 'location': array([18.97646531,  0.86655059,  0.        ]),
      'metadata': IndividualMetadata(pedigree_id=136095, age=4, population=1, sex=-1, flags=0),
      'nodes': array([862, 863], dtype=int32), 'population': 1, 'time': 504.0}
   {'id': 1276, 'flags': 131072, 'location': array([15.02692262, 14.84464791,  0.        ]),
      'metadata': IndividualMetadata(pedigree_id=137001, age=1, population=1, sex=-1, flags=0),
      'nodes': array([1274, 1275], dtype=int32), 'population': 1, 'time': 501.0}
   {'id': 1499, 'flags': 131072, 'location': array([ 3.52153941, 11.4547264 ,  0.        ]),
      'metadata': IndividualMetadata(pedigree_id=137262, age=0, population=1, sex=-1, flags=0),
      'nodes': array([1720, 1721], dtype=int32), 'population': 1, 'time': 500.0}

Notice that among other things, each `individual` carries around a list of their `node` IDs,
i.e., their genomes.
We need to put these all in a list of lists,
so that, for instance, the first element of the list will have the node IDs of all the genomes
of the individuals in the "topleft" group.
And, since we kept the individual IDs in a dict, which are unordered,
we'll have to do some extra work to make sure we keep track of order.

.. code-block:: python


   sampled_nodes = [[] for _ in corners]
   for j, k in enumerate(group_order):
      for ind in corners[k]:
         sampled_nodes[j].extend(ts.individual(ind).nodes)


Let's do a sanity check: the number of nodes in each element of this list
should be twice the number of individuals in the corresponding list.

.. code-block:: none

   >>> [len(corners[k]) for k in corners]
   [51, 16, 22, 35, 46, 5]

   >>> [len(u) for u in sampled_nodes]
   [102, 32, 44, 70, 92, 10]
   
So, in the 'topleft' corner there are 51 diploids. That checks out.   

Now, we can compute the matrix of pairwise mean sequence divergences
between and within these sets.
This is done using the :meth:`ts.divergence <tskit.TreeSequence.divergence>` method.

.. code-block:: python

   pairs = [(i, j) for i in range(6) for j in range(6)]
   group_div = ts.divergence(sampled_nodes, indexes=pairs)[0].reshape((6, 6))

   print("\t" + "\t".join(group_names))
   for i, group in enumerate(group_order):
      print(f"{group_order[i]}:\t" + "\t".join(map(str, np.round(group_div[i], 7))))


.. code-block:: none

   topleft:	1.41e-05	5.68e-05	4.69e-05	4.67e-05	4.42e-05	6.05e-05
   topright:	5.68e-05	1e-07	5.48e-05	5.53e-05	5.67e-05	6.28e-05
   bottomleft:	4.69e-05	5.48e-05	1.91e-05	2.93e-05	3.35e-05	5.88e-05
   bottomright:	4.67e-05	5.53e-05	2.93e-05	1.5e-05	2.68e-05	5.99e-05
   center:	4.42e-05	5.67e-05	3.35e-05	2.68e-05	2.98e-05	6.08e-05
   ancient:	6.05e-05	6.28e-05	5.88e-05	5.99e-05	6.08e-05	5.34e-05

That's nice, but to look at isolation by distance,
we should actually separate out the individuals.
To do that, we need to create a list of lists of nodes
whose j-th entry is the nodes belonging to the j-th individual,
and to keep track of which group each one belongs to.

.. code-block:: python

   ind_nodes = []
   ind_group = []
   for j, group in enumerate(group_order):
      for ind in corners[group]:
         ind_nodes.append(ts.individual(ind).nodes)
         ind_group.append(group_order[j])

   nind = len(ind_nodes)
   pairs = [(i, j) for i in range(nind) for j in range(nind) if i <= j]
   ind_div = ts.divergence(ind_nodes, indexes=pairs)[0]

Here we've only computed divergences in the *upper triangle* of the pairwise divergence matrix,
with heterozygosities on the diagonal.
We'll also need pairwise geographic distances:

.. code-block:: python

   geog_dist = np.repeat(0.0, len(pairs))
   locs = ts.individual_locations
   for k, (i, j) in enumerate(pairs):
      geog_dist[k] = np.sqrt(np.sum((locs[i, :] - locs[j, :])**2))


Let's check that makes sense: distances of individuals from themselves should be zero.

.. code-block:: python

   for (i, j), x in zip(pairs, geog_dist):
     if i == j:
       assert(x == 0)

Python does not complain, which is good.
Now let's plot genetic distance against geographic distance.

.. code-block:: python

   pair_colors = np.repeat(0, len(pairs))
   for k, (i, j) in enumerate(pairs):
      if ind_group[i] == "ancient" or ind_group[j] == "ancient":
         pair_colors[k] = 1

   fig = plt.figure(figsize=(6, 6))
   ax = fig.add_subplot(111)
   ax.scatter(geog_dist, 1e3 * ind_div, s=20, alpha=0.5,
              c=pair_colors)
   ax.set_xlabel("geographic distance")
   ax.set_ylabel("genetic distance (diffs/Kb)")
   fig.savefig("spatial_sim_ibd.png")


.. image:: _static/spatial_sim_ibd.png
   :width: 600px
   :alt: Geographic and genetic distances in the simulation.


Since we multiplied `ind_div` by 1,000,
the units of genetic distance are in mean number of nucleotide differences per kilobase.
There is *not* strong IBD in this noisy and relatively small simulation,
but notice that the "ancient" samples are more deeply diverged from modern samples (in yellow)
than most modern ones are from each other.



*****************************
Simplification and VCF output
*****************************

Now we want to write out these data for analysis with other programs.
Recall that the tree sequence contains information about *everyone* in the population.
Before we output, we will remove this extra information.
We could have done this earlier (e.g., before statistic computation)
but since simplification reorders nodes and individuals,
we would have added an extra layer of bookkeeping.
However, it's necessary here because at this point,
:meth:`ts.write_vcf <tskit.TreeSequence.write_vcf>` only outputs *all* sample genomes,
so we need to remove the other ones,
and make sure they are in the order we want.

To do this, and make sure that everything stays nicely cross-referenced,
we're going to loop through our individuals, writing their information to a file,
while at the same time constructing a list of node IDs, in the same order.
We'll use this as the `samples` argument to `simplify`,
because in the simplified tree sequence the sample nodes are these, in order;
and in the VCF output each two samples, in order, are represented as a diploid individual.
In the resulting VCF file, the individuals will be called `msp_0`, `msp_1`, etcetera;
so we'll record that in the individual file also.

.. code-block:: python

   nodelist = []
   with open("spatial_sim_individuals.txt", "w") as indfile:
      indfile.writelines("\t".join(["vcf_label", "slim_id", "birth_time_ago", "age", "x", "y"]) + "\n")
      j = 0
      for group in group_order:
         for i in corners[group]:
            vcf_label = f"msp_{j}"
            ind = ts.individual(i)
            nodelist.extend(ind.nodes)
            data = [vcf_label, str(ind.id), str(ind.time), str(ind.metadata.age),
                    str(ind.location[0]), str(ind.location[1])]
            indfile.writelines("\t".join(data) + "\n")
            j += 1

   sts = ts.simplify(nodelist)
   with open("spatial_sim_genotypes.vcf", "w") as vcffile:
      sts.write_vcf(vcffile, ploidy=2)


****************
More information
****************

1. We plan to use individual information when writing out to VCF,
   but `this is not yet implemented <https://github.com/tskit-dev/tskit/issues/73>`_.

2. The distinction between "nodes" (i.e., genomes) and "individuals" can be confusing,
   as well as the idea of "samples".
   Please see the
   `tskit documentation <https://tskit.readthedocs.io/en/latest/data-model.html>`_
   for more explanation about these concepts.

3. The general interface for computing statistics (explaining, for instance, the "indexes"
   argument above) is described in
   `the tskit documentation <https://tskit.readthedocs.io/en/latest/stats.html>`_
   also.


