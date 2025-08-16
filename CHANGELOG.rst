***************************
[1.1.0] - 2025-08-16
***************************

Major update release to support the release of SLiM v5.0. The main update here
is support for multiple chromosomes, which entailed a change to the
node table metadata that records "vacancy" for chromosomes not diploid
in a given individual. Changes are minimal if you're simulating a single,
diploid chromosome; for more information see
https://tskit.dev/pyslim/docs/latest/previous_versions.html

**Breaking changes**:

- SLiM tree sequence file version number has changed to 0.9.

- This is a change in SLiM, really, but top-level SLiM metadata now requires
  a `"this_chromosome"` entry.

- Similarly, node metadata no longer has `genome_type` or `is_null`; instead
  they have `is_vacant`, and the chromosome type is in top-level metadata,
  under `metadata['SLiM']['this_chromosome']`.

- `pyslim.slim_metadata_schemas["node"]` has an appropriate metadata schema
    for a single-chromosome simulation (not a breaking change), but users
    wishing to set up a multi-chromosome simulation should use instead
    `pyslim.slim_node_metadata_schema` (which has the appropriate value in
    `["properties"]["is_vacant"]["length"]`)
    (:user:`petrelharp`, :pr:`367`).

- Previously, `pyslim.annotate` would leave existing node and individual
    metadata, even if this metadata came from a different schema. This
    could silently create garbage metadata. Now, `annotate` removes
    any existing metadata, and warns if this occurs
    (:user:`petrelharp`, :pr:`390`).

**Notable changes**:

- `pyslim.individual_parents` and `pyslim.has_individual_parents` are no
    longer needed, as this information can be more reliably extracted from
    the `parents` column/attribute of individuals, so they are now deprecated
    and will be removed in the future.

**New features**:

- `pyslim.set_slim_state` will adjust times and "alive" flags so
    that when the tree sequence is loaded into SLiM it will have
    a specified set of individuals alive at a particular time.
    (:user:`petrelharp`, :pr:`384`)

- Functions `pyslim.node_is_vacant` and `pyslim.has_vacant_samples`
    test for vacancy in the current chromosome.
    (:user:`petrelharp`, :pr:`367`)

- `pyslim.remove_vacant` and `pyslim.restore_vacant`, respectively,
    remove and restore vacant samples nodes, necessary for recapitation
    and other operations. Also, corresponding `pyslim.remove_vacant_tables`
    and `pyslim.restore_vacant_tables`.
    (:user:`petrelharp`, :pr:`367`)

- `pyslim.recapitate` by default removes the sample flags on vacant nodes,
    but this behavior is controlled by an argument, `keep_vacant`.
    (:user:`petrelharp`, :pr:`367`)

- `pyslim.set_metadata_schemas` now includes a `num_chromosomes` argument
    (:user:`petrelharp`, :pr:`367`)


**Bugfixes**:

- Recapitation on tree sequences with null genomes would attempt to simulate
    the history of those null genomes; this would in all but exceptional cases
    fail with an error ("not all roots are at the time expected"). Now, null
    genomes are "vacant" (see above) and `recapitate` removes their
    sample flags before recapitating (and optionally puts them back)
    as described in `pyslim.remove_vacant` (:user:`petrelharp`, :pr:`367`).

- Previously, recapitation would require the roots of all trees to be
    at the same time (roughly) as the 'tick' stored in the top-level metadata;
    however, this would not be the case if the first population was added
    later than the first tick. The requirement has therefore been removed.
    (:user:`petrelharp`, :pr:`382`)

- The `generated_nucleotides` method now sets the `nucleotide_based` entry
    in top-level metadata to True. (:user:`petrelharp`, :pr:`385`)

- The individual flags `INDIVIDUAL_ALIVE`, `INDIVIDUAL_REMEMBERED`,
    and `INDIVIDUAL_RETAINED` were signed integers, but the flags in the
    individual table they apply to are unsigned, so using the
    bitwise negation operator `~` could result in an error. Now,
    they are np.uint32 values. (:user:`petrelharp`, :pr:`378`)


***************************
[1.0.4] - 2023-08-01
***************************

**Bugfixes**:

- The last bugfix introduced a small bug: recapitation on a tree sequence
    whose roots are at least 100,000 ticks ago would produce an msprime error:
    "Attempt to sample a lineage from an inactive population". Reported by
    Meaghan Clark. (:user:`petrelharp`, :pr:`322`)

***************************
[1.0.3] - 2023-06-21
***************************

**Bugfixes**:

- From 1.0.1 back to 0.700, there was a bug in `recapitate` when using the
    `ancestral_Ne` parameter that introduced a bottleneck to diploid size Ne=1
    for each SLiM subpopulation for 1 or 2 generations *unless* either (a) it
    was a WF simulation, with calls to addSubPop() in first() or early() and
    treeSeqOutput() in late(), or (b) it was a nonWF simulation, with calls to
    addSubPop() in first() and treeSeqOutput() in early() or late(). The fix
    correctly starts the msprime population with effective size `ancestral_Ne`
    at the time of the roots, which might be at the value of
    `ts.metadata['SLiM']['tick']`, this value minus 1, or this value minus 2.
    Furthermore, `recapitate` now throws an error if any roots of any trees
    are not at the same time as the others. (:user:`petrelharp`, :pr:`308`)


***************************
[1.0.2] - 2023-06-20
***************************

This was a bugfix release that was pushed out without the actual bug fix.
Please don't use this one.

***************************
[1.0.1] - 2022-09-23
***************************

- Documentation of how to empirically measure generation time
    and check that it is correct
    (:user:`silastittes`, :user:`petrelharp`, :pr:`301`, :pr:`293`).

- Minor modifications to `convert_alleles` and `generate_nucleotides`
    so that they run in a reasonable amount of time
    (:user:`petrelharp`, :pr:`299`).

- Addition of method to find the next SLiM mutation ID,
    `pyslim.next_slim_id` (:user:`mufernando`, :pr:`290`).


***************************
[1.0] - 2022-08-12
***************************

**Breaking changes**:

- Removed `SlimTreeSequence` class entirely (it was previously deprecated).
    All its methods are either available in `tskit.TreeSequence`
    or are now called by `pyslim.fn(ts, ...)` instead of `ts.fn(...)`.

- TODO: Deprecated `util.unique_labels_by_group`.

- Moved some methods of `SlimTreeSequence` to pyslim:
    * instead of `slim_ts.slim_time(t)` do `pyslim.slim_time(ts, t)`
    * instead of `slim_ts.individuals_alive_at(t)` do `pyslim.individuals_alive_at(ts, t)`
    * instead of `slim_ts.individuals_parents(t)` do `pyslim.individuals_parents(ts, t)`
    * instead of `slim_ts.individuals_ages(t)` do `pyslim.individuals_ages(ts, t)`

- The methods `slim_ts.mutation_at( )` and `slim_ts.nucleotide_at( )`
    are now methods of pyslim, whose first argument is the tree sequence.

- In SLiM v4 "generation" has been renamed to "tick", and so corresponding things
  in pyslim have been renamed: top-level metadata now has `ts.metadata["SLiM"]["tick"]`
  instead of `ts.metadata["SLiM"]["generation"]`

- Renamed `pyslim.annotate_defaults()` to `pyslim.annotate()`, with slight
  changes in behavior: since msprime.sim_ancestry() now simulates individuals
  by default, annotation does not set up individuals: if you have a tree
  sequence without individuals (e.g., produced by msprime.simulate()) then you
  need to set up those individuals yourself.

- To update a tree sequence produced by an old version of SLiM to the current one,
  use `pyslim.update( )`. (However, note that reading it in to SLiM and
  writing it out again might be even easier.)

- The method `pyslim.set_tree_sequence_metadata` now has arguments `tick` and `cycle`
  instead of `generation`.

- Removed `pyslim.make_slim_provenance_dict`.

**Other notable changes**:

- Top-level metadata now has a `tick` attribute that is (for now) a synonym
    for `generation`; the latter will be deprecated at some point in the future.

- Methods for getting time, population, and location information about individuals
  are now in tskit:
    * `SlimTreeSequence.individual_times` is now `TreeSequence.individuals_time()`
    * `SlimTreeSequence.individual_populations` is now `TreeSequence.individuals_population()`
    * `SlimTreeSequence.individual_locations` is now `TreeSequence.individuals_location()`
  However, this will be invisible to the user. In each case note the the
  location of the "s" has moved (to "individual*s* time" instead of "individual
  time*s*"), but the original version remains an undocumented alias.

**New features**:

- Methods like `pyslim.individuals_alive_at( )` now deal with the new `stage="first"`.


********************
[0.700] - 2021-12-20
********************

**Breaking changes**:

- `pyslim.recapitate` is updated to use new demography features in msprime 1.0,
    and differs from `SlimTreeSequence.recapitate()` (now deprecated). Since
    the backend is now `msprime.sim_ancestry()` instead of `msprime.simulate()`,
    the argument `Ne` should be replaced with `ancestral_Ne`.

- `reference_sequence` is now a tskit TreeSequence attribute, no longer managed
    by pyslim. It is no longer mutable on tree sequences (only TableCollections),
    and previous calls to `ts.reference_sequence` to get the actual sequence
    should be replaced by `ts.reference_sequence.data`.

- Old-style "legacy" metadata (previously deprecated) has been removed.
    See `the documentation <https://tskit.dev/pyslim/docs/previous_versions.html>`_
    for instructions on migrating your code.


**New features**:

- Added `pyslim.population_size( )` to compute an array giving numbers of
    individuals across a grid of space and time bins. ({user}giliapatterson)


********************
[0.600] - 2021-02-24
********************

**New features**:

- Added `ts.individual_parents()`, a way to get the IDs of individual's parents
    when both of them are present in the tree sequence. :user:@petrelharp

- Added and documented `TSK_INDIVIDUAL_RETAINED` flag to reflect the additional
    of "retained" individuals in SLiM v3.6. :user:@hyanwong, :user:@petrelharp

**Bugfix**:

- Modified `recaptiate` to not error with the current msprime 1.0 alpha release.

********************
[0.501] - 2020-12-08
********************

**Bugfix**:

- Making `.slim_generation` derive from the tree sequence's top-level metadata
    had the unanticipated consequence that it could not be modified, which some
    people were doing. This restores the previous behavior, but in the future,
    modifying `.slim_generation` on a tree sequence will be deprecated - instead,
    this should be modified in the metadata of the TableCollection.

********************
[0.500] - 2020-12-07
********************

**Breaking changes**:

- "First generation" individuals no longer need to be retained by SLiM to recapitate,
  thanks to the "keep_input_roots" argument to simplify (new in tskit 0.3.0).
  The FIRST_GEN flag and `.first_generation_individuals()` methods are now deprecated,
  and if you want these to remain in the tree sequence you must explicitly Remember them.
  (However, their *nodes* will remain if necessary for recapitation.)
  If you wish to simplify an un-recapitated tree sequence you now can, but you must
  pass `keep_input_roots=True`. This should only cause breakages if you made explicit
  use of the first generation individuals, without explicitly Remembering them.

- Information about the tree sequence is now stored in *top-level metadata*,
  accessible through `ts.metadata['SLiM']`. Previous interfaces remain: for instance,
  `ts.slim_generation` is now redundant with `ts.metadata['SLiM']['generation']`.
  This should not cause breakages, but will cause warnings where none were previously:
  for instance, `pyslim.SlimTreeSequence(msprime.mutate(ts))` may throw a warning
  because `msprime.mutate( )` does not preserve top-level metadata, and so SLiM-relevant
  information is retrieved from provenance (as in previous file versions).

**Notable changes**:

- Switched to using tskit native encoding/decoding of metadata via schemas.
- added to conda-forge (@winni2k)

**New features**:

- added `samples_only` and `population` arguments to `ts.individuals_alive_at()`
- added the `ts.slim_time()` method
- enabled dumping the reference sequence for nucleotide models

********************
[0.403] - 2020-08-27
********************

BUGFIX: if a tree had all first generation individuals removed
   (e.g., if it had been simplified) then individuals_alive_at( ) failed.

********************
[0.402] - 2020-08-27
********************


This is a compatibility release, for the tskit 0.3.0 release.


**New features**:

- added has_individual_parents, a method to find individuals with all parents
  are also recorded as individuals
- Provenance handling:
   * added the `.slim_provenances` property to return all SLiM provenance entries
   * added the `slim_provenance_version` and `parse_provenance` methods to tell if
      provenance entries come from SLiM and to parse them

- documentation for recapitation with a nonuniform map by :user:@TeresaPegan

**Bug fixes**:

- fixed differential time offset for tree sequences saved out in early versus late:
   prior to this, mutation_at and nucleotides_at would have been sometimes wrong if the tree sequence
   was saved out during late

- initialises correctly to work with tskit 0.3.0

********************
[0.401] - 2020-03-27
********************

**Bug fixes**:

- checks for the ability to simulate with a discrete recombination map
   in the available version of msprime, and sets the default flat
   recombination map in recapitate appropriately

********************
[0.400] - 2020-03-24
********************

**New features**:

- updated to take and output SLiM file version 0.4, which only differs from 0.3
   in minor aspects of provenance

********************
[0.314] - 2019-10-31
********************

**New features**:

- allows passing in of a recombination map to recapitate (:user:`mufernando`)
- added first_generation_individuals() function
- defined individual ages for WF ages
- added mutation_at() and fixed up nucleotide_at() functions

