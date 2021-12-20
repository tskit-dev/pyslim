***************************
[UPCOMING.X.X] - XXXX-XX-XX
***************************

********************
[0.700] - 2021-02-24
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

