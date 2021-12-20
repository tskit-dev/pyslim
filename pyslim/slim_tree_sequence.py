import attr
import struct
import msprime
import tskit
import kastore
import json
from collections import OrderedDict
import warnings
import numpy as np

from ._version import *
from .slim_metadata import *
from .provenance import *
from .util import *
from .slim_metadata import _old_metadata_schema

INDIVIDUAL_ALIVE = 2**16
INDIVIDUAL_REMEMBERED = 2**17
INDIVIDUAL_RETAINED = 2**18
# deprecated but keep it around for backwards compatibility
# (also, it means effectively the same thing as RETAINED)
INDIVIDUAL_FIRST_GEN = INDIVIDUAL_RETAINED

# A nucleotide k in mutation metadata actually means
# something that in reference_sequence is NUCLEOTIDES[k]
NUCLEOTIDES = ['A', 'C', 'G', 'T']


def _deprecation_warning(w):
    warnings.warn(
        "The SlimTreeSequence class is being phased out, "
        "as most important functionality is provided by tskit. "
        "Please see the "
        "`documentation <https://tskit.dev/pyslim/latest/previous_versions.html>`_. "
        f"{w}",
        FutureWarning
    )


def load(path, legacy_metadata=False):
    '''
    Load the SLiM-compatible tree sequence found in the .trees file at ``path``.

    :param string path: The path to a .trees file.
    '''
    if legacy_metadata:
        raise ValueError(
            "It looks like you're trying to use pyslim legacy metadata tools, "
            "which are no longer supported. See `the documentation "
            "<https://tskit.dev/pyslim/pyslim/latest/metadata.html#legacy-metadata>`_"
            "for how to update your script."
    )

    ts = SlimTreeSequence.load(path)
    return ts


def load_tables(tables, **kwargs):
    '''
    See :func:`SlimTreeSequence.load_tables`.

    :param TableCollection tables: A set of tables.
    '''
    ts = SlimTreeSequence.load_tables(tables, **kwargs)
    return ts


def annotate_defaults(ts, **kwargs):
    '''
    Takes a tree sequence (as produced by msprime, for instance), and adds in the
    information necessary for SLiM to use it as an initial state, filling in
    mostly default values. Returns a :class:`SlimTreeSequence`.

    :param TreeSequence ts: A :class:`TreeSequence`.
    :param string model_type: SLiM model type: either "WF" or "nonWF".
    :param int slim_generation: What generation number in SLiM correponds to
        ``time=0`` in the tree sequence.
    :param str reference_sequence: A reference sequence of length
        equal to ts.sequence_length.
    :param bool annotate_mutations: Whether to replace mutation metadata
        with defaults. (If False, the mutation table is unchanged.)
    '''
    tables = ts.dump_tables()
    annotate_defaults_tables(tables, **kwargs)
    return SlimTreeSequence.load_tables(tables)


def annotate_defaults_tables(tables, model_type, slim_generation, reference_sequence=None,
        annotate_mutations=True):
    '''
    Does the work of :func:`annotate_defaults()`, but modifies the tables in place: so,
    takes tables as produced by ``msprime``, and makes them look like the
    tables as output by SLiM. See :func:`annotate_defaults` for details.
    '''
    if (type(slim_generation) is not int) or (slim_generation < 1):
        raise ValueError("SLiM generation must be an integer and at least 1.")
    # set_nodes must come before set_populations
    if model_type == "WF":
        default_ages = -1
    elif model_type == "nonWF":
        default_ages = 0
    else:
        raise ValueError("Model type must be 'WF' or 'nonWF'")
    if not np.allclose(tables.sites.position, np.floor(tables.sites.position)):
        raise ValueError(
                "Site positions in this tree sequence are not at integer values, "
                "but must be for loading into SLiM: generate mutations with "
                "sim_mutations(..., discrete_genome=True), not simulate()."
        )
    top_metadata = default_slim_metadata('tree_sequence')['SLiM']
    top_metadata['model_type'] = model_type
    top_metadata['generation'] = slim_generation
    set_tree_sequence_metadata(tables, **top_metadata)
    set_metadata_schemas(tables)
    _set_nodes_individuals(tables, age=default_ages)
    _set_populations(tables)
    if annotate_mutations:
        _set_sites_mutations(tables)
    if reference_sequence is not None:
        tables.reference_sequence.data = reference_sequence


class MetadataDictWrapper(dict):
    '''
    A simple wrapper around metadata dicts that will throw an informative error
    message if ``md.X`` is used instead of ``md["X"]``, and (for mutation metadata)
    if ``md[k]`` is used instead of ``md["mutation_list"][k]``.
    '''

    def __getattr__(self, name):
        if name in self.keys():
            raise AttributeError(
                f"'dict' object has no attribute '{name}'. "
                "It looks like you're trying to use the legacy "
                "metadata interface: see the "
                "`documentation <https://tskit.dev/pyslim/latest/previous_versions.html#legacy-metadata>`_ "
                "for how to switch over your script"
            )
        else:
            raise AttributeError(f"'dict' object has no attribute '{name}'")

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError as e:
            if isinstance(key, int):
                msg = e.args[0]
                e.args = (
                    f"{msg}: It looks like you're trying to use the legacy "
                    "metadata interface: see the "
                    "`documentation <https://tskit.dev/pyslim/latest/previous_versions.html#legacy-metadata>`_ "
                    "for how to switch over your script",
                )
            raise e


class SlimTreeSequence(tskit.TreeSequence):
    '''
    This is just like a :class:`tskit.TreeSequence`, with a few more properties
    and methods, including:

    - :attr:`.slim_generation` - the SLiM "generation" at the end of the simulation
    - :attr:`.reference_sequence` - if using a nucleotide model, the reference sequence
    - :attr:`.individual_locations` - numpy array of individual locations
    - :attr:`.individual_ages` - numpy array of individiual ages
    - :attr:`.individual_times` - numpy array of how long ago each individual was born
    - :attr:`.individual_populations` - numpy array of individual's populations

    All mutable properties of individuals (e.g., age) are as they were last recorded
    in the tree sequence: either the last time they were Remembered, or at the end
    of the simulation, if they are still alive then.

    However, `.individual_populations` gives individuals' birth populations;
    for their final location use the subpopulation attribute of metadata.

    You can create a :class:`.SlimTreeSequence` using one of

    - :meth:`.SlimTreeSequence.load_tables`, :meth:`.SlimTreeSequence.load`,
    - :func:`load`, or :func:`load_tables`.

    :ivar reference_sequence: None, or an string of length equal to the sequence
        length that gives the entire reference sequence for nucleotide models.
    '''

    def __init__(self, ts, legacy_metadata=False):
        if legacy_metadata:
            raise ValueError(
                "It looks like you're trying to use pyslim legacy metadata tools, "
                "which are no longer supported. See `the documentation "
                "<https://tskit.dev/pyslim/latest/metadata.html#legacy-metadata>`_"
                "for how to update your script."
        )

        if not (isinstance(ts.metadata, dict) and 'SLiM' in ts.metadata
                and ts.metadata['SLiM']['file_version'] in compatible_slim_file_versions):
            tables = ts.dump_tables()
            update_tables(tables)
            ts = tables.tree_sequence()
        super().__init__(ts._ll_tree_sequence)
        # need this for backwards compatibility
        self._slim_generation = ts.metadata['SLiM']['generation']
        # pre-extract individual metadata
        self.individual_locations = ts.tables.individuals.location
        self.individual_locations.shape = (int(len(self.individual_locations)/3), 3)
        self.individual_ages = np.zeros(ts.num_individuals, dtype='int')
        if self.metadata['SLiM']['model_type'] != "WF":
            self.individual_ages = np.fromiter(map(lambda ind: ind.metadata['age'], ts.individuals()), dtype='int64')

        self.individual_times = np.zeros(ts.num_individuals)
        self.individual_populations = np.repeat(np.int32(-1), ts.num_individuals)
        if not np.all(unique_labels_by_group(ts.tables.nodes.individual,
                                              ts.tables.nodes.population)):
            raise ValueError("Individual has nodes from more than one population.")
        if not np.all(unique_labels_by_group(ts.tables.nodes.individual,
                                              ts.tables.nodes.time)):
            raise ValueError("Individual has nodes from more than one time.")
        has_indiv = (ts.tables.nodes.individual >= 0)
        which_indiv = ts.tables.nodes.individual[has_indiv]
        # if we did not do the sanity check above then an individual with nodes in more than one pop
        # would get the pop of their last node in the list
        self.individual_populations[which_indiv] = ts.tables.nodes.population[has_indiv]
        self.individual_times[which_indiv] = ts.tables.nodes.time[has_indiv]

    def __getstate__(self):
        return {
            'tables': self.dump_tables(),
        }

    def __setstate__(self, state):
        self.__init__(
            state['tables'].tree_sequence(),
        )

    @property
    def slim_generation(self):
        _deprecation_warning("Please access ts.metadata['SLiM']['generation'] instead.")
        return self._slim_generation

    @slim_generation.setter
    def slim_generation(self, value):
        _deprecation_warning("Please access ts.metadata['SLiM']['generation'] instead.")
        self._slim_generation = value

    @property
    def model_type(self):
        _deprecation_warning("Please access ts.metadata['SLiM']['model_type'] instead.")
        return self.metadata['SLiM']['model_type']

    @classmethod
    def load(cls, path, legacy_metadata=False):
        '''
        Load a :class:`SlimTreeSequence` from a .trees file on disk.

        :param string path: The path to a .trees file.
        :rtype SlimTreeSequence:
        '''
        ts = tskit.load(path)
        return cls(ts, legacy_metadata=legacy_metadata)

    @classmethod
    def load_tables(cls, tables, **kwargs):
        '''
        Creates the :class:`SlimTreeSequence` defined by the tables.

        :param TableCollection tables: A set of tables, as produced by SLiM
            or by :func:`annotate_defaults`.
        :rtype SlimTreeSequence:
        '''
        # a roundabout way to copy the tables
        ts = tables.tree_sequence()
        return cls(ts, **kwargs)

    def dump(self, path, **kwargs):
        '''
        Dumps the tree sequence to the path specified. This is only a wrapper
        :func:`tskit.TreeSequence.dump`.

        :param str path: The file path to write the TreeSequence to.
        :param kwargs: Additional keyword args to pass to tskit.TreeSequence.dump
        '''
        # temporary until we remove support for setting slim_generation
        if self._slim_generation != self.metadata['SLiM']['generation']:
            tables = self.dump_tables()
            md = tables.metadata
            md['SLiM']['generation'] = self._slim_generation
            tables.metadata = md
            tables.dump(path, **kwargs)
        else:
            super(SlimTreeSequence, self).dump(path, **kwargs)

    def simplify(self, *args, **kwargs):
        '''
        This is a wrapper for :meth:`tskit.TreeSequence.simplify`.
        The only difference is that this method returns the
        derived class :class:`.SlimTreeSequence`.

        If you have not yet recapitated your SlimTreeSequence, you probably want to
        pass ``keep_input_roots=True``, so that recapitation is possible in the future.

        :rtype SlimTreeSequence:
        '''
        # temporary until we remove support for setting slim_generation
        if self._slim_generation != self.metadata['SLiM']['generation']:
            tables = self.dump_tables()
            md = tables.metadata
            md['SLiM']['generation'] = self._slim_generation
            tables.metadata = md
            # note we have to go to a tree sequence to get the map_nodes argument
            sts = tables.tree_sequence().simplify(*args, **kwargs)
        else:
            sts = super(SlimTreeSequence, self).simplify(*args, **kwargs)
        if (type(sts) == tuple):
            ret = (SlimTreeSequence(sts[0]), sts[1])
        else:
            ret = SlimTreeSequence(sts)

        return ret

    def population(self, id_):
        '''
        Returns the population whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.population`, but with additional attributes:
        `slim_id`, `selfing_fraction`, `female_cloning_fraction`,
        `male_cloning_fraction`, `sex_ratio`,
        `bounds_x0`, `bounds_x1`, `bounds_y0`, `bounds_y1`, `bounds_z0`, `bounds_z1`,
        and `migration_records`.
        These are all recorded by SLiM in the metadata.

        Note that SLiM populations are usually indexed starting from 1,
        but in tskit from zero, so there may be populations (e.g., with id_=0)
        that have no metadata and were not created by SLiM.

        :param int id_: The ID of the population (i.e., its index).
        '''
        pop = super(SlimTreeSequence, self).population(id_)
        if pop.metadata is not None:
            pop.metadata = MetadataDictWrapper(pop.metadata)
        return pop

    def individual(self, id_):
        '''
        Returns the individual whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.individual`, but with additional attributes:
        `time`, `pedigree_id`, `age`, `slim_population`, `sex`, and `slim_flags`.
        The `time` and `population` properties are extracted from the nodes,
        and an error will be thrown if the individual's nodes derive from
        more than one population or more than one time.

        :param int id_: The ID of the individual (i.e., its index).
        '''
        ind = super(SlimTreeSequence, self).individual(id_)
        ind.population = self.individual_populations[id_]
        ind.time = self.individual_times[id_]
        ind.metadata = MetadataDictWrapper(ind.metadata)
        return ind

    def node(self, id_):
        '''
        Returns the node whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.node`, but with additional attributes:
        `slim_id`, `is_null`, and `genome_type`.
        These are all recorded by SLiM in the metadata.

        :param int id_: The ID of the node (i.e., its index).
        '''
        node = super(SlimTreeSequence, self).node(id_)
        if node.metadata is not None:
            node.metadata = MetadataDictWrapper(node.metadata)
        return node

    def mutation(self, id_):
        '''
        Returns the mutation whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.mutation`, but with additional attributes:
        `mutation_type`, `selection_coeff`, `population`, `slim_time`,
        and `nucleotide`.
        These are all recorded by SLiM in the metadata.

        :param int id_: The ID of the mutation (i.e., its index).
        '''
        mut = super(SlimTreeSequence, self).mutation(id_)
        mut.metadata = MetadataDictWrapper(mut.metadata)
        return mut

    def recapitate(self,
                   recombination_rate=None,
                   population_configurations=None,
                   recombination_map=None, **kwargs):
        '''
        Returns a "recapitated" tree sequence, by using msprime to run a
        coalescent simulation from the "top" of this tree sequence, i.e.,
        allowing any uncoalesced lineages to coalesce.

        .. warning::

            This method is deprecated: please use pyslim.recapitate( ) instead.

        To allow recapitation to be done correctly, the nodes of the
        first generation of the SLiM simulation from whom all samples inherit
        are still present in the tree sequence, but are not marked as samples.
        If you simplify the tree sequence before recapitating you must ensure
        these are not removed, which you do by passing the argument
        ``keep_input_roots=True`` to :meth:`.simplify()`.

        Note that ``Ne`` is not set automatically, so defaults to ``1.0``; you probably
        want to set it explicitly.  Similarly, migration is not set up
        automatically, so that if there are uncoalesced lineages in more than
        one population, you will need to pass in a demography that allows coalescence.

        In general, all defaults are whatever the defaults of
        {ref}`msprime.simulate` are; this includes recombination rate, so
        that if neither ``recombination_rate`` or a ``recombination_map`` are
        provided, there will be *no* recombination.

        :param float recombination_rate: A (constant) recombination rate,
            in units of crossovers per nucleotide per unit of time.
        :param list population_configurations: See :meth:`msprime.simulate` for
            this argument; if not provided, each population will have zero growth rate
            and the same effective population size.
        :type recombination_map: :class`msprime.RecombinationMap`
        :param recombination_map: The recombination map, or None,
            if recombination_rate is specified.
        :param dict kwargs: Any other arguments to :meth:`msprime.simulate`.
        '''
        if "keep_first_generation" in kwargs:
            raise ValueError("The keep_first_generation argument is deprecated: "
                             "the FIRST_GEN flag is no longer used.")
        _deprecation_warning("Please use pyslim.recapitate( ) instead.")

        if recombination_map is None:
            recombination_map = msprime.RecombinationMap(
                       positions = [0.0, self.sequence_length],
                       rates = [recombination_rate, 0.0],
                       num_loci = int(self.sequence_length))

        if population_configurations is None:
            population_configurations = [msprime.PopulationConfiguration()
                                         for _ in range(self.num_populations)]
        # temporary until we remove support for setting slim_generation
        if self._slim_generation != self.metadata['SLiM']['generation']:
            tables = self.dump_tables()
            md = tables.metadata
            md['SLiM']['generation'] = self._slim_generation
            tables.metadata = md
            ts = tables.tree_sequence()
        else:
            ts = self

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", msprime.IncompletePopulationMetadataWarning)
            recap = msprime.simulate(
                        from_ts = ts,
                        population_configurations = population_configurations,
                        recombination_map = recombination_map,
                        start_time = self.metadata['SLiM']['generation'],
                        **kwargs)
        return SlimTreeSequence(recap)

    def mutation_at(self, node, position, time=None):
        '''
        Finds the mutation present in the genome of ``node`` at ``position``,
        returning -1 if there is no such mutation recorded in the tree
        sequence.  Warning: if ``node`` is not actually in the tree sequence
        (e.g., not ancestral to any samples) at ``position``, then this
        function will return -1, possibly erroneously.  If `time` is provided,
        returns the last mutation at ``position`` inherited by ``node`` that
        occurred at or before ``time`` ago.

        :param int node: The index of a node in the tree sequence.
        :param float position: A position along the genome.
        :param int time: The time ago that we want the nucleotide, or None,
            in which case the ``time`` of ``node`` is used.

        :returns: Index of the mutation in question, or -1 if none.
        '''
        if position < 0 or position >= self.sequence_length:
            raise ValueError("Position {} not valid.".format(position))
        if node < 0 or node >= self.num_nodes:
            raise ValueError("Node {} not valid.".format(node))
        if time is None:
            time = self.node(node).time
        tree = self.at(position)
        site_pos = self.tables.sites.position
        out = tskit.NULL
        if position in site_pos:
            site_index = np.where(site_pos == position)[0][0]
            site = self.site(site_index)
            mut_nodes = []
            # look for only mutations that occurred before `time`
            # not strictly necessary if time was None
            for mut in site.mutations:
                if mut.time >= time:
                    mut_nodes.append(mut.node)
            n = node
            while n > -1 and n not in mut_nodes:
                n = tree.parent(n)
            if n >= 0:
                # do careful error checking here
                for mut in site.mutations:
                    if mut.node == n:
                        assert(out == tskit.NULL or out == mut.parent)
                        out = mut.id
        return out

    def nucleotide_at(self, node, position, time=None):
        '''
        Finds the nucleotide present in the genome of ``node`` at ``position``.
        Warning: if ``node`` is not actually in the tree sequence (e.g., not
        ancestral to any samples) at ``position``, then this function will
        return the reference sequence nucleotide, possibly erroneously.  If
        `time` is provided, returns the last nucletide produced by a mutation
        at ``position`` inherited by ``node`` that occurred at or before
        ``time`` ago.

        :param int node: The index of a node in the tree sequence.
        :param float position: A position along the genome.
        :param int time: The time ago that we want the nucleotide, or None,
            in which case the ``time`` of ``node`` is used.

        :returns: Index of the nucleotide in ``NUCLEOTIDES`` (0=A, 1=C, 2=G, 3=T).
        '''
        if not self.has_reference_sequence():
            raise ValueError("This tree sequence has no reference sequence.")
        mut_id = self.mutation_at(node, position, time)
        if mut_id == tskit.NULL:
            out = NUCLEOTIDES.index(self.reference_sequence.data[int(position)])
        else:
            mut = self.mutation(mut_id)
            k = np.argmax([u["slim_time"] for u in mut.metadata["mutation_list"]])
            out = mut.metadata["mutation_list"][k]["nucleotide"]
        return out

    @property
    def slim_provenance(self):
        '''
        Returns model type, slim generation, and remembered node count from
        the last entry in the provenance table that is tagged with "program"="SLiM".

        NOTE: you probably want to use the ``.metadata`` property instead.

        :rtype ProvenanceMetadata:
        '''
        warnings.warn("The 'slim_provenance' attribute is deprecated: get information from "
                      "ts.metadata['SLiM'] instead.", FutureWarning)
        return get_provenance(self, only_last=True)

    @property
    def slim_provenances(self):
        '''
        Returns model type, slim generation, and remembered node count from *all*
        entries in the provenance table that are tagged with "program"="SLiM".

        :rtype ProvenanceMetadata:
        '''
        return get_provenance(self, only_last=False)

    def individuals_alive_at(self, time, stage='late', remembered_stage=None,
                             population=None, samples_only=False):
        """
        Returns an array giving the IDs of all individuals that are known to be
        alive at the given time ago.  This is determined using their birth time
        ago (given by their `time` attribute) and, for nonWF models,
        their `age` attribute (which is equal to their age at the last time
        they were Remembered). See also :meth:`.individual_ages_at`.

        In WF models, birth occurs after "early()", so that individuals are only
        alive during "late()" for the time step when they have age zero,
        while in nonWF models, birth occurs before "early()", so they are alive
        for both stages.
        
        In both WF and nonWF models, mortality occurs between
        "early()" and "late()", so that individuals are last alive during the
        "early()" stage of the time step of their final age, and if individuals
        are alive during "late()" they will also be alive during "early()" of the
        next time step. This means it is important to know during which stage
        individuals were Remembered - for instance, if the call to
        sim.treeSeqRememberIndividuals() was made during "early()" of a given time step,
        then those individuals might not have survived until "late()" of that
        time step. Since SLiM does not record the stage at which individuals
        were Remembered, you can specify this by setting ``remembered_stages``:
        it should be the stage during which *all* calls to
        ``sim.treeSeqRememberIndividuals()``,  as well as to ``sim.treeSeqOutput()``,
        were made.

        Note also that in nonWF models, birth occurs before "early()", so the
        possible parents in a given time step are those that are alive in
        "early()" and have age greater than zero, or, equivalently, are alive in
        "late()" during the previous time step.
        In WF models, birth occurs after "early()", so possible parents in a
        given time step are those that are alive during "early()" of that time
        step or are alive during "late()" of the previous time step.

        :param float time: The number of time steps ago.
        :param str stage: The stage in the SLiM life cycle that we are inquiring
            about (either "early" or "late"; defaults to "late").
        :param str remembered_stage: The stage in the SLiM life cycle
            during which individuals were Remembered (defaults to the stage the
            tree sequence was recorded at, stored in metadata).
        :param int population: If given, return only individuals in the
            population(s) with these population ID(s).
        :param bool samples_only: Whether to return only individuals who have at
            least one node marked as samples.
        """
        if stage not in ("late", "early"):
            raise ValueError(f"Unknown stage '{stage}': "
                              "should be either 'early' or 'late'.")

        if remembered_stage is None:
            remembered_stage = self.metadata['SLiM']['stage']

        if remembered_stage not in ("late", "early"):
            raise ValueError(f"Unknown remembered_stage '{remembered_stage}': "
                              "should be either 'early' or 'late'.")
        if remembered_stage != self.metadata['SLiM']['stage']:
            warnings.warn(f"Provided remembered_stage '{remembered_stage}' does not"
                          " match the stage at which the tree sequence was saved"
                          f" ('{self.metadata['SLiM']['stage']}'). This is not necessarily"
                          " an error, but mismatched stages will lead to inconsistencies:"
                          " make sure you know what you're doing.")

        # birth_time is the time ago that they were first alive in 'late'
        # in a nonWF model they are alive for the same time step's 'early'
        # but in a WF model the first 'early' they were alive for is one more recent
        birth_time = self.individual_times
        # birth_time - birth_offset is the first time ago they were alive
        # during stage 'stage'
        if stage == "early" and self.metadata['SLiM']['model_type'] == "WF":
            birth_offset = 1
        else:
            birth_offset = 0
        # ages is the number of complete life cycles they are known to have lived through,
        # and so individuals have lived through at least 'age + 1' of both stages.
        # In nonWF models, they live for one more 'early' than 'late',
        # but this is only reflected in their age if Remembered in 'early'.
        ages = self.individual_ages
        # ages + age_offset + 1 is the number of 'stage' stages they are known
        # to have lived through
        if (self.metadata['SLiM']['model_type'] == "WF"
                or stage == remembered_stage):
            age_offset = 0
        else:
            if (remembered_stage == "early"
                    and stage == "late"):
                age_offset = -1
            else:
                age_offset = 1
        # if adjusted age=0 then they are be alive at exactly one time step
        alive_bool = np.logical_and(
                birth_time >= time + birth_offset,
                birth_time - ages <= time + birth_offset + age_offset)

        if population is not None:
            alive_bool &= np.isin(self.individual_populations, population)
        if samples_only:
            alive_bool &= (0 < np.bincount(1 + self.tables.nodes.individual,
                                           self.tables.nodes.flags & tskit.NODE_IS_SAMPLE,
                                           minlength=1 + self.num_individuals)[1:])

        return np.where(alive_bool)[0]

    def individual_ages_at(self, time, stage="late", remembered_stage="late"):
        """
        Returns the `ages` of each individual at the corresponding time ago,
        which will be ``nan`` if the individual is either not born yet or dead.
        This is computed as the time ago the individual was born (found by the
        `time` associated with the the individual's nodes) minus the `time`
        argument; while "death" is inferred from the individual's ``age``,
        recorded in metadata. These values are the same as what would be shown
        in SLiM during the corresponding time step and stage.

        Since age increments at the end of each time step,
        the age is the number of time steps ends the individual has lived
        through, so if they were born in time step `time`, then their age
        will be zero.

        In a WF model, this method does not provide any more information than
        does :meth:`.individuals_alive_at`, but for consistency, non-nan ages
        will be 0 in "late" and 1 in "early".
        See :meth:`.individuals_alive_at` for further discussion.

        :param float time: The reference time ago.
        :param str stage: The stage in the SLiM life cycle used to determine who
            is alive (either "early" or "late"; defaults to "late").
        :param str remembered_stage: The stage in the SLiM life cycle during which
            individuals were Remembered.
        """
        ages = np.repeat(np.nan, self.num_individuals)
        alive = self.individuals_alive_at(time, stage=stage,
                                          remembered_stage=remembered_stage)
        ages[alive] = self.individual_times[alive] - time
        return ages

    def slim_time(self, time, stage="late"):
        """
        Converts the given "tskit times" (i.e., in units of time before the end
        of the simulation) to SLiM times (those recorded by SLiM, usually in units
        of generations since the start of the simulation). Although the latter are
        always integers, these will not be if the provided times are not integers.

        When the tree sequence is written out, SLiM records the value of its
        current generation, which can be found in the metadata:
        ``ts.metadata['SLiM']['generation']``. In most cases, the “SLiM time”
        referred to by a time ago in the tree sequence (i.e., the value that would
        be reported by sim.generation within SLiM at the point in time thus
        referenced) can be obtained by subtracting that time ago from
        ``ts.metadata['SLiM']['generation']``. However, in WF models, birth
        happens between the “early()” and “late()” stages, so if the tree
        sequence was written out using sim.treeSeqOutput() during “early()” in
        a WF model, the tree sequence’s times measure time before the last set
        of individuals are born, i.e., before SLiM time step
        ``ts.metadata['SLiM']['generation'] - 1``.

        In some situations (e.g., mutations added during early() in WF models)
        this may not return what you expect. See :ref:`sec_metadata_converting_times`
        for more discussion.

        :param array time: An array of times to be converted.
        :param string stage: The stage of the SLiM life cycle that the SLiM time
            should be computed for.
        """
        slim_time = self.metadata['SLiM']['generation'] - time
        if self.metadata['SLiM']['model_type'] == "WF":
            if (self.metadata['SLiM']['stage'] == "early" and stage == "late"):
                slim_time -= 1
            if (self.metadata['SLiM']['stage'] == "late" and stage == "early"):
                slim_time += 1
        return slim_time

    def first_generation_individuals(self):
        """
        Returns the IDs of the individuals remembered as part of the first SLiM generation,
        as determined by their flags.

        .. warning::
        
            This method is deprecated, because from SLiM version 3.5
            the first generation individuals are no longer marked as such:
            only tree sequences from older versions of SLiM will have
            these individuals.
        """
        warnings.warn(
                "This method is deprecated: SLiM no longer marks individuals as "
                "'first generation' any longer - you must explicitly Remember them "
                "to keep them in the tree sequence. Returning *retained* individuals "
                "instead.", FutureWarning)
        return np.where(self.tables.individuals.flags & INDIVIDUAL_RETAINED > 0)[0]

    def _do_individual_parents_stuff(self, return_parents=False):
        # Helper for has_individual_parents and individual_parents,
        # which share a lot of machinery.
        edges = self.tables.edges
        nodes = self.tables.nodes
        edge_parent_indiv = nodes.individual[edges.parent]
        edge_child_indiv = nodes.individual[edges.child]
        # nodes whose parent nodes are all in the same individual
        unique_parent_nodes = unique_labels_by_group(
                edges.child,
                edge_parent_indiv,
                minlength=nodes.num_rows)
        unique_parent_edges = unique_parent_nodes[edges.child]
        # edges describing relationships between individuals
        indiv_edges = np.logical_and(
                np.logical_and(edge_parent_indiv != tskit.NULL,
                                     edge_child_indiv != tskit.NULL),
                unique_parent_edges)
        # individual edges where the parent was alive during "late"
        # of the time step before the child is born
        child_births = self.individual_times[edge_child_indiv[indiv_edges]]
        parent_births = self.individual_times[edge_parent_indiv[indiv_edges]]
        alive_edges = indiv_edges.copy()
        if self.metadata['SLiM']['model_type'] == "WF":
            alive_edges[indiv_edges] = (child_births + 1 == parent_births)
        else:
            parent_deaths = parent_births - self.individual_ages[edge_parent_indiv[indiv_edges]]
            alive_edges[indiv_edges] = (child_births + 1 >= parent_deaths)
        edge_spans = edges.right - edges.left
        parental_span = np.bincount(edge_child_indiv[alive_edges],
                weights=edge_spans[alive_edges], minlength=self.num_individuals)
        # we could also check for edges without individual parents terminating
        # in this individual, but this is unnecessary as the entire genome is
        # accounted for
        has_all_parents = (parental_span == 2 * self.sequence_length)
        if return_parents:
            full_parent_edges = np.logical_and(
                    alive_edges,
                    has_all_parents[edge_child_indiv])
            parents = np.unique(np.column_stack(
                            [edge_parent_indiv[full_parent_edges],
                             edge_child_indiv[full_parent_edges]]
                            ), axis=0)
            return parents
        else:
            return has_all_parents

    def individual_parents(self):
        '''
        Finds all parent-child relationships in the tree sequence (as far as we
        can tell). The output will be a two-column array with row [i,j]
        indicating that individual i is a parent of individual j.  See
        :meth:`.has_individual_parents` for exactly which parents are returned.

        See :meth:`.individuals_alive_at` for further discussion about how
        this is determined based on when the individuals were Remembered.

        :return: An array of individual IDs, with row [i, j] if individual i is
            a parent of individual j.
        '''
        return self._do_individual_parents_stuff(return_parents=True)

    def has_individual_parents(self):
        '''
        Finds which individuals have both their parent individuals also present
        in the tree sequence, as far as we can tell. To do this, we return a
        boolean array with True for those individuals for which:

        - all edges terminating in that individual's nodes are in individuals,
        - each of the individual's nodes inherit from a single individual only,
        - those parental individuals were alive when the individual was born,
        - the parental individuals account for two whole genomes.

        This returns a boolean array indicating for each individual whether all
        these are true. Note in particular that individuals with only *one*
        recorded parent are *not* counted as "having parents".

        See :meth:`.individuals_alive_at` for further discussion about how
        this is determined based on when the individuals were Remembered.

        :return: A boolean array of length equal to ``targets``.
        '''
        return self._do_individual_parents_stuff(return_parents=False)


def _set_nodes_individuals(tables, age):
    '''
    Adds to a TableCollection the information relevant to individuals required
    for SLiM to load in a tree sequence, that is found in Node and Individual
    tables.  This will replace any existing Individual table, and will replace
    any information already in the individual, metadata, and population columns
    of the Node table.

    This is designed to make it easy to assign default values:
    - (node_ind) the 2*j-th and (2*j+1)-st `sample` nodes to individual j
    - (location) individual locations to (0, 0, 0)
    - (age) individual age to 0
    - (ind_id) SLiM individual pedigree IDs to sequential integers starting from 0
    - (ind_population) individual populations to 0
    - (node_id) SLiM genome IDs to sequential integers starting with samples from 0
    - (node_is_null) genomes to be non-null
    - (node_type) genome type to 0 (= autosome)
    - (ind_flags) INDIVIDUAL_ALIVE

    If you have other situations, like non-alive "remembered" individuals, you
    will need to edit the tables by hand, afterwards.
    '''
    samples = np.where(tables.nodes.flags & tskit.NODE_IS_SAMPLE)[0]
    if (len(samples) % 2) != 0:
        raise ValueError("There must be an even number of sampled nodes,"\
                         + "since organisms are diploid.")

    num_individuals = int(len(samples) / 2)
    node_ind = np.repeat(tskit.NULL, tables.nodes.num_rows).astype("int32")
    node_ind[samples] = np.arange(len(samples)) // 2
    ind_id = np.arange(num_individuals)
    slim_node_id = np.repeat(tskit.NULL, tables.nodes.num_rows)
    slim_node_id[samples] = np.arange(len(samples))

    ind_population = np.repeat(tskit.NULL, num_individuals)
    ind_population[node_ind[samples]] = tables.nodes.population[samples]

    if not np.all(unique_labels_by_group(node_ind,
                                          tables.nodes.population)):
        raise ValueError("Individual has nodes from more than one population.")
    if not np.all(unique_labels_by_group(node_ind,
                                          tables.nodes.time)):
        raise ValueError("Individual has nodes from more than one time.")

    loc_vec = np.zeros(num_individuals * 3).astype("float64")
    loc_off = 3 * np.arange(num_individuals + 1).astype("uint32")
    ind_flags = np.repeat(INDIVIDUAL_ALIVE, num_individuals).astype("uint32")

    default_ind = default_slim_metadata("individual")
    sex = default_ind['sex']
    slim_flag = default_ind['flags']
    pedigree_p1 = default_ind['pedigree_p1']
    pedigree_p2 = default_ind['pedigree_p2']

    ims = tables.individuals.metadata_schema
    individual_metadata = [
            ims.encode_row({
                'pedigree_id': iid,
                'age': age,
                'subpopulation': int(pop),
                'sex': sex,
                'flags': slim_flag,
                'pedigree_p1': pedigree_p1,
                'pedigree_p2': pedigree_p2,
            })
            for (iid, pop) in
            zip(ind_id, ind_population)]
    imb, imo = tskit.pack_bytes(individual_metadata)
    tables.individuals.set_columns(
            flags=ind_flags, location=loc_vec, location_offset=loc_off,
            metadata=imb, metadata_offset=imo)
    assert(tables.individuals.num_rows == num_individuals)
    
    default_node = default_slim_metadata("node")
    node_is_null = default_node["is_null"]
    node_type = default_node["genome_type"]
    nms = tables.nodes.metadata_schema
    node_metadata = [b'' for _ in range(tables.nodes.num_rows)]
    for j in samples:
        node_metadata[j] = nms.encode_row({'slim_id': slim_node_id[j],
                               'is_null': node_is_null,
                               'genome_type': node_type
                               })
    nmb, nmo = tskit.pack_bytes(node_metadata)
    tables.nodes.set_columns(flags=tables.nodes.flags, time=tables.nodes.time,
                             population=tables.nodes.population, individual=node_ind,
                             metadata=nmb, metadata_offset=nmo)


def _set_populations(tables):
    '''
    Adds to a TableCollection the information about populations required for SLiM
    to load a tree sequence. This will replace anything already in the Population
    table.
    '''
    num_pops = max(tables.nodes.population) + 1
    tables.populations.clear()
    pop_id = np.arange(num_pops)
    for j in range(num_pops):
        md = default_slim_metadata("population")
        md["slim_id"] = j
        md["name"] = f"pop_{j}"
        tables.populations.add_row(metadata=md)


def _set_sites_mutations(tables):
    '''
    Adds to a TableCollection the information relevant to mutations required
    for SLiM to load in a tree sequence. This means adding to the metadata column
    of the Mutation table,  It will also
    - give SLiM IDs to each mutation
    - replace ancestral states with ""
    This will replace any information already in the metadata or derived state
    columns of the Mutation table. We set slim_time in metadata so that
    - slim_generation = floor(tskit time) + slim_time
    '''
    num_mutations = tables.mutations.num_rows
    default_mut = default_slim_metadata("mutation_list_entry")
    dsb, dso = tskit.pack_bytes([str(j).encode() for j in range(num_mutations)])
    slim_time = tables.metadata["SLiM"]["generation"] - np.floor(tables.mutations.time).astype("int")
    mms = tables.mutations.metadata_schema
    mutation_metadata = [
            mms.encode_row(
                {"mutation_list": 
                 [{"mutation_type": default_mut["mutation_type"],
                  "selection_coeff": default_mut["selection_coeff"],
                  "subpopulation": default_mut["subpopulation"],
                  "slim_time": st,
                  "nucleotide": default_mut["nucleotide"]
                   }]
                })
            for st in slim_time]
    mdb, mdo = tskit.pack_bytes(mutation_metadata)
    tables.mutations.set_columns(
            site=tables.mutations.site,
            node=tables.mutations.node,
            time=tables.mutations.time,
            derived_state=dsb,
            derived_state_offset=dso,
            parent=tables.mutations.parent,
            metadata=mdb,
            metadata_offset=mdo)
    tables.sites.set_columns(
            position=tables.sites.position,
            ancestral_state=np.array([], dtype='int8'),
            ancestral_state_offset=np.zeros(tables.sites.num_rows + 1, dtype='uint32'))

def update_tables(tables):
    """
    Update tables produced by a previous verion of SLiM to the current file version.
    """
    # First we ensure we can find the file format version number
    # in top-level metadata. Then we proceed to fix up the tables as necessary.
    if not (isinstance(tables.metadata, dict) and 'SLiM' in tables.metadata):
        # Old versions kept information in provenance, not top-level metadata.
        # Note this uses defaults on keys not present in provenance,
        # which prior to 0.5 was everything but generation and model_type.
        values = default_slim_metadata('tree_sequence')['SLiM']
        prov = None
        file_version = 'unknown'
        # use only the last SLiM provenance
        for p in tables.provenances:
            is_slim, this_file_version = slim_provenance_version(p) 
            if is_slim:
                prov = p
                file_version = this_file_version
        values['file_version'] = file_version
        try:
            record = json.loads(prov.record)
            if file_version == "0.1":
                values['model_type'] = record['model_type']
                values['generation'] = record['generation']
            else:
                for k in values:
                    if k in record['parameters']:
                        values[k] = record['parameters'][k]
                values['generation'] = record['slim']['generation']
        except:
            raise ValueError("Failed to obtain metadata from provenance.")
        set_tree_sequence_metadata(tables, **values)

    file_version = tables.metadata['SLiM']['file_version']
    if file_version != slim_file_version:
        warnings.warn("This is a version {} SLiM tree sequence.".format(file_version) +
                      " When you write this out, " +
                      "it will be converted to version {}.".format(slim_file_version))

        # the only tables to have metadata schema changed thus far
        # are populations, individuals, and mutations:
        old_schema = _old_metadata_schema("population", file_version)
        if old_schema is not None:
            pops = tables.populations.copy()
            tables.populations.clear()
            pops.metadata_schema = old_schema
            new_schema = slim_metadata_schemas["population"]
            tables.populations.metadata_schema = new_schema
            defaults = default_slim_metadata("population")
            # just needs recoding
            for pop in pops:
                tables.populations.append(pop)

        old_schema = _old_metadata_schema("individual", file_version)
        if old_schema is not None:
            inds = tables.individuals.copy()
            tables.individuals.clear()
            inds.metadata_schema = old_schema
            new_schema = slim_metadata_schemas["individual"]
            tables.individuals.metadata_schema = new_schema
            defaults = default_slim_metadata("individual")
            d = {}
            for k in ["pedigree_p1", "pedigree_p2"]:
                d[k] = defaults[k]
            for ind in inds:
                md = ind.metadata
                md.update(d)
                tables.individuals.append(ind.replace(metadata=md))

        old_schema = _old_metadata_schema("mutation", file_version)
        if old_schema is not None:
            muts = tables.mutations.copy()
            tables.mutations.clear()
            muts.metadata_schema = old_schema
            tables.mutations.metadata_schema = slim_metadata_schemas["mutation"]
            for mut in muts:
                md = mut.metadata
                for ml in md['mutation_list']:
                    ml['nucleotide'] = -1
                tables.mutations.append(mut.replace(metadata=md))

        if file_version == "0.1":
            # shift times
            slim_generation = tables.metadata['SLiM']['generation']
            node_times = tables.nodes.time + slim_generation
            tables.nodes.set_columns(
                    flags=tables.nodes.flags,
                    time=node_times,
                    population=tables.nodes.population,
                    individual=tables.nodes.individual,
                    metadata=tables.nodes.metadata,
                    metadata_offset=tables.nodes.metadata_offset)
            migration_times = tables.migrations.time + slim_generation
            tables.migrations.set_columns(
                    left=tables.migrations.left,
                    right=tables.migrations.right,
                    node=tables.migrations.node,
                    source=tables.migrations.source,
                    dest=tables.migrations.dest,
                    time=migration_times)

        new_record = {
                    "schema_version": "1.0.0",
                    "software": {
                        "name": "pyslim",
                        "version": pyslim_version,
                        },
                    "parameters": {
                        "command": ["updrade_tables"],
                        "old_file_version": file_version,
                        "new_file_version": slim_file_version,
                        },
                    "environment": get_environment(),
                }
        tskit.validate_provenance(new_record)
        tables.provenances.add_row(json.dumps(new_record))

        set_metadata_schemas(tables)
        md = tables.metadata
        md['SLiM']['file_version'] = slim_file_version
        tables.metadata = md
