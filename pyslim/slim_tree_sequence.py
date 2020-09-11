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
from .slim_metadata import _decode_mutation_pre_nucleotides, _set_metadata_schemas

INDIVIDUAL_ALIVE = 2**16
INDIVIDUAL_REMEMBERED = 2**17
# no longer used but keep for a while
INDIVIDUAL_FIRST_GEN = 2**18

# A nucleotide k in mutation metadata actually means
# something that in reference_sequence is NUCLEOTIDES[k]
NUCLEOTIDES = ['A', 'C', 'G', 'T']


def load(path, legacy_metadata=False):
    '''
    Load the SLiM-compatible tree sequence found in the .trees file at ``path``.

    :param string path: The path to a .trees file.
    :param bool legacy_metadata: If True, then the resulting tree sequence will
        provide old-style metadata: as objects instead of dictionaries. This
        option is deprecated and will dissappear at some point in the future.
    '''
    ts = SlimTreeSequence.load(path, legacy_metadata=legacy_metadata)
    return ts


def load_tables(tables, **kwargs):
    '''
    See :func:`SlimTreeSequence.load_tables`.

    :param TableCollection tables: A set of tables.
    '''
    ts = SlimTreeSequence.load_tables(tables, **kwargs)
    return ts


def annotate_defaults(ts, model_type, slim_generation, reference_sequence=None):
    '''
    Takes a tree sequence (as produced by msprime, for instance), and adds in the
    information necessary for SLiM to use it as an initial state, filling in
    mostly default values. Returns a :class:`SlimTreeSequence`.

    :param TreeSequence ts: A :class:`TreeSequence`.
    :param string model_type: SLiM model type: either "WF" or "nonWF".
    :param int slim_generation: What generation number in SLiM correponds to
        ``time=0`` in the tree sequence.
    '''
    tables = ts.dump_tables()
    annotate_defaults_tables(tables, model_type, slim_generation)
    return SlimTreeSequence.load_tables(tables,
                reference_sequence=reference_sequence)


def annotate_defaults_tables(tables, model_type, slim_generation):
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
    top_metadata = default_slim_metadata['tree_sequence']['SLiM']
    top_metadata['model_type'] = model_type
    top_metadata['generation'] = slim_generation
    set_tree_sequence_metadata(tables, **top_metadata)
    _set_nodes_individuals(tables, age=default_ages)
    _set_populations(tables)
    _set_sites_mutations(tables)


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
                    "metadata interface: see "
                    "`the documentation <https://pyslim.readthedocs.io/en/latest/metadata.html#legacy-metadata>`_ "
                    "for how to switch over your script")
        else:
            raise AttributeError(f"'dict' object has no attribute '{name}'")

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError as e:
            if isinstance(key, int):
                msg = e.args[0]
                e.args = (f"{msg}: It looks like you're trying to use the legacy "
                           "metadata interface: see "
                           "`the documentation <https://pyslim.readthedocs.io/en/latest/metadata.html#legacy-metadata>`_ "
                           "for how to switch over your script",)
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

    All mutable properties of individuals (e.g., age) is as it was recorded during
    the individual's last time step alive (or at the end of the simulation, if they
    are still alive).

    You can create a :class:`.SlimTreeSequence` using one of

    - :meth:`.SlimTreeSequence.load_tables` :meth:`.SlimTreeSequence.load`,
    - :func:`.load`, or :func:`.load_tables`.

    :ivar slim_generation: The generation that the SLiM simulation was at upon writing;
        will be read from metadata if not provided.
    :ivar reference_sequence: None, or an string of length equal to the sequence
        length that gives the entire reference sequence for nucleotide models.
    :ivar legacy_metadata: Whether this tree sequence returns metadata in objects
        (as in older versions of pyslim) rather than dicts: see
        `the documentation <https://pyslim.readthedocs.io/en/latest/metadata.html#legacy-metadata>`_.
        This option is deprecated and will disappear at some point.
    :vartype slim_generation: int
    :vartype reference_sequence: string
    '''

    def __init__(self, ts, reference_sequence=None, legacy_metadata=False):
        self.legacy_metadata = legacy_metadata
        if not (isinstance(ts.metadata, dict) and 'SLiM' in ts.metadata
                and ts.metadata['SLiM']['file_version'] == slim_file_version):
            tables = ts.dump_tables()
            if not (isinstance(tables.metadata, dict) and 'SLiM' in tables.metadata):
                _set_metadata_from_provenance(tables)
            if tables.metadata['SLiM']['file_version'] != slim_file_version:
                _upgrade_old_tables(tables)
                # cannot assign directly to keys of metadata
                md = tables.metadata
                md['SLiM']['file_version'] = slim_file_version
                tables.metadata = md
            ts = tables.tree_sequence()
        slim_generation = ts.metadata['SLiM']['generation']
        self.model_type = ts.metadata['SLiM']['model_type']
        super().__init__(ts._ll_tree_sequence)
        self.slim_generation = slim_generation
        self.reference_sequence = reference_sequence
        # pre-extract individual metadata
        self.individual_locations = ts.tables.individuals.location
        self.individual_locations.shape = (int(len(self.individual_locations)/3), 3)
        self.individual_ages = np.zeros(ts.num_individuals, dtype='int')
        if self.model_type != "WF":
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

    @classmethod
    def load(cls, path, legacy_metadata=False):
        '''
        Load a :class:`SlimTreeSequence` from a .trees file on disk.

        :param string path: The path to a .trees file.
        :rtype SlimTreeSequence:
        '''
        ts = tskit.load(path)
        # extract the reference sequence from the kastore
        kas = kastore.load(path)
        if 'reference_sequence/data' in kas:
            int_rs = kas['reference_sequence/data']
            reference_sequence = int_rs.tostring().decode('ascii')
        else:
            reference_sequence = None
        return cls(ts, reference_sequence=reference_sequence, legacy_metadata=legacy_metadata)

    @classmethod
    def load_tables(cls, tables, **kwargs):
        '''
        Creates the :class:`SlimTreeSequence` defined by the tables.

        :param TableCollection tables: A set of tables, as produced by SLiM
            or by annotate_defaults().
        :param TableCollection reference_sequence: An optional string of ACGT giving
            the reference sequence.
        :rtype SlimTreeSequence:
        '''
        # a roundabout way to copy the tables
        ts = tables.tree_sequence()
        return cls(ts, **kwargs)

    def dump(self, path, **kwargs):
        '''
        Dumps the tree sequence to the path specified. This is mostly just a wrapper for
        tskit.TreeSequence.dump(), but also writes out the reference sequence.

        :param str path: The file path to write the TreeSequence to.
        :param kwargs: Additional keyword args to pass to tskit.TreeSequence.dump
        '''
        super(SlimTreeSequence, self).dump(path, **kwargs)
        if self.reference_sequence is not None:
            # to convert to a kastore store we need to reload from a file,
            # and for it to be mutable we need to make it a dict
            kas = dict(kastore.load(path))
            kas['reference_sequence/data'] = np.frombuffer(self.reference_sequence.encode(),
                                                           dtype=np.int8)
            kastore.dump(kas, path)

    def simplify(self, *args, **kwargs):
        '''
        This is a wrapper for :meth:`tskit.TreeSequence.simplify`.
        The only difference is that this method returns the
        derived class :class:`.SlimTreeSequence`.

        If you have not yet recapitated your SlimTreeSequence, you probably want to
        pass ``keep_input_roots=True``, so that recapitation is possible in the future.

        :rtype SlimTreeSequence:
        '''
        sts = super(SlimTreeSequence, self).simplify(*args, **kwargs)
        if (type(sts) == tuple):
            ret = (SlimTreeSequence(sts[0]), sts[1])
            ret[0].reference_sequence = self.reference_sequence
        else:
            ret = SlimTreeSequence(sts)
            ret.reference_sequence = self.reference_sequence

        return ret

    def population(self, id_):
        '''
        Returns the population whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.population`, but with additional attributes::

            slim_id, selfing_fraction, female_cloning_fraction,
            male_cloning_fraction, sex_ratio,
            bounds_x0, bounds_x1, bounds_y0, bounds_y1, bounds_z0, bounds_z1,
            migration_records.

        These are all recorded by SLiM in the metadata.

        Note that SLiM populations are usually indexed starting from 1,
        but in tskit from zero, so there may be populations (e.g., with id_=0)
        that have no metadata and are not used by SLiM.

        :param int id_: The ID of the population (i.e., its index).
        '''
        pop = super(SlimTreeSequence, self).population(id_)
        if self.legacy_metadata:
            try:
                pop.metadata = PopulationMetadata.fromdict(pop.metadata)
            except:
                pass
        else:
            if pop.metadata is not None:
                pop.metadata = MetadataDictWrapper(pop.metadata)
        return pop

    def individual(self, id_):
        '''
        Returns the individual whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.individual`, but with additional attributes::

          time, pedigree_id, age, slim_population, sex, slim_flags.

        The `time` and `population` properties are extracted from the nodes,
        and an error will be thrown if the individual's nodes derive from
        more than one population or more than one time.

        :param int id_: The ID of the individual (i.e., its index).
        '''
        ind = super(SlimTreeSequence, self).individual(id_)
        ind.population = self.individual_populations[id_]
        ind.time = self.individual_times[id_]
        if self.legacy_metadata:
            try:
                ind.metadata = IndividualMetadata.fromdict(ind.metadata)
            except:
                pass
        else:
            ind.metadata = MetadataDictWrapper(ind.metadata)
        return ind

    def node(self, id_):
        '''
        Returns the node whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.node`, but with additional attributes::

           slim_id, is_null, genome_type.

        These are all recorded by SLiM in the metadata.

        :param int id_: The ID of the node (i.e., its index).
        '''
        node = super(SlimTreeSequence, self).node(id_)
        if self.legacy_metadata:
            try:
                node.metadata = NodeMetadata.fromdict(node.metadata)
            except:
                pass
        else:
            if node.metadata is not None:
                node.metadata = MetadataDictWrapper(node.metadata)
        return node

    def mutation(self, id_):
        '''
        Returns the mutation whose ID is given by `id_`, as documented in
        :meth:`tskit.TreeSequence.mutation`, but with additional attributes::

           mutation_type, selection_coeff, population, slim_time, nucleotide.

        These are all recorded by SLiM in the metadata.

        :param int id_: The ID of the mutation (i.e., its index).
        '''
        mut = super(SlimTreeSequence, self).mutation(id_)
        if self.legacy_metadata:
            try:
                mut.metadata = [MutationMetadata.fromdict(x) for x in mut.metadata['mutation_list']]
            except:
                pass
        else:
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

        To allow recapitation to be done correctly, the nodes of the
        first generation of the SLiM simulation from whom all samples inherit
        are still present in the tree sequence, but are not marked as samples.
        If you simplify the tree sequence before recapitating you must ensure
        these are not removed, which you do by passing the argument
        ``keep_input_roots=True`` to :meth:`.simplify()`.

        Note that ``Ne`` is not set automatically, so defaults to ``1.0``; you probably
        want to set it explicitly.  Similarly, migration is not set up
        automatically, so that if there are uncoalesced lineages in more than
        one population, you will need to pass in a migration matrix to allow
        coalescence. In both cases, remember that population IDs in ``tskit`` begin
        with 0, so that if your SLiM simulation has populations ``p1`` and ``p2``,
        then the tree sequence will have three populations (but with no nodes
        assigned to population 0), so that migration rate of 1.0 between ``p1`` and
        ``p2`` needs a migration matrix of::

           [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]

        In general, all defaults are whatever the defaults of ``msprime.simulate`` are;
        this includes recombination rate, so that if neither ``recombination_rate``
        or a ``recombination_map`` are provided, there will be *no* recombination.

        However, if ``recombination_rate`` *is* provided, then recapitation will
        use a constant rate of recombination on a discretized map -- in other words,
        recombinations in the coalescent portion of the simulation will only occur
        at integer locations, just as in SLiM. If you do not want this to happen,
        you need to construct a ``recombination_map`` explicitly.

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
            raise ValueError("The keep_first_generation argument is deprecated:"
                             "the FIRST_GEN flag is no longer used.")

        # toggle for hacks below to deal with old msprime
        discrete_msprime = hasattr(msprime, "RateMap")

        if recombination_rate is not None:
            if recombination_map is not None:
                raise ValueError("Cannot specify length/recombination_rate along with a recombination map")
            if discrete_msprime:
                recombination_map = msprime.RecombinationMap(positions = [0.0, self.sequence_length],
                                                             rates = [recombination_rate, 0.0])
            else:
                recombination_map = msprime.RecombinationMap(positions = [0.0, self.sequence_length],
                                                             rates = [recombination_rate, 0.0],
                                                             num_loci = int(self.sequence_length))
        if discrete_msprime and ('discrete_genome' not in kwargs):
            kwargs['discrete_genome'] = True

        if population_configurations is None:
            population_configurations = [msprime.PopulationConfiguration()
                                         for _ in range(self.num_populations)]

        recap = msprime.simulate(
                    from_ts = self,
                    population_configurations = population_configurations,
                    recombination_map = recombination_map,
                    start_time = self.slim_generation,
                    **kwargs)
        # HACK to deal with older msprime that doesn't retain metadata
        # by copying over the metadata
        if recap.metadata == b'':
            tables = recap.tables
            tables.metadata = self._ll_tree_sequence.get_metadata()
            tables.metadata_schema = self.metadata_schema
            _set_metadata_schemas(tables)
            recap = tables.tree_sequence()
        return SlimTreeSequence(recap, reference_sequence=self.reference_sequence)

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
        if self.reference_sequence is None:
            raise ValueError("This tree sequence has no reference sequence.")
        mut_id = self.mutation_at(node, position, time)
        if mut_id == tskit.NULL:
            out = NUCLEOTIDES.index(self.reference_sequence[int(position)])
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
        warnings.warn("This is deprecated: get information from "
                      "ts.metadata['SLiM'] instead.", FutureWarning)
        return get_provenance(self, only_last=True)

    @property
    def slim_provenances(self):
        '''
        Returns model type, slim generation, and remembered node count from *all*
        entries in the provenance table that is tagged with "program"="SLiM"

        :rtype ProvenanceMetadata:
        '''
        return get_provenance(self, only_last=False)

    def _mark_first_generation(self):
        '''
        Mark all 'first generation' individuals' nodes as samples, and return
        the corresponding tree sequence.

        DEPRECATED
        '''
        tables = self.dump_tables()
        first_gen_nodes = ((tables.nodes.individual > 0)
                           & ((tables.individuals.flags[tables.nodes.individual]
                               & INDIVIDUAL_FIRST_GEN) > 0))
        if sum(first_gen_nodes) == 0:
            warnings.warn("Tree sequence does not have the initial generation; " +
                          " did you simplify it after output from SLiM?")
        flags = tables.nodes.flags
        flags[first_gen_nodes] = (flags[first_gen_nodes] | tskit.NODE_IS_SAMPLE)
        tables.nodes.set_columns(flags=flags, population=tables.nodes.population,
                individual=tables.nodes.individual, time=tables.nodes.time,
                metadata=tables.nodes.metadata,
                metadata_offset=tables.nodes.metadata_offset)
        ts = load_tables(tables)
        ts.reference_sequence = self.reference_sequence
        return ts

    def individuals_alive_at(self, time, stage='late', remembered_stage=None):
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
        it should be the stage during which *all* calls to sim.treeSeqRememberIndividuals,
        as well as to sim.treeSeqOutput(), were made.

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
            that individuals were Remembered during (defaults to the stage the
            tree sequence was recorded at, stored in metadata).
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
        return np.where(alive_bool)[0]

    def individual_ages_at(self, time, stage="late", remembered_stage="late"):
        """
        Returns the *ages* of each individual at the corresponding time ago,
        which will be `nan` if the individual is either not born yet or dead.
        This is computed as the time ago the individual was born (found by the
        `time` associated with the the individual's nodes) minus the `time`
        argument; while "death" is inferred from the individual's `age`,
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
        :param str remembered_stage: The stage in the SLiM life cycle that
            individuals were Remembered during.
        """
        ages = np.repeat(np.nan, self.num_individuals)
        alive = self.individuals_alive_at(time, stage=stage,
                                          remembered_stage=remembered_stage)
        ages[alive] = self.individual_times[alive] - time
        return ages

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
                "to retain them in the tree sequence.", FutureWarning)
        return np.where(self.tables.individuals.flags & INDIVIDUAL_FIRST_GEN > 0)[0]

    def has_individual_parents(self):
        '''
        Finds which individuals have both their parent individuals also present
        in the tree sequence, as far as we can tell. To do this, we return the
        IDs of individuals for which:

        - all edges terminating in that individual's nodes are in individuals,
        - each of the individual's nodes inherit from a single individual only,
        - those parental individuals were alive when the individual was born,
        - the parental individuals account for two whole genomes.

        This returns a boolean array indicating for each individual whether all
        these are true.

        See :meth:`.individuals_alive_at` for further discussion about how
        this is determined based on when the individuals were Remembered.

        :return: A boolean array of length equal to ``targets``.
        '''
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
        # total genome inherited from parents
        edge_spans = edges.right - edges.left
        parental_span = np.bincount(edge_child_indiv[alive_edges],
                weights=edge_spans[alive_edges], minlength=self.num_individuals)
        # we could also check for edges without individual parents terminating
        # in this individual, but this is unnecessary as the entire genome is
        # accounted for
        has_all_parents = (parental_span == 2 * self.sequence_length)
        return has_all_parents


def _set_metadata_from_provenance(tables):
    # note this uses defaults on keys not present in provenance,
    # which prior to 0.5 was everything but generation and model_type
    values = default_slim_metadata['tree_sequence']['SLiM']
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


def _upgrade_old_tables(tables):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        provenance = get_provenance(tables)
    file_version = provenance.file_version
    slim_generation = provenance.slim_generation
    warnings.warn("This is an version {} SLiM tree sequence.".format(file_version) +
                  " When you write this out, " +
                  "it will be converted to version {}.".format(slim_file_version))
    if file_version == "0.1" or file_version == "0.2":
        # add empty nucleotide slots to metadata
        mut_bytes = tskit.unpack_bytes(tables.mutations.metadata,
                                       tables.mutations.metadata_offset)
        mut_metadata = [_decode_mutation_pre_nucleotides(md)
                        for md in mut_bytes]
        metadata, metadata_offset = tskit.pack_bytes(mut_metadata)
        tables.mutations.set_columns(
                site=tables.mutations.site,
                node=tables.mutations.node,
                parent=tables.mutations.parent,
                derived_state=tables.mutations.derived_state,
                derived_state_offset=tables.mutations.derived_state_offset,
                metadata=metadata,
                metadata_offset=metadata_offset)
    if file_version == "0.1":
        # shift times
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
                    "command": ["_upgrade_old_tables"],
                    "old_file_version": file_version,
                    "new_file_version": slim_file_version,
                    },
                "environment": get_environment(),
            }
    tskit.validate_provenance(new_record)
    tables.provenances.add_row(json.dumps(new_record))


def _set_nodes_individuals(
        tables, node_ind=None, location=(0, 0, 0), age=0, ind_id=None,
        ind_population=None, ind_sex=INDIVIDUAL_TYPE_HERMAPHRODITE,
        ind_flags=INDIVIDUAL_ALIVE, slim_ind_flags=0, node_id=None,
        node_is_null=False, node_type=GENOME_TYPE_AUTOSOME):
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
    samples = list(filter(lambda j: tables.nodes.flags[j] & tskit.NODE_IS_SAMPLE,
                          range(tables.nodes.num_rows)))
    if (len(samples) % 2) != 0:
        raise ValueError("There must be an even number of sampled nodes,"\
                         + "since organisms are diploid.")

    if node_ind is None:
        node_ind = [tskit.NULL for _ in range(tables.nodes.num_rows)]
        for j, k in enumerate(samples):
            node_ind[j] = int(k/2)

    num_individuals = max(node_ind) + 1
    num_nodes = tables.nodes.num_rows

    if type(location) is tuple:
        location = [location for _ in range(num_individuals)]
    assert(len(location) == num_individuals)

    if type(age) is int or type(age) is float:
        age = [age for _ in range(num_individuals)]
    assert(len(age) == num_individuals)

    if ind_id is None:
        ind_id = list(range(num_individuals))
    assert(len(ind_id) == num_individuals)

    if type(ind_sex) is int:
        ind_sex = [ind_sex for _ in range(num_individuals)]
    assert(len(ind_sex) == num_individuals)

    if type(slim_ind_flags) is int:
        slim_ind_flags = [slim_ind_flags for _ in range(num_individuals)]
    assert(len(slim_ind_flags) == num_individuals)

    if type(ind_flags) is int:
        ind_flags = [ind_flags for _ in range(num_individuals)]
    assert(len(ind_flags) == num_individuals)

    if node_id is None:
        node_id = [-1 for _ in range(num_nodes)]
        for j, k in enumerate(list(samples)
                              + sorted(list(set(range(num_nodes))
                                            - set(samples)))):
            node_id[k] = j
    assert(len(node_id) == num_nodes)

    if type(node_is_null) is bool:
        node_is_null = [node_is_null for _ in range(num_nodes)]
    assert(len(node_is_null) == num_nodes)

    if type(node_type) is int:
        node_type = [node_type for _ in range(num_nodes)]
    assert(len(node_type) == tables.nodes.num_rows)

    if ind_population is None:
        # set the individual populations based on what's in the nodes
        ind_population = [tskit.NULL for _ in range(num_individuals)]
        for j, u in enumerate(node_ind):
            if u >= 0:
                ind_population[u] = tables.nodes.population[j]
    assert(len(ind_population) == num_individuals)

    # check for consistency: every individual has two nodes, and populations agree
    ploidy = [0 for _ in range(num_individuals)]
    for j in samples:
        u = node_ind[j]
        assert(u >= 0)
        ploidy[u] += 1
        if tables.nodes.population[j] != ind_population[u]:
            raise ValueError("Inconsistent populations: nodes and individuals do not agree.")

    if any([p != 2 for p in ploidy]):
        raise ValueError("Not all individuals have two assigned nodes.")

    loc_vec, loc_off = tskit.pack_bytes(location)

    ims = tables.individuals.metadata_schema
    individual_metadata = [
            {'pedigree_id': iid, 'age': a, 'subpopulation': int(pop), 'sex': sex, 'flags': f}
            for (iid, a, pop, sex, f) in
            zip(ind_id, age, ind_population, ind_sex, slim_ind_flags)]
    assert(len(individual_metadata) == num_individuals)
    individual_metadata, individual_metadata_offset = tskit.pack_bytes(
            [ims.encode_row(r) for r in individual_metadata])
    tables.individuals.set_columns(
            flags=ind_flags, location=loc_vec, location_offset=loc_off,
            metadata=individual_metadata,
            metadata_offset=individual_metadata_offset)
    assert(tables.individuals.num_rows == num_individuals)
    
    node_metadata = [None for _ in range(num_nodes)]
    for j in samples:
        node_metadata[j] = {'slim_id': node_id[j],
                            'is_null': node_is_null[j],
                            'genome_type': node_type[j]
                            }
    nms = tables.nodes.metadata_schema
    node_metadata, node_metadata_offset = tskit.pack_bytes(
            [nms.encode_row(r) for r in node_metadata])
    tables.nodes.set_columns(flags=tables.nodes.flags, time=tables.nodes.time,
                             population=tables.nodes.population, individual=node_ind,
                             metadata=node_metadata,
                             metadata_offset=node_metadata_offset)


def _set_populations(
        tables, pop_id=None, selfing_fraction=0.0, female_cloning_fraction=0.0,
        male_cloning_fraction=0.0, sex_ratio=0.5, bounds_x0=0.0, bounds_x1=0.0,
        bounds_y0=0.0, bounds_y1=0.0, bounds_z0=0.0, bounds_z1=0.0,
        migration_records=None):
    '''
    Adds to a TableCollection the information about populations required for SLiM
    to load a tree sequence. This will replace anything already in the Population
    table.
    '''
    num_pops = max(tables.nodes.population) + 1
    for ind in tables.individuals:
        md = ind.metadata
        if not (isinstance(md, dict) and 'subpopulation' in md):
            raise ValueError("Individuals do not have valid metadata: "
                    "need to run set_nodes_individuals() first?")
        if md['subpopulation'] >= num_pops:
            raise ValueError("Bad population in individual metadata.")
    if pop_id is None:
        pop_id = list(range(num_pops))
    assert(len(pop_id) == num_pops)

    if type(selfing_fraction) is float:
        selfing_fraction = [selfing_fraction for _ in range(num_pops)]
    assert(len(selfing_fraction) == num_pops)

    if type(female_cloning_fraction) is float:
        female_cloning_fraction = [female_cloning_fraction for _ in range(num_pops)]
    assert(len(female_cloning_fraction) == num_pops)

    if type(male_cloning_fraction) is float:
        male_cloning_fraction = [male_cloning_fraction for _ in range(num_pops)]
    assert(len(male_cloning_fraction) == num_pops)

    if type(sex_ratio) is float:
        sex_ratio = [sex_ratio for _ in range(num_pops)]
    assert(len(sex_ratio) == num_pops)

    if type(bounds_x0) is float:
        bounds_x0 = [bounds_x0 for _ in range(num_pops)]
    assert(len(bounds_x0) == num_pops)

    if type(bounds_x1) is float:
        bounds_x1 = [bounds_x1 for _ in range(num_pops)]
    assert(len(bounds_x1) == num_pops)

    if type(bounds_y0) is float:
        bounds_y0 = [bounds_y0 for _ in range(num_pops)]
    assert(len(bounds_y0) == num_pops)

    if type(bounds_y1) is float:
        bounds_y1 = [bounds_y1 for _ in range(num_pops)]
    assert(len(bounds_y1) == num_pops)

    if type(bounds_z0) is float:
        bounds_z0 = [bounds_z0 for _ in range(num_pops)]
    assert(len(bounds_z0) == num_pops)

    if type(bounds_z1) is float:
        bounds_z1 = [bounds_z1 for _ in range(num_pops)]
    assert(len(bounds_z1) == num_pops)

    if migration_records is None:
        migration_records = [[] for _ in range(num_pops)]
    assert(len(migration_records) == num_pops)
    for mrl in migration_records:
        for mr in mrl:
            assert(isinstance(mr, dict))

    population_metadata = [
            {
                "slim_id" : pid,
                "selfing_fraction" : sf,
                "female_cloning_fraction" : fcf,
                "male_cloning_fraction" : mcf,
                "sex_ratio" : sr,
                "bounds_x0" : x0, "bounds_x1" : x1,
                "bounds_y0" : y0, "bounds_y1" : y1,
                "bounds_z0" : z0, "bounds_z1" : z1,
                "migration_records" : mr,
            } for (pid, sf, fcf, mcf, sr, x0, x1, y0, y1, z0, z1, mr) in
            zip(pop_id, selfing_fraction, female_cloning_fraction,
                male_cloning_fraction, sex_ratio, bounds_x0,
                bounds_x1, bounds_y0, bounds_y1, bounds_z0, bounds_z1,
                migration_records)]
    ms = tables.populations.metadata_schema
    tables.populations.packset_metadata(
            [ms.encode_row(r) for r in population_metadata])


def _set_sites_mutations(
        tables, mutation_id=None, mutation_type=1, selection_coeff=0.0,
        population=tskit.NULL, slim_time=None):
    '''
    Adds to a TableCollection the information relevant to mutations required
    for SLiM to load in a tree sequence. This means adding to the metadata column
    of the Mutation table,  It will also
    - give SLiM IDs to each mutation
    - round Site positions to integer values
    - stack any mutations that end up at the same position as a result
    - replace ancestral states with ""
    This will replace any information already in the metadata or derived state
    columns of the Mutation table.
    '''
    num_mutations = tables.mutations.num_rows

    if mutation_id is None:
        mutation_id = list(range(num_mutations))
    assert(len(mutation_id) == num_mutations)

    if type(mutation_type) is int:
        mutation_type = [mutation_type for _ in range(num_mutations)]
    assert(len(mutation_type) == num_mutations)

    if type(selection_coeff) is float:
        selection_coeff = [selection_coeff for _ in range(num_mutations)]
    assert(len(selection_coeff) == num_mutations)

    if type(population) is int:
        population = [population for _ in range(num_mutations)]
    assert(len(population) == num_mutations)

    if slim_time is None:
        ## This may *not* make sense because we have to round:
        # slim_time = [(-1) * int(tables.nodes.time[u]) for u in tables.mutations.node]
        slim_time = [0 for _ in range(num_mutations)]
    assert(len(slim_time) == num_mutations)

    mutation_metadata = [
            {"mutation_list": 
                [{"mutation_type": mt,
                  "selection_coeff": sc,
                  "subpopulation": pop,
                  "slim_time": st,
                  "nucleotide": -1
                  }]
            }
            for (mt, sc, pop, st) in
                         zip(mutation_type, selection_coeff, population, slim_time)]
    ms = tables.mutations.metadata_schema
    tables.mutations.packset_metadata(
            [ms.encode_row(r) for r in mutation_metadata])
