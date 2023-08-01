import msprime
import tskit
import warnings
import numpy as np

from .slim_tree_sequence import *
from .slim_metadata import *
from .provenance import *
from .util import *


def recapitate(ts,
               ancestral_Ne=None,
               **kwargs
    ):
    '''
    Returns a "recapitated" tree sequence, by using msprime to run a
    coalescent simulation from the "top" of this tree sequence, i.e.,
    allowing any uncoalesced lineages to coalesce.

    To allow recapitation to be done correctly, the nodes of the
    first generation of the SLiM simulation from whom all samples inherit
    are still present in the tree sequence, but are not marked as samples.
    If you simplify the tree sequence before recapitating you must ensure
    these are not removed, which you do by passing the argument
    ``keep_input_roots=True`` to :meth:`simplify <tskit.TreeSequence.simplify>`.

    If you specify an ``ancestral_Ne``, then the recapitated portion of the
    tree sequence will be simulated in a single population with this
    (diploid) size. In other words, all lineages are moved to a single
    population of this size (named "ancestral" if this name is not already
    taken), and coalescence is allowed to happen.

    You may control the ancestral demography by passing in a ``demography``
    argument: see :func:`msprime.sim_ancestry`.
    
    In general, all defaults are whatever the defaults of
    {meth}`msprime.sim_ancestry` are; this includes recombination rate, so
    that if neither ``recombination_rate`` or a ``recombination_map`` are
    provided, there will be *no* recombination.

    :param tskit.TreeSequence ts: The tree sequence to transform.
    :param float ancestral_Ne: If specified, then will simulate from a single
        ancestral population of this size. It is an error to specify this
        as well as ``demography``.
    :param dict kwargs: Any other arguments to :func:`msprime.sim_ancestry`.
    '''
    is_current_version(ts, _warn=True)
    if ancestral_Ne is not None:
        if "demography" in kwargs:
            raise ValueError("You cannot specify both `demography` and `ancestral_Ne`.")
        recap_time = ts.metadata['SLiM']['tick']
        # In various circumstances depending on the stage in which the simulation was
        # started and when the tree sequence was written out, the time of the roots
        # might be one or even two less than the "tick" value. We have access
        # to the stage the tree sequence was written out, but not the time it
        # was *started*. So, we'll just check if we need to subtract one or two.
        # Consistency checking this requires looping over all the trees, unfortunately,
        # but it avoids some common errors.
        root_times = list(set([ts.node(n).time for t in ts.trees() for n in t.roots]))
        for adj in (1, 2):
            if np.abs(recap_time - root_times[0] - adj) < 1e-8:
                recap_time -= adj
        if len(root_times) > 1 or not np.abs(root_times[0] - recap_time) < 1e-8:
            message = (
                "Not all roots of the provided tree sequence are at the time expected "
                "by recapitate(). This could happen if you've simplified in "
                "python before recapitating (fix: don't simplify first). "
                "If could also happen in other situations, e.g., "
                "you added new individuals without parents in SLiM "
                "during the course of the simulation with sim.addSubPop(), "
                "in which case you will probably need to recapitate with "
                "msprime.sim_ancestry(initial_state=ts, ...). "
                f"(Expected root time: {recap_time}; "
                f"Observed root times: {root_times})"
            )
            raise ValueError(message)
        demography = msprime.Demography.from_tree_sequence(ts)
        # must set pop sizes to >0 even though we merge immediately
        for pop in demography.populations:
            pop.initial_size=1.0
        ancestral_name = "ancestral"
        derived_names = [pop.name for pop in demography.populations]
        while ancestral_name in derived_names:
            ancestral_name = (ancestral_name + "_ancestral")
        demography.add_population(
                name=ancestral_name,
                description="ancestral population simulated by msprime",
                initial_size=ancestral_Ne,
        )
        # the split has to come slightly longer ago than slim's tick
        # since that's when all the linages are at, and otherwise the event
        # won't apply to them
        demography.add_population_split(
                np.nextafter( recap_time, 2 * recap_time),
                derived=derived_names,
                ancestral=ancestral_name,
        )
        kwargs["demography"] = demography

    recap = msprime.sim_ancestry(
                initial_state = ts,
                **kwargs)

    return recap


def convert_alleles(ts):
    """
    Returns a modified tree sequence in which alleles have been replaced by
    their corresponding nucleotides. For sites, SLiM-produced tree sequences
    have "" (the empty string) for the ancestral state at each site; this method
    will replace this with the corresponding nucleotide from the reference sequence.
    For mutations, SLiM records the 'derived state' as a SLiM mutation ID; this
    method will this with the nucleotide from the mutation's metadata.

    This operation is not reversible: since SLiM mutation IDs are lost, the tree
    sequence will not be able to be read back into SLiM.

    The main purpose of this method is for output: for instance, this code will produce
    a VCF file with nucleotide alleles:

    .. code-block:: python

        nts = pyslim.convert_alleles(ts)
        with open('nucs.vcf', 'w') as f:
            nts.write_vcf(f)

    This method will produce an error if the tree sequence does not have a
    valid reference sequence or if any mutations do not have nucleotides: to first
    generate these, see :func:`.generate_nucleotides`.

    :param tskit.TreeSequence ts: The tree sequence to transform.
    """
    tables = ts.dump_tables()
    has_refseq = (
            ts.has_reference_sequence()
            and len(ts.reference_sequence.data) > 0
    )
    if not has_refseq:
        raise ValueError("Tree sequence must have a valid reference sequence.")
    # unfortunately, nucleotide mutations may be stacked (e.g., substitutions
    # will appear this way) and they don't appear in any particular order;
    # so we must guess which is the most recent, by choosing the one that
    # has the largest SLiM time, doesn't appear in the parent list, or has
    # the lagest SLiM ID.
    nuc_inds = tables.mutations.metadata_vector(['mutation_list', 0, 'nucleotide'], dtype='int')
    num_stacked = np.array([len(m.metadata['mutation_list']) for m in ts.mutations()])
    for k in np.where(num_stacked > 1)[0]:
        mut = ts.mutation(k)
        if mut.parent == tskit.NULL:
            pids = []
        else:
            pids = ts.mutation(mut.parent).derived_state.split(",")
        x = [
            (
                md['slim_time'],
                i not in pids,
                int(i),
                j
            ) for j, (i, md) in
            enumerate(
                zip(mut.derived_state.split(","), mut.metadata['mutation_list'])
            )
        ]
        x.sort()
        j = x[-1][3]
        nuc_inds[k] = mut.metadata['mutation_list'][j]['nucleotide']
    if np.any(nuc_inds == -1):
        raise ValueError("All mutations must be nucleotide mutations.")
    da = np.array(NUCLEOTIDES)[nuc_inds]
    tables.mutations.packset_derived_state(da)
    k = tables.sites.position.astype('int')
    aa = np.frombuffer(ts.reference_sequence.data.encode('utf-8'), dtype='S1')[k]
    tables.sites.packset_ancestral_state(aa.tobytes().decode('utf-8'))

    return tables.tree_sequence()


def generate_nucleotides(ts, reference_sequence=None, keep=True, seed=None):
    """
    Returns a modified tree sequence in which mutations have been randomly assigned nucleotides
    and (optionally) a reference sequence has been randomly generated.

    If ``reference_sequence`` is a string of nucleotides (A, C, G, and T) of
    length equal to the sequence length, this is used for the reference
    sequence. If no reference sequence is given, the reference_sequence property
    of ts is used if present; if not then a sequence of independent and
    uniformly random nucleotides is generated.

    SLiM stores the nucleotide as an integer in the mutation metadata, with -1 meaning "not
    a nucleotide mutation". This method assigns nucleotides by stepping through
    each mutation and picking a random nucleotide uniformly out of the three
    possible nucleotides that differ from the parental state (i.e., the derived
    state of the parental mutation, or the ancestral state if the mutation has
    no parent). If ``keep=True`` (the default), the mutations that already have a 
    nucleotide (i.e., an integer 0-3 in metadata) will not be modified.

    Technical note: in the case of stacked mutations, the SLiM mutation that
    determines the nucleotide state of the (tskit) mutation is the one with the largest
    slim_time attribute. This method tries to assign nucleotides so that each mutation
    differs from the previous state, but this is not always possible in certain
    unlikely cases.

    :param tskit.TreeSequence ts: The tree sequence to transform.
    :param bool reference_sequence: A reference sequence, or None to use an existing reference,
        or to randomly generate one.
    :param bool keep: Whether to leave existing nucleotides in mutations that already have one.
    :param int seed: The random seed for generating new alleles.
    """
    rng = np.random.default_rng(seed=seed)
    if reference_sequence is None:
        if not ts.has_reference_sequence():
            reference_sequence = rng.choice(
                    np.array([65, 67, 71, 84], dtype=np.int8),
                    int(ts.sequence_length),
                    replace=True,
            ).tobytes().decode('ascii')
    else:
        if len(reference_sequence) != ts.sequence_length:
            raise ValueError("Reference sequence must have length equal to sequence_length.")
        if len([x for x in reference_sequence if x not in NUCLEOTIDES]) > 0:
            raise ValueError("Reference sequence must be a string of A, C, G, and T only.")

    tables = ts.dump_tables()
    if reference_sequence is not None:
        tables.reference_sequence.data = reference_sequence
    tables.mutations.clear()
    sets = [[k for k in range(4) if k != i] for i in range(4)]
    states = np.full((ts.num_mutations,), -1)
    k = tables.sites.position.astype('int')
    aa_list = np.searchsorted(
            NUCLEOTIDES,
            np.frombuffer(tables.reference_sequence.data.encode('utf-8'), dtype='S1')[k],
    )
    for site in ts.sites():
        aa = aa_list[site.id]
        muts = {}
        for mut in site.mutations:
            if mut.parent == tskit.NULL:
                pa = aa
                pds = []
            else:
                pa = states[mut.parent]
                pds = ts.mutation(mut.parent).derived_state.split(",")
            this_da = pa
            ml = mut.metadata
            max_time = -np.inf
            for i, md in zip(mut.derived_state.split(","), ml['mutation_list']):
                da = md['nucleotide']
                if da == -1 or not keep:
                    if i in muts:
                        da = muts[i]
                    else:
                        da = sets[pa][rng.integers(3)]
                    md['nucleotide'] = da
                muts[i] = da
                # the official nucleotide state is from the SLiM mutation with
                # the largest slim_time attribute that was not present in the parent
                # mutation
                if md["slim_time"] >= max_time and i not in pds:
                    this_da = da
                    max_time = md["slim_time"]
            states[mut.id] = this_da
            tables.mutations.append(mut.replace(metadata=ml))

    return tables.tree_sequence()


def individual_ages(ts):
    """
    Returns the ages of all individuals in the tree sequence, extracted
    from metadata. The result is a array of length equal to the number of
    individuals, with k-th entry equal to ``ts.individual(k).metadata["age"]``.

    :return: An array of ages of individuals.
    """
    if ts.metadata['SLiM']['model_type'] != "WF":
        ages = ts.tables.individuals.metadata_vector("age")
    else:
        ages = np.zeros(ts.num_individuals, dtype='int')
    return ages


def individuals_alive_at(ts, time, stage='late', remembered_stage=None,
                         population=None, samples_only=False):
    """
    Returns an array giving the IDs of all individuals that are known to be
    alive at the given time ago.  This is determined using their birth time
    ago (given by their `time` attribute) and, for nonWF models,
    their `age` attribute (which is equal to their age at the last time
    they were Remembered). See also {func}`.individual_ages_at`.

    In WF models, birth occurs after "early()", so that individuals are only
    alive during "late()" for the time step when they have age zero,
    while in nonWF models, birth occurs before "early()", so they are alive
    for both stages.
    
    In both WF and nonWF models, mortality occurs between
    "early()" and "late()", so that individuals are last alive during the
    "early()" stage of the time step of their final age, and if individuals
    are alive during "late()" they will also be alive during "first()" and
    "early()" of the next time step. This means it is important to know during
    which stage individuals were Remembered - for instance, if the call to
    sim.treeSeqRememberIndividuals() was made during "early()" of a given time step,
    then those individuals might not have survived until "late()" of that
    time step. Since SLiM does not record the stage at which individuals
    were Remembered, you can specify this by setting ``remembered_stages``:
    it should be the stage during which *all* calls to
    ``sim.treeSeqRememberIndividuals()``,  as well as to ``sim.treeSeqOutput()``,
    were made.

    Note also that in nonWF models, birth occurs between "first()" and
    "early()", so the possible parents in a given time step are those that are
    alive in "early()" and have age greater than zero, or, equivalently, are
    alive in "late()" during the previous time step.
    In WF models, birth occurs after "early()", so possible parents in a
    given time step are those that are alive during "early()" of that time
    step or are alive during "late()" of the previous time step.

    Since individuals may be created not during the usual 'birth' stage
    by `addSubPop( )`, and the stage at which they are created is not stored,
    results may be unreliable for the first-generation individuals.

    :param tskit.TreeSequence ts: A tree sequence.
    :param float time: The number of ticks (i.e., time steps) ago.
    :param str stage: The stage in the SLiM life cycle that we are inquiring
        about (either "first", "early" or "late"; defaults to "late").
    :param str remembered_stage: The stage in the SLiM life cycle
        during which individuals were Remembered (defaults to the stage the
        tree sequence was recorded at, stored in metadata).
    :param int population: If given, return only individuals in the
        population(s) with these population ID(s).
    :param bool samples_only: Whether to return only individuals who have at
        least one node marked as samples.
    """
    is_current_version(ts, _warn=True)
    if stage not in ("late", "early", "first"):
        raise ValueError(f"Unknown stage '{stage}': "
                          "should be either 'first', 'early' or 'late'.")

    if remembered_stage is None:
        remembered_stage = ts.metadata['SLiM']['stage']

    if remembered_stage not in ("late", "early", "first"):
        raise ValueError(f"Unknown remembered_stage '{remembered_stage}': "
                          "should be either 'first', 'early' or 'late'.")
    if remembered_stage != ts.metadata['SLiM']['stage']:
        warnings.warn(f"Provided remembered_stage '{remembered_stage}' does not"
                      " match the stage at which the tree sequence was saved"
                      f" ('{ts.metadata['SLiM']['stage']}'). This is not necessarily"
                      " an error, but mismatched stages will lead to inconsistencies:"
                      " make sure you know what you're doing.")

    # An individual's tskit time is the tskit counter at the time of their birth.
    # The tskit counter clicks once at the start of each reproduction bout.
    # In a WF model, individuals are alive for all stages with the same value
    # of the tskit counter. In a nonWF model, an individual of age `a` will be
    # alive for `a + 1` 'early' events, and `a` 'late' and 'first' events,
    # starting with the one having the same value of the tskit counter.
    # Also note that if remembered_stage is 'late', then they are recorded
    # before their age counter clicks, so they will actually have age one greater.
    # Since `time` is in "number of time steps ago",
    # we are inquiring about SLiM tick `n - time`, where `n` is
    # the value of `ts.metadata["SLiM"]["tick"]`. To convert this along with
    # our stage information into a tskit time `t`,
    # let x = 1 if the stage is 'first' or (is 'early' and WF)
    # and y = 1 if remembered stage is 'late' or (is 'early' and nonWF);
    # then t = time + x + y - 1 .
    is_wf = (ts.metadata['SLiM']['model_type'] == "WF")
    x = (stage == "first" or (stage == "early" and is_wf))
    y = (remembered_stage == "late" or (remembered_stage == "early" and not is_wf))
    t = time + x + y - 1
    birth_times = ts.individuals_time
    ages = individual_ages(ts)
    if is_wf:
        alive_bool = (birth_times == t)
    else:
        age_offset = (stage == "early") + (remembered_stage == "late")
        alive_bool = np.logical_and(
                birth_times >= t,
                birth_times - ages < t + age_offset
        )
    if population is not None:
        alive_bool &= np.isin(
                ts.individuals_population,
                population
        )
    if samples_only:
        nodes = ts.tables.nodes
        alive_bool &= (
                0 < np.bincount(1 + nodes.individual,
                                nodes.flags & tskit.NODE_IS_SAMPLE,
                                minlength=1 + ts.num_individuals)[1:]
        )
    return np.where(alive_bool)[0]


def individual_ages_at(ts, time, stage="late", remembered_stage="late"):
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
    does {func}`.individuals_alive_at`, but for consistency, non-nan ages
    will be 0 in "late" and 1 in "first" and "early".
    See {func}`.individuals_alive_at` for further discussion.

    :param tskit.TreeSequence ts: A tree sequence.
    :param float time: The reference time ago.
    :param str stage: The stage in the SLiM life cycle used to determine who
        is alive (either "early" or "late"; defaults to "late").
    :param str remembered_stage: The stage in the SLiM life cycle during which
        individuals were Remembered.
    """
    ages = np.repeat(np.nan, ts.num_individuals)
    alive = individuals_alive_at(
            ts,
            time,
            stage=stage,
            remembered_stage=remembered_stage
    )
    # to convert individuals_time to number of ticks ago we subtract (y - 1), so
    is_wf = (ts.metadata['SLiM']['model_type'] == "WF")
    y = (remembered_stage == "late" or (remembered_stage == "early" and not is_wf))
    t = time + y - 1
    ages[alive] = ts.individuals_time[alive] - t
    return ages


def slim_time(ts, time, stage="late"):
    """
    Converts the given "tskit times" (i.e., in units of time before the end
    of the simulation) to SLiM times (those recorded by SLiM, usually in units
    of ticks since the start of the simulation). Although the latter are
    always integers, these will not be if the provided times are not integers.

    When the tree sequence is written out, SLiM records the current
    current tick in the metadata:
    ``ts.metadata['SLiM']['tick']``. In most cases, the "SLiM time"
    referred to by a time ago in the tree sequence (i.e., the value that would
    be reported by community.tick within SLiM at the point in time thus
    referenced) can be obtained by subtracting that time ago from
    ``ts.metadata['SLiM']['tick']``. However, in WF models, birth
    happens between the “early()” and “late()” stages, so if the tree
    sequence was written out using sim.treeSeqOutput() during “early()” in
    a WF model, the tree sequence’s times measure time before the last set
    of individuals are born, i.e., before SLiM time step
    ``ts.metadata['SLiM']['tick'] - 1``. The same thing applies to
    the "first" stage for both WF and nonWF models.

    In some situations (e.g., mutations added during early() in WF models)
    this may not return what you expect. See :ref:`sec_metadata_converting_times`
    for more discussion.

    :param tskit.TreeSequence ts: A SLiM-compatible TreeSequence.
    :param numpy.ndarray time: An array of times to be converted.
    :param str stage: The stage of the SLiM life cycle that the SLiM time
        should be computed for.
    """
    is_current_version(ts, _warn=True)
    is_wf = (ts.metadata['SLiM']['model_type'] == "WF")
    remembered_stage = ts.metadata['SLiM']['stage']
    x = (stage == "first" or (stage == "early" and is_wf))
    y = (remembered_stage == "late" or (remembered_stage == "early" and not is_wf))
    slim_time = ts.metadata['SLiM']['tick'] - time + x + y - 1
    return slim_time


def _do_individual_parents_stuff(ts, return_parents=False):
    # Helper for has_individual_parents and individual_parents,
    # which share a lot of machinery.
    tables = ts.tables
    edges = tables.edges
    nodes = tables.nodes
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
    ind_times = ts.individuals_time
    ind_ages = individual_ages(ts)
    child_births = ind_times[edge_child_indiv[indiv_edges]]
    parent_births = ind_times[edge_parent_indiv[indiv_edges]]
    alive_edges = indiv_edges.copy()
    if ts.metadata['SLiM']['model_type'] == "WF":
        alive_edges[indiv_edges] = (child_births + 1 == parent_births)
    else:
        parent_deaths = parent_births - ind_ages[edge_parent_indiv[indiv_edges]]
        alive_edges[indiv_edges] = (child_births + 1 >= parent_deaths)
    edge_spans = edges.right - edges.left
    parental_span = np.bincount(edge_child_indiv[alive_edges],
            weights=edge_spans[alive_edges], minlength=ts.num_individuals)
    # we could also check for edges without individual parents terminating
    # in this individual, but this is unnecessary as the entire genome is
    # accounted for
    has_all_parents = (parental_span == 2 * ts.sequence_length)
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


def individual_parents(ts):
    '''
    Finds all parent-child relationships in the tree sequence (as far as we
    can tell). The output will be a two-column array with row [i,j]
    indicating that individual i is a parent of individual j.  See
    {func}`.has_individual_parents` for exactly which parents are returned.

    See {func}`.individuals_alive_at` for further discussion about how
    this is determined based on when the individuals were Remembered.

    :param tskit.TreeSequence ts: A :class:`tskit.TreeSequence`.
    :return: An array of individual IDs, with row [i, j] if individual i is
        a parent of individual j.
    '''
    return _do_individual_parents_stuff(ts, return_parents=True)


def has_individual_parents(ts):
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

    See {func}`.individuals_alive_at` for further discussion about how
    this is determined based on when the individuals were Remembered.

    :param tskit.TreeSequence ts: A :class:`tskit.TreeSequence`.
    :return: A boolean array of length equal to ``targets``.
    '''
    return _do_individual_parents_stuff(ts, return_parents=False)


def annotate(ts, **kwargs):
    '''
    Takes a tree sequence (as produced by msprime, for instance), and adds in the
    information necessary for SLiM to use it as an initial state, filling in
    mostly default values. Returns a :class:`tskit.TreeSequence`.

    :param tskit.TreeSequence ts: A :class:`tskit.TreeSequence`.
    :param str model_type: SLiM model type: either "WF" or "nonWF".
    :param int tick: What tick number in SLiM correponds to
        ``time=0`` in the tree sequence.
    :param int cycle: What cycle number in SLiM correponds to
        ``time=0`` in the tree sequence (default: equal to ``tick``).
    :param int stage: What stage in SLiM's cycle has the tree sequence
        been written out in (defaults to "early"; must be "early" or "late").
    :param str reference_sequence: A reference sequence of length
        equal to ts.sequence_length.
    :param bool annotate_mutations: Whether to replace mutation metadata
        with defaults. (If False, the mutation table is unchanged.)
    '''
    tables = ts.dump_tables()
    annotate_tables(tables, **kwargs)
    return tables.tree_sequence()


def annotate_tables(tables, model_type, tick, cycle=None, stage="early", reference_sequence=None,
        annotate_mutations=True):
    '''
    Does the work of :func:`annotate`, but modifies the tables in place: so,
    takes tables as produced by ``msprime``, and makes them look like the
    tables as output by SLiM. See :func:`annotate` for details.
    '''
    if stage not in ("early", "late"):
        raise ValueError(f"stage must be 'early' or 'late' (provided {stage})")
    if (type(tick) is not int) or (tick < 1):
        raise ValueError("tick must be an integer and at least 1.")
    # set_nodes must come before set_populations
    if model_type == "WF":
        default_ages = -1
    elif model_type == "nonWF":
        default_ages = 0
    else:
        raise ValueError("Model type must be 'WF' or 'nonWF'")
    if cycle is None:
        cycle = tick
    if not np.allclose(tables.sites.position, np.floor(tables.sites.position)):
        raise ValueError(
                "Site positions in this tree sequence are not at integer values, "
                "but must be for loading into SLiM: generate mutations with "
                "sim_mutations(..., discrete_genome=True), not simulate()."
        )
    top_metadata = default_slim_metadata('tree_sequence')['SLiM']
    top_metadata['model_type'] = model_type
    top_metadata['tick'] = tick
    top_metadata['cycle'] = cycle
    top_metadata['stage'] = stage
    set_tree_sequence_metadata(tables, **top_metadata)
    set_metadata_schemas(tables)
    _annotate_nodes_individuals(tables, age=default_ages)
    _annotate_populations(tables)
    if annotate_mutations:
        _annotate_sites_mutations(tables)
    if reference_sequence is not None:
        tables.reference_sequence.data = reference_sequence

def next_slim_mutation_id(ts):
    '''
    Returns the next SLiM mutation ID for this tree sequence. This is useful 
    because if you want to add more mutations to your SLiM tree sequence using 
    :func:`msprime.sim_mutations`, you may need to specify the parameter 
    `next_id` in your :class:`msprime.SLiMMutationModel` to be larger than any 
    existing mutation IDs. Setting `next_id` equal to the output of this 
    function will allow the mutated tree sequence to be read in by SLiM.
    To do this, recall that the "derived state" of SLiM's mutations are
    comma-separated strings of mutation IDs; this function just parses all derived
    states and returns one larger than the largest integer found. It will return an error
    if it encounters derived states that are not comma-separated strings of integers.
    '''
    max_id = 0
    for mut in ts.mutations():
        ds = mut.derived_state
        if len(ds) > 0:
            for d in ds.split(","):
                try:
                    max_id = max(max_id, int(d))
                except ValueError:
                    raise ValueError(f"The derived state of a mutation ({ds}) in the tree "
                                     "sequence is not a comma-separated list of values "
                                     "coercible to int. This is not a valid SLiM tree sequence.")
    return max_id + 1


def _annotate_nodes_individuals(tables, age):
    '''
    Adds to a TableCollection the information relevant to individuals required
    for SLiM to load in a tree sequence, that is found in Node and Individual
    tables.  This will replace the metadata in those tables. For this to work,
    assumes all sample nodes are in individuals.

    This is designed to make it easy to assign default values:
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
    ind_population = np.full(tables.individuals.num_rows, -1, dtype="int")
    ind_slim_id = np.full(tables.individuals.num_rows, 0, dtype='int')
    nid = 0
    node_metadata = []
    for j, n in enumerate(tables.nodes):
        if n.flags & tskit.NODE_IS_SAMPLE > 0:
            i = n.individual
            if i == tskit.NULL:
                raise ValueError("To annotate, samples must be in individuals.")
            else:
                ind_population[i] = n.population
                ind_slim_id[i] = 1
            md = default_slim_metadata("node")
            md["slim_id"] = nid
            nid += 1
        else:
            md = n.metadata
        node_metadata.append(md)

    nms = tables.nodes.metadata_schema
    tables.nodes.packset_metadata([
        nms.validate_and_encode_row(x)
        for x in node_metadata
    ])

    slim_ind = (ind_slim_id != 0)
    ind_slim_id = np.cumsum(ind_slim_id) - 1

    ind_metadata = []
    ind_flags = tables.individuals.flags
    for j, ind in enumerate(tables.individuals):
        if slim_ind[j]:
            md = default_slim_metadata("individual")
            md["pedigree_id"] = int(ind_slim_id[j])
            md["subpopulation"] = int(ind_population[j])
            md["age"] = age
            # no packset analogue for flags, so do it here
            # and we don't have to resize for flags,
            # so no big deal
            ind_flags[j] |= INDIVIDUAL_ALIVE
        else:
            md = ind.metdadata
        ind_metadata.append(md)
    tables.individuals.set_columns(
        flags=ind_flags,
        parents=tables.individuals.parents,
        parents_offset=tables.individuals.parents_offset,
    )
    ims = tables.individuals.metadata_schema
    tables.individuals.packset_metadata([
            ims.validate_and_encode_row(x)
            for x in ind_metadata
    ])
    tables.individuals.packset_location(
            [[0.0] * 3 if si else [] for si in slim_ind]
    )


def _annotate_populations(tables):
    '''
    Adds to a TableCollection the information about populations required for SLiM
    to load a tree sequence. This will replace anything already in the Population
    table for populations referenced by nodes alive at time zero.
    '''
    alive_ind = (tables.individuals.flags & INDIVIDUAL_ALIVE > 0)
    do_pops = np.unique(
        tables.nodes.population[
            np.logical_and(
                tables.nodes.individual >= 0,
                alive_ind[tables.nodes.individual]
            )
        ]
    )
    for j, p in enumerate(tables.populations):
        if j in do_pops:
            if "slim_id" not in p.metadata:
                md = default_slim_metadata("population")
                md["slim_id"] = j
                md["name"] = f"pop_{j}"
                tables.populations[j] = p.replace(metadata=md)


def _annotate_sites_mutations(tables):
    '''
    Adds to a TableCollection the information relevant to mutations required
    for SLiM to load in a tree sequence. This means adding to the metadata column
    of the Mutation table,  It will also
    - give SLiM IDs to each mutation
    - replace ancestral states with ""
    This will replace any information already in the metadata or derived state
    columns of the Mutation table. We set slim_time in metadata so that
    - tick = floor(tskit time) + slim_time
    '''
    if len(tables.mutations.metadata) > 0:
        warnings.warn(
                "The provided tree sequence already has some mutations with "
                "metadata; this metadata will be overwritten."
        )
    num_mutations = tables.mutations.num_rows
    default_mut = default_slim_metadata("mutation_list_entry")
    dsb, dso = tskit.pack_bytes([str(j).encode() for j in range(num_mutations)])
    slim_time = tables.metadata["SLiM"]["tick"] - np.floor(tables.mutations.time).astype("int")
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
