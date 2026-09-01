import functools
import warnings

import msprime
import numpy as np
import tskit

from pyslim import INDIVIDUAL_ALIVE, NODE_IS_VACANT_SAMPLE, NUCLEOTIDES

from .slim_metadata import (
    default_slim_metadata,
    is_current_version,
    set_metadata_schemas,
    set_tree_sequence_metadata,
)
from .slim_tree_sequence import mutation_metadata
from .util import unique_labels_by_group


def _mark_samples(tables, nodes):
    # Modifies tables in place.
    flags = tables.nodes.flags
    flags[nodes] |= np.uint32(tskit.NODE_IS_SAMPLE)
    tables.nodes.set_columns(
        flags=flags,
        time=tables.nodes.time,
        population=tables.nodes.population,
        individual=tables.nodes.individual,
        metadata=tables.nodes.metadata,
        metadata_offset=tables.nodes.metadata_offset,
    )


def _mark_not_samples(tables, nodes):
    # Modifies tables in place.
    flags = tables.nodes.flags
    flags[nodes] &= ~np.uint32(tskit.NODE_IS_SAMPLE)
    tables.nodes.set_columns(
        flags=flags,
        time=tables.nodes.time,
        population=tables.nodes.population,
        individual=tables.nodes.individual,
        metadata=tables.nodes.metadata,
        metadata_offset=tables.nodes.metadata_offset,
    )


def _chromosome_index(ts_metadata):
    """
    For a tree sequence produced by a multichromosome simulation, returns
    the index of the chromosome whose information is stored in this tree sequence
    in the list of all chromosomes. In a single-chromosome simulation,
    returns 0. This is extracted from
    ``ts.metadata['SLiM']['this_chromosome']``, and provides the index of this
    chromosome into `ts.metadata['SLiM']['chromosomes']``, if present.

    :param dict ts_metadata: The top-level metadata from a tree sequence
        or table collection.
    """
    if not (
        isinstance(ts_metadata, dict)
        and "SLiM" in ts_metadata
        and "this_chromosome" in ts_metadata["SLiM"]
    ):
        raise ValueError(
            "The tree sequence does not have the necessary "
            "information in top-level metadata."
        )
    k = ts_metadata["SLiM"]["this_chromosome"]["index"]
    return k


def _is_chrom_vacant(k, b):
    """
    :param int k: The index of the chromosome.
    :param int b: The list of ints containing the metadata.
    """
    # From the SLiM manual on is_vacant:
    # M bytes (uint8_t): a series of bytes comprising a bitfield of is_vacant
    # values, true (1) if this node represents a null haplosome for a given
    # chromosome, false (0) otherwise. For chromosomes with indices 0...N−1, the
    # chromosome with index k has its is_vacant bit in bit k%8 of byte k/8, where
    # byte 0 is the first byte in the series of bytes provided, and bit 0 is the
    # least-significant bit, the one with value 0x01 (hexadecimal 1). The number
    # of bytes present, M, is equal to (N+7)/8, the minimum number of bytes
    # necessary. The operators / and % here are integer divide (rounding down)
    # and integer modulo, respectively.
    assert len(b) >= (1 + k) / 8
    b = b[int(k / 8)]
    i = k % 8
    return (b >> i & 1) > 0


def has_vacant_samples(ts, ts_metadata=None):
    """
    Returns whether the tree sequence has vacant sample nodes.
    See :func:`remove_vacant`.

    :param tskit.TreeSequence ts: The tree sequence.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    out = False
    k = _chromosome_index(ts_metadata)
    for n in ts.samples():
        md = ts.node(n).metadata
        if md is not None:
            if _is_chrom_vacant(k, md["is_vacant"]):
                out = True
                break
    return out


def nodes_vacant(ts):
    """
    Evaluates which nodes in the tree sequence are vacant: returns a boolean
    vector whose k-th element is True if the k-th node is labelled as *vacant*
    in the node's metadata recorded by SLiM. A vacant node represents a blank
    placeholder in SLiM: either a "null haplosome" (used as placeholders for
    sex chromosomes and other chromosome types not of consistent ploidy in all
    individuals) or simply an unused node for haploid chromosome types. See
    :func:`remove_vacant`.


    :param tskit.TreeSequence ts: The tree sequence.
    :return boolean ndarray:
    """
    # not using chrom_index here because we expect people to call this on lots of nodes
    k = ts.metadata["SLiM"]["this_chromosome"]["index"]
    out = np.array(
        [
            node.metadata is not None and _is_chrom_vacant(k, node.metadata["is_vacant"])
            for node in ts.nodes()
        ],
        dtype="bool",
    )
    return out


def node_is_vacant(ts, node):
    """
    **DEPRECATED:** use :func:`.nodes_vacant` instead. This function requires
    top-level metadata access, which can be costly, so it is much better to do,
    for instance:

    .. code-block:: python

        vacant = nodes_vacant(ts)
        for node in ts.nodes():
            # instead of node_is_vacant(ts, node), use:
            vacant[node.id]

    Returns True if the node is labelled as *vacant* in the node's metadata
    recorded by SLiM. A vacant node represents a blank placeholder in SLiM:
    either a "null haplosome" (used as placeholders for sex chromosomes and other
    chromosome types not of consistent ploidy in all individuals) or simply an
    unused node for haploid chromosome types. See :func:`remove_vacant`.

    :param tskit.TreeSequence ts: The tree sequence.
    :param tskit.Node node: The node object.
    """
    warnings.warn(
        "The node_is_vacant method is deprecated: changes in SLiM v6 "
        " means that repeated use of this method will be unreasonably slow, "
        "so it will be removed in a future version of pyslim: "
        "obtain this information from pyslim.vacant_nodes( ) instead.",
        FutureWarning,
    )
    # not using chrom_index here because we expect people to call this on lots of nodes
    k = ts.metadata["SLiM"]["this_chromosome"]["index"]
    return node.metadata is not None and _is_chrom_vacant(k, node.metadata["is_vacant"])


def _record_vacant_tables(tables, ts_metadata):
    """
    Sets the NODE_IS_VACANT_SAMPLE flag for all vacant, sample nodes.
    See :func:`remove_vacant`.

    :param tskit.TableCollection tables: The table collection.
    """
    if np.any(tables.nodes.flags & NODE_IS_VACANT_SAMPLE):
        warnings.warn(
            "Some nodes are already flagged as vacant samples, and these "
            "flags are being overwritten; this may mean you've already run "
            "remove_vacant and so don't need to run it again."
        )
    k = _chromosome_index(ts_metadata)

    dn = tables.nodes.asdict()
    dn["flags"] &= ~NODE_IS_VACANT_SAMPLE
    samples = np.where(tables.nodes.flags & tskit.NODE_IS_SAMPLE > 0)[0]
    for s in samples:
        n = tables.nodes[s]
        if _is_chrom_vacant(k, n.metadata["is_vacant"]):
            dn["flags"][s] |= NODE_IS_VACANT_SAMPLE
    tables.nodes.set_columns(**dn)


def _remove_vacant_sample_flags(tables):
    """
    Set all NODE_IS_VACANT_SAMPLE flags off.
    """
    dn = tables.nodes.asdict()
    dn["flags"] &= ~NODE_IS_VACANT_SAMPLE
    tables.nodes.set_columns(**dn)


def remove_vacant(ts, ts_metadata=None):
    """
    Remove sample flags from all vacant nodes.

    In SLiM's internal state, there are two nodes per individual, even on
    chromosomes for which the individual is not diploid.  Thus, some sample
    nodes may be placeholders, not actually representing physical haplosomes,
    and are called "vacant" nodes. In the tree sequence these nodes have no
    ancestry, and thus have 'missing' data; however, their presence can
    cause problems with methods not designed for missing data.

    This method returns a copy of the tree sequence for which all vacant nodes
    have the sample flag removed; these nodes will thus not affect
    :func:`recapitate`, tree sequence statistics, etc. This also sets
    the :data:`NODE_IS_VACANT_SAMPLE` flag on these nodes,
    so they can be restored with :func:`restore_vacant`.  You probably don't
    want to run this method again on the output, since these flags will be
    overwritten and :func:`restore_vacant` will no longer work as expected.

    :param tskit.TreeSequence ts: The tree sequence.
    :return tskit.TreeSequence: A copy of the tree sequence with vacant nodes
        not marked as samples.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    tables = ts.dump_tables()
    remove_vacant_tables(tables, ts_metadata)
    return tables.tree_sequence()


def remove_vacant_tables(tables, ts_metadata=None):
    """
    Does the work of :func:`remove_vacant`, modifying ``tables`` in place.

    :param tskit.TableCollection tables: The tables underlying a tree sequence.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = tables.metadata
    _record_vacant_tables(tables, ts_metadata)
    is_vacant = np.where(tables.nodes.flags & NODE_IS_VACANT_SAMPLE > 0)[0]
    _mark_not_samples(tables, is_vacant)


def restore_vacant(ts, ts_metadata=None):
    """
    The inverse of :func:`remove_vacant`.

    This method returns a copy of the tree sequence for which all nodes
    with the :data:`NODE_IS_VACANT_SAMPLE` flag set have their sample flags
    set also and the :data:`NODE_IS_VACANT_SAMPLE` removed. If these nodes
    are not vacant, an error will be raised.

    :param tskit.TreeSequence ts: The tree sequence.
    :return tskit.TreeSequence: A copy of the tree sequence with vacant nodes
        marked as samples.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    tables = ts.dump_tables()
    restore_vacant_tables(tables, ts_metadata)
    return tables.tree_sequence()


def restore_vacant_tables(tables, ts_metadata=None):
    """
    Does the work of :func:`restore_vacant`, modifying ``tables`` in place.

    :param tskit.TableCollection tables: The tables underlying a tree sequence.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = tables.metadata
    is_vacant = np.where(tables.nodes.flags & NODE_IS_VACANT_SAMPLE > 0)[0]
    k = _chromosome_index(ts_metadata)
    for j in is_vacant:
        n = tables.nodes[j]
        if n.metadata is None:
            raise ValueError(
                "Something is wrong: node that is flagged as "
                "a vacant sample node has no metadata."
            )
        if not _is_chrom_vacant(k, n.metadata["is_vacant"]):
            raise ValueError(
                "Something is wrong: node that is flagged as "
                "a vacant sample node is not vacant."
            )
    _remove_vacant_sample_flags(tables)
    _mark_samples(tables, is_vacant)


def recapitate(ts, ancestral_Ne=None, *, keep_vacant=False, **kwargs):
    """
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
    :func:`msprime.sim_ancestry` are; this includes recombination rate, so
    that if neither ``recombination_rate`` or a ``recombination_map`` are
    provided, there will be *no* recombination.

    Sample flags from vacant nodes will be removed before recapitating:
    see :func:`remove_vacant`. To restore these, use ``keep_vacant=True``.
    You only need to do this if some individuals are not diploid and
    you will be loading the tree sequence back into SLiM.

    Precisely, this function will: (1) :func:`.remove_vacant` nodes, if any;
    (2) create a demography in which all populations split from
    a new population called "ancestral" of the desired size;
    (3) simulate with :func:`msprime.sim_ancestry`; and
    (4) :func:`.restore_vacant` nodes, if any.

    :param tskit.TreeSequence ts: The tree sequence to transform.
    :param float ancestral_Ne: If specified, then will simulate from a single
        ancestral population of this size. It is an error to specify this
        as well as ``demography``.
    :param bool keep_vacant: Whether to restore the sample flags on any
        vacant sample nodes. Default: False.
    :param dict kwargs: Any other arguments to :func:`msprime.sim_ancestry`.
    :return tskit.TreeSequence: A copy of the tree sequence with additional ancestral
        history added.
    """
    ts_metadata = ts.metadata
    is_current_version(ts_metadata, _warn=True)

    # we need to ask msprime to *not* simulate from any 'vacant' haplosomes;
    # which we do by marking these as not samples; note that `initial_state`
    # can take a TableCollection, not just a TreeSequence
    has_vacant = has_vacant_samples(ts, ts_metadata)
    if has_vacant:
        ts = remove_vacant(ts, ts_metadata)

    if ancestral_Ne is not None:
        if "demography" in kwargs:
            raise ValueError("You cannot specify both `demography` and `ancestral_Ne`.")
        # Since the tick can be set manually, the tick value in metadata has no
        # relationship to the time of the roots. Even in the case where the simulation
        # starts at tick=1, in various circumstances depending on the stage in
        # which the simulation was started and when the tree sequence was
        # written out, the time of the roots might be one or even two less than
        # the "tick" value. We have access to the stage the tree sequence was
        # written out, but not the time it was *started*. So, we'll just check
        # if we need to subtract one or two.  Consistency checking this
        # requires looping over all the trees, unfortunately, but it avoids
        # some common errors.
        root_times = set([ts.node(n).time for t in ts.trees() for n in t.roots])
        if len(root_times) > 1:
            message = (
                "Not all roots of the provided tree sequence are at the time expected "
                "by recapitate(). This could happen if you've simplified in "
                "python before recapitating (fix: don't simplify first, or "
                "pass keep_input_roots=True to simplify). "
                "It could also happen in other situations, e.g., "
                "you added new individuals without parents in SLiM "
                "during the course of the simulation with sim.addSubpop(), "
                "in which case you will probably need to recapitate with "
                "msprime.sim_ancestry(initial_state=ts, ...). "
                f"(Observed root times: {', '.join(map(str, list(root_times)[:5]))}"
                f"{', ...' if len(root_times) > 5 else ''})"
            )
            raise ValueError(message)
        recap_time = root_times.pop()
        demography = msprime.Demography.from_tree_sequence(ts)
        # must set pop sizes to >0 even though we merge immediately
        for pop in demography.populations:
            pop.initial_size = 1.0
        ancestral_name = "ancestral"
        derived_names = [pop.name for pop in demography.populations]
        while ancestral_name in derived_names:
            ancestral_name = ancestral_name + "_ancestral"
        demography.add_population(
            name=ancestral_name,
            description="ancestral population simulated by msprime",
            initial_size=ancestral_Ne,
        )
        # the split has to come slightly longer ago than slim's tick
        # since that's when all the linages are at, and otherwise the event
        # won't apply to them
        demography.add_population_split(
            np.nextafter(recap_time, 2 * recap_time),
            derived=derived_names,
            ancestral=ancestral_name,
        )
        kwargs["demography"] = demography

    recap = msprime.sim_ancestry(initial_state=ts, **kwargs)

    if has_vacant and keep_vacant:
        recap = restore_vacant(recap, ts_metadata)

    return recap


def add_mutation_metadata(ts, mutation_type=0, remove_unused=False):
    """
    Returns a new tree sequence with default information added to the top-level metadata
    for each mutation in the tree sequence that does not already have this information.
    To do this, mutations must be in SLiM format, as for instance produced by
    :class:`msprime.SLiMv6MutationModel`. Any information about SLiM mutations already
    in top-level metadata will remain unchanged.

    To do this, this method looks for all SLiM IDs that are found in the derived
    state of some mutation but are not represented in the top-level metadata
    (see :func:`.mutation_metadata`). This function then adds entries to that top-level
    metadata with default values (see :func:`.default_slim_metadata`),
    except that (a) the ``mutation_type`` can be specified;
    and (b) the ``slim_time`` is set using the ``tick`` value in top-level metadata
    and the ``time`` of the oldest tskit mutation in which the SLiM mutation occurs.

    :param tskit.TreeSequence ts: The tree sequence to transform.
    :param int mutation_type: The numeric ID of the mutation type in SLiM.
    :param bool remove_unused: Whether to also remove from metadata information about any
        mutations not seen in the derived states of the tree sequence.
    :return tskit.TreeSequence: A copy of the tree sequence with mutation information in
        metadata.
    """
    tables = ts.dump_tables()
    add_mutation_metadata_tables(
        tables, mutation_type=mutation_type, remove_unused=remove_unused
    )
    return tables.tree_sequence()


def add_mutation_metadata_tables(tables, mutation_type=0, remove_unused=False):
    """
    Does the work of :func:`.add_mutation_metadata`, modifying ``tables`` in place.

    :param tskit.TableCollection tables: The table collection to be modified.
    :param int mutation_type: The numeric ID of the mutation type in SLiM.
    :param bool remove_unused: Whether to also remove from metadata information about any
        mutations not seen in the derived states of the tree sequence.
    """
    ts_metadata = tables.metadata
    if (
        not isinstance(ts_metadata, dict)
        or "SLiM" not in ts_metadata
        or "SLiM_mutation_list" not in ts_metadata
    ):
        raise ValueError(
            "Top-level metadata schema is not correct: "
            "do you need to run pyslim.annotate()?"
        )
    num_traits = len(ts_metadata["SLiM"]["traits"])
    existing_muts = {x["mutation_id"] for x in ts_metadata["SLiM_mutation_list"]}
    mut_ids = [
        (int(j), mut.time)
        for mut in tables.mutations
        for j in mut.metadata["derived_states"]
    ]
    if len(mut_ids) > 0:
        mut_ids.sort()
        mut_ids = np.array(mut_ids, dtype="int")  # floors times
        # remove duplicate IDs, keeping the last (most recent)
        keep = np.full(len(mut_ids), True, dtype="bool")
        keep[np.where(np.diff(mut_ids[:, 0]) == 0)[0]] = False
        mut_ids = mut_ids[keep, :]
        mut_ids[:, 1] = slim_time(
            tables, mut_ids[:, 1], stage="late", ts_metadata=ts_metadata
        )
        # this assumes mutations were added in late(), which is what SLiM does
        ts_metadata["SLiM_mutation_list"].extend(
            [
                default_slim_metadata(
                    "mutation_list_entry",
                    num_traits=num_traits,
                    mutation_id=int(j),
                    mutation_type=mutation_type,
                    slim_time=int(t),
                )
                for j, t in mut_ids
                if j not in existing_muts
            ]
        )
    if remove_unused and len(mut_ids) < len(ts_metadata["SLiM_mutation_list"]):
        ts_metadata["SLiM_mutation_list"] = [
            x for x in ts_metadata["SLiM_mutation_list"] if x["mutation_id"] in mut_ids
        ]
    tables.metadata = ts_metadata


def convert_alleles(ts):
    """
    Returns a modified tree sequence in which alleles have been replaced by
    their corresponding nucleotides. For sites, SLiM-produced tree sequences
    have "" (the empty string) for the ancestral state at each site; this method
    will replace this with the corresponding nucleotide from the reference sequence.
    For mutations, SLiM records the 'derived state' as a SLiM mutation ID; this
    method will replace this with the nucleotide from the mutation's metadata.

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
    :return tskit.TreeSequence: A copy of the tree sequence with modified
        ancestral and derived states.
    """
    tables = ts.dump_tables()
    has_refseq = ts.has_reference_sequence() and len(ts.reference_sequence.data) > 0
    if not has_refseq:
        raise ValueError("Tree sequence must have a valid reference sequence.")
    # unfortunately, nucleotide mutations may be stacked (e.g., substitutions
    # will appear this way) and they don't appear in any particular order;
    # so we must guess which is the most recent, by choosing the one that
    # has the largest SLiM time, doesn't appear in the parent list, or has
    # the lagest SLiM ID.
    mut_metadata = mutation_metadata(ts)
    mut_ids = np.array([x["mutation_id"] for x in mut_metadata.values()], dtype="int")
    alleles = np.array([x["nucleotide"] for x in mut_metadata.values()], dtype="int")
    # mut_inds will map from tskit-mutations to slim-mutations
    mut_inds = np.array([mut.metadata["derived_states"][0] for mut in ts.mutations()])
    num_stacked = np.array(
        [len(mut.metadata["derived_states"]) for mut in ts.mutations()]
    )
    mut_inds[num_stacked > 0] = -1
    for k in np.where(num_stacked > 0)[0]:
        mut = ts.mutation(k)
        if mut.parent == tskit.NULL:
            pids = []
        else:
            pids = ts.mutation(mut.parent).metadata["derived_states"]
        x = [
            (mut_metadata[i]["slim_time"], i not in pids, i, j)
            for j, i in enumerate(mut.metadata["derived_states"])
        ]
        x.sort()
        mut_inds[k] = x[-1][2]
    assert np.all(mut_inds >= 0), "This should not occur: please file a bug report."
    nuc_inds = alleles[np.searchsorted(mut_ids, mut_inds)]
    if np.any(nuc_inds == -1):
        raise ValueError("All mutations must be nucleotide mutations.")
    da = np.array(NUCLEOTIDES)[nuc_inds]
    tables.mutations.packset_derived_state(da)
    k = tables.sites.position.astype("int")
    aa = np.frombuffer(ts.reference_sequence.data.encode("utf-8"), dtype="S1")[k]
    tables.sites.packset_ancestral_state(aa.tobytes().decode("utf-8"))
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
    :return tskit.TreeSequence: A copy of the tree sequence with nucleotides.
    """
    rng = np.random.default_rng(seed=seed)
    if reference_sequence is None:
        if not ts.has_reference_sequence():
            reference_sequence = (
                rng.choice(
                    np.array([65, 67, 71, 84], dtype=np.int8),
                    int(ts.sequence_length),
                    replace=True,
                )
                .tobytes()
                .decode("ascii")
            )
    else:
        if len(reference_sequence) != ts.sequence_length:
            raise ValueError(
                "Reference sequence must have length equal to sequence_length."
            )
        if len([x for x in reference_sequence if x not in NUCLEOTIDES]) > 0:
            raise ValueError(
                "Reference sequence must be a string of A, C, G, and T only."
            )
    ts_metadata = ts.metadata
    mut_info = mutation_metadata(ts, ts_metadata=ts_metadata)
    tables = ts.dump_tables()
    if reference_sequence is not None:
        tables.reference_sequence.data = reference_sequence
    sets = [[k for k in range(4) if k != i] for i in range(4)]
    states = np.full((ts.num_mutations,), -1)
    k = tables.sites.position.astype("int")
    aa_list = np.searchsorted(
        NUCLEOTIDES,
        np.frombuffer(tables.reference_sequence.data.encode("utf-8"), dtype="S1")[k],
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
                pds = ts.mutation(mut.parent).metadata["derived_states"]
            this_da = pa
            max_time = -np.inf
            for i in mut.metadata["derived_states"]:
                md = mut_info[i]
                da = md["nucleotide"]
                if da == -1 or not keep:
                    if i in muts:
                        da = muts[i]
                    else:
                        da = sets[pa][rng.integers(3)]
                    md["nucleotide"] = da
                muts[i] = da
                # the official nucleotide state is from the SLiM mutation with
                # the largest slim_time attribute that was not present in the parent
                # mutation
                if md["slim_time"] >= max_time and i not in pds:
                    this_da = da
                    max_time = md["slim_time"]
            states[mut.id] = this_da
    ts_metadata["SLiM"]["nucleotide_based"] = True
    ts_metadata["SLiM_mutation_list"] = list(mut_info.values())
    tables.metadata = ts_metadata
    return tables.tree_sequence()


def individual_ages(ts, ts_metadata=None):
    """
    Returns the ages of all individuals in the tree sequence, as extracted
    from metadata. The result is an array of length equal to the number of
    individuals, with k-th entry equal to ``ts.individual(k).metadata["age"]``.

    These are the ages of the indivdiuals when they were recorded in the tree sequence:
    either when the tree sequence was saved (if they are alive) or when they were
    last Remembered. See also :func:`.individual_ages_at`.

    :param tskit.TreeSequence ts: The tree sequence.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    :return: An array of ages of individuals.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    if ts_metadata["SLiM"]["model_type"] != "WF":
        ages = ts.tables.individuals.metadata_vector("age")
    else:
        ages = np.zeros(ts.num_individuals, dtype="int")
    return ages


def individuals_alive_at(
    ts,
    time,
    stage="late",
    remembered_stage=None,
    population=None,
    samples_only=False,
    ts_metadata=None,
):
    """
    Returns an array giving the IDs of all individuals that are known to be
    alive at the given time ago.  This is determined using their birth time
    ago (given by their `time` attribute) and, for nonWF models,
    their `age` attribute (which is equal to their age at the last time
    they were Remembered). See also :func:`individual_ages_at`.

    In WF models, birth occurs after "early()", so that individuals are only
    alive during "late()" for the time step when they have age zero,
    while in nonWF models, birth occurs before "early()", so they are alive
    for both stages.

    In both the WF and nonWF life cycles, mortality occurs between
    "early()" and "late()", so that individuals are last alive during the
    "early()" stage of the time step of their final age, and if individuals
    are alive during "late()" they will also be alive during "first()" and
    "early()" of the next time step (unless the user calls `killIndividuals()`).
    This means it is important to know during
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
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if stage not in ("late", "early", "first"):
        raise ValueError(
            f"Unknown stage '{stage}': should be either 'first', 'early' or 'late'."
        )
    if ts_metadata is None:
        ts_metadata = ts.metadata
    is_current_version(ts_metadata, _warn=True)
    if remembered_stage is None:
        remembered_stage = ts_metadata["SLiM"]["stage"]

    if remembered_stage not in ("late", "early", "first"):
        raise ValueError(
            f"Unknown remembered_stage '{remembered_stage}': "
            "should be either 'first', 'early' or 'late'."
        )
    if remembered_stage != ts_metadata["SLiM"]["stage"]:
        warnings.warn(
            f"Provided remembered_stage '{remembered_stage}' does not"
            " match the stage at which the tree sequence was saved"
            f" ('{ts_metadata['SLiM']['stage']}'). This is not necessarily"
            " an error, but mismatched stages will lead to inconsistencies:"
            " make sure you know what you're doing."
        )

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
    is_wf = ts_metadata["SLiM"]["model_type"] == "WF"
    x = stage == "first" or (stage == "early" and is_wf)
    y = remembered_stage == "late" or (remembered_stage == "early" and not is_wf)
    t = time + x + y - 1
    birth_times = ts.individuals_time
    ages = individual_ages(ts, ts_metadata)
    if is_wf:
        alive_bool = birth_times == t
    else:
        age_offset = (stage == "early") + (remembered_stage == "late")
        alive_bool = np.logical_and(
            birth_times >= t, birth_times - ages < t + age_offset
        )
    if population is not None:
        alive_bool &= np.isin(ts.individuals_population, population)
    if samples_only:
        nodes = ts.tables.nodes
        alive_bool &= (
            0
            < np.bincount(
                1 + nodes.individual,
                nodes.flags & tskit.NODE_IS_SAMPLE,
                minlength=1 + ts.num_individuals,
            )[1:]
        )
    return np.where(alive_bool)[0]


def individual_ages_at(
    ts, time, stage="late", remembered_stage="late", ts_metadata=None
):
    """
    Returns the `ages` of each individual at the corresponding time ago,
    which will be ``nan`` if the individual is either not born yet or dead.
    This is computed as the time ago the individual was born (found by the
    `time` associated with the the individual's nodes) minus the `time`
    argument; while "death" is inferred from the individual's ``age``,
    recorded in metadata (see :func:`.individual_ages`). These values should
    be the same as what would be shown in SLiM during the corresponding time
    step and stage.

    Since age increments at the end of each time step,
    the age is the number of time steps ends the individual has lived
    through, so if they were born in time step `time`, then their age
    will be zero.

    In a WF model, this method does not provide any more information than
    does :func:`individuals_alive_at`, but for consistency, non-nan ages
    will be 0 in "late" and 1 in "first" and "early".
    See :func:`individuals_alive_at` for further discussion.

    :param tskit.TreeSequence ts: A tree sequence.
    :param float time: The reference time ago.
    :param str stage: The stage in the SLiM life cycle used to determine who
        is alive (either "early" or "late"; defaults to "late").
    :param str remembered_stage: The stage in the SLiM life cycle during which
        individuals were Remembered.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``. If
        this does not match the actual top-level metadata, incorrect values may result.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    ages = np.repeat(np.nan, ts.num_individuals)
    alive = individuals_alive_at(
        ts, time, stage=stage, remembered_stage=remembered_stage, ts_metadata=ts_metadata
    )
    # to convert individuals_time to number of ticks ago we subtract (y - 1), so
    is_wf = ts_metadata["SLiM"]["model_type"] == "WF"
    y = remembered_stage == "late" or (remembered_stage == "early" and not is_wf)
    t = time + y - 1
    ages[alive] = ts.individuals_time[alive] - t
    return ages


def slim_time(ts, time, stage="late", ts_metadata=None):
    """
    Converts the given "tskit times" (i.e., in units of time before the end
    of the simulation) to SLiM times (those recorded by SLiM, usually in units
    of ticks since the start of the simulation). Although times in SLiM are
    always integers, the returned values will not be integers if the values
    in `time` are not.

    When the tree sequence is written out, SLiM records the current
    current tick in the metadata:
    ``ts.metadata['SLiM']['tick']``. In most cases, the "SLiM time"
    referred to by a time ago in the tree sequence (i.e., the value that would
    be reported by community.tick within SLiM at the point in time thus
    referenced) can be obtained by subtracting that time ago from
    ``ts.metadata['SLiM']['tick']``. However, in WF models, birth
    happens between the "early()" and "late()" stages, so if the tree
    sequence was written out using sim.treeSeqOutput() during "early()" in
    a WF model, the tree sequence's times measure time before the last set
    of individuals are born, i.e., before SLiM time step
    ``ts.metadata['SLiM']['tick'] - 1``. The same thing applies to
    the "first" stage for both WF and nonWF models.

    In some situations (e.g., mutations added during early() in WF models)
    this may not return what you expect. See :ref:`sec_metadata_converting_times`
    for more discussion.

    This method accesses top-level metadata, which may be a costly operation,
    so if this method will be called many times, it is recommended to
    extract this to a variable (e.g., ``ts_metadata = ts.metadata``) and pass it
    to this method (as ``ts_metadata``). However, beware: if ``ts_metadata``
    is not in sync with the actual top-level metadata, incorrect values may result.

    :param tskit.TreeSequence ts: A SLiM-compatible TreeSequence.
    :param numpy.ndarray time: An array of times to be converted.
    :param str stage: The stage of the SLiM life cycle that the SLiM time
        should be computed for.
    :param dict ts_metadata: Optionally, the top-level metadata for ``ts``.
    """
    if ts_metadata is None:
        ts_metadata = ts.metadata
    is_current_version(ts_metadata, _warn=True)
    is_wf = ts_metadata["SLiM"]["model_type"] == "WF"
    remembered_stage = ts_metadata["SLiM"]["stage"]
    x = stage == "first" or (stage == "early" and is_wf)
    y = remembered_stage == "late" or (remembered_stage == "early" and not is_wf)
    slim_time = ts_metadata["SLiM"]["tick"] - time + x + y - 1
    return slim_time


def _shift_times(ts, dt):
    """
    Add `dt` to the times ago of all nodes and mutations.
    (Note: migrations not currently supported by simplify,
    so not used here.)
    """
    tables = ts.dump_tables()
    if dt != 0:
        d = tables.nodes.asdict()
        d["time"] += dt
        tables.nodes.set_columns(**d)
        d = tables.mutations.asdict()
        d["time"] += dt
        tables.mutations.set_columns(**d)
    return tables


def set_slim_state(ts, time=0, individuals=None):
    """
    Returns a new tree sequence for which
    the information stored in metadata has been changed so that when loaded into
    SLiM, the current time will be ``time`` units ago (i.e., at tskit time ``time``)
    and the alive individuals will be ``individuals``. The time in SLiM
    (i.e., the value of the tick counter) will also be ``time`` units earlier.

    This is useful, for instance, if you Remember a set of individuals at some
    point partway through a SLiM simulation, and would like to load *those*
    individuals into SLiM, rather than the final generation. (These could be
    used to start a new simulation, for instance.)

    Appropriate individuals might be found, for instance, with
    ``pyslim.individuals_alive_at(ts, time)``.

    To do this, the "tick" in top-level metadata is changed; individual flags
    are reset so that only ``individuals`` have the :data:`INDIVIDUAL_ALIVE`
    flag set; and tskit times have ``time`` subtracted from them (so they measure
    time ago relative to the new tick).

    As a result, some individuals in the tree sequence may have negative times,
    i.e., have lived "in the future".  Since these individuals will not be
    "alive", they will be ignored by SLiM, even when their birth times arrive,
    so that any future history will also be present, unchanged, in the tree
    sequence that results from additional simulation. If additional simulation
    is performed then contradictions may arise: for instance, if one of the
    ``individuals`` had offspring for an additional ten ticks past ``time`` in
    the original simulation, but in the new simulation they die immediately,
    then the resulting tree sequence will still be valid but the history as
    recorded by SLiM in metadata will not make sense. To avoid such issues, one
    can in the original SLiM script Remember and then immediately kill the
    individuals to be later used in this way.

    This also subtracts ``time`` from the value of "cycle" in top-level
    metadata; if this is not desired (or the value in "tick" needs to be
    adjusted), change the metadata directly.

    :param tskit.TreeSequence ts: A SLiM-compatible TreeSequence.
    :param int time: The number of time units (ticks) into the past to shift
        times.  (Default: zero; can be positive or negative.)
    :param numpy.ndarray individuals: An array of the tskit IDs of the individuals
        that should be marked as alive (all others will be not alive).
        (Default: leave unchanged.)
    :return tskit.TreeSequence: A copy of the tree sequence, modified.
    """
    tables = _shift_times(ts, -time)
    if time != 0:
        md = tables.metadata
        md["SLiM"]["tick"] -= time
        md["SLiM"]["cycle"] -= time
        tables.metadata = md
    if individuals is not None:
        d = tables.individuals.asdict()
        d["flags"] &= ~INDIVIDUAL_ALIVE
        d["flags"][individuals] |= INDIVIDUAL_ALIVE
        tables.individuals.set_columns(**d)
    return tables.tree_sequence()


def _do_individual_parents_stuff(ts, return_parents=False):
    warnings.warn(
        "The individual_parents( ) and has_individual_parents( ) methods are "
        "no longer needed and will be removed in a future version of pyslim: "
        "obtain this information from the `parents' property of individuals "
        "instead.",
        FutureWarning,
    )

    # Helper for has_individual_parents and individual_parents,
    # which share a lot of machinery.
    tables = ts.tables
    edges = tables.edges
    nodes = tables.nodes
    edge_parent_indiv = nodes.individual[edges.parent]
    edge_child_indiv = nodes.individual[edges.child]
    # nodes whose parent nodes are all in the same individual
    unique_parent_nodes = unique_labels_by_group(
        edges.child, edge_parent_indiv, minlength=nodes.num_rows
    )
    unique_parent_edges = unique_parent_nodes[edges.child]
    # edges describing relationships between individuals
    indiv_edges = np.logical_and(
        np.logical_and(edge_parent_indiv != tskit.NULL, edge_child_indiv != tskit.NULL),
        unique_parent_edges,
    )
    # individual edges where the parent was alive during "late"
    # of the time step before the child is born
    ind_times = ts.individuals_time
    ind_ages = individual_ages(ts)
    child_births = ind_times[edge_child_indiv[indiv_edges]]
    parent_births = ind_times[edge_parent_indiv[indiv_edges]]
    alive_edges = indiv_edges.copy()
    if ts.metadata["SLiM"]["model_type"] == "WF":
        alive_edges[indiv_edges] = child_births + 1 == parent_births
    else:
        parent_deaths = parent_births - ind_ages[edge_parent_indiv[indiv_edges]]
        alive_edges[indiv_edges] = child_births + 1 >= parent_deaths
    edge_spans = edges.right - edges.left
    parental_span = np.bincount(
        edge_child_indiv[alive_edges],
        weights=edge_spans[alive_edges],
        minlength=ts.num_individuals,
    )
    # we could also check for edges without individual parents terminating
    # in this individual, but this is unnecessary as the entire genome is
    # accounted for
    has_all_parents = parental_span == 2 * ts.sequence_length
    if return_parents:
        full_parent_edges = np.logical_and(
            alive_edges, has_all_parents[edge_child_indiv]
        )
        parents = np.unique(
            np.column_stack(
                [
                    edge_parent_indiv[full_parent_edges],
                    edge_child_indiv[full_parent_edges],
                ]
            ),
            axis=0,
        )
        return parents
    else:
        return has_all_parents


def individual_parents(ts):
    """
    **DEPRECATED:** now SLiM records `parents` directly in the individual
    table (see for instance `ind.parents`).

    Finds all parent-child relationships in the tree sequence (as far as we
    can tell). The output will be a two-column array with row [i,j]
    indicating that individual i is a parent of individual j.  See
    :func:`has_individual_parents` for exactly which parents are returned.

    See :func:`individuals_alive_at` for further discussion about how
    this is determined based on when the individuals were Remembered.

    :param tskit.TreeSequence ts: A :class:`tskit.TreeSequence`.
    :return: An array of individual IDs, with row [i, j] if individual i is
        a parent of individual j.
    """
    return _do_individual_parents_stuff(ts, return_parents=True)


def has_individual_parents(ts):
    """
    **DEPRECATED:** now SLiM records `parents` directly in the individual
    table (see for instance `ind.parents`).

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

    See :func:`individuals_alive_at` for further discussion about how
    this is determined based on when the individuals were Remembered.

    :param tskit.TreeSequence ts: A :class:`tskit.TreeSequence`.
    :return: A boolean array of length equal to ``targets``.
    """
    return _do_individual_parents_stuff(ts, return_parents=False)


def annotate(
    ts,
    model_type,
    tick,
    cycle=None,
    stage="early",
    reference_sequence=None,
    annotate_mutations=True,
):
    """
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
        with defaults. (If False, information about mutations is unchanged.)
    """
    tables = ts.dump_tables()
    annotate_tables(
        tables,
        model_type=model_type,
        tick=tick,
        cycle=cycle,
        stage=stage,
        reference_sequence=reference_sequence,
        annotate_mutations=annotate_mutations,
    )
    return tables.tree_sequence()


def annotate_tables(
    tables,
    model_type,
    tick,
    cycle=None,
    stage="early",
    reference_sequence=None,
    annotate_mutations=True,
):
    """
    Does the work of :func:`annotate`, but modifies the tables in place: so,
    takes tables as produced by ``msprime``, and makes them look like the
    tables as output by SLiM. See :func:`annotate` for details.
    """
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
    top_metadata = default_slim_metadata("tree_sequence")["SLiM"]
    top_metadata["model_type"] = model_type
    top_metadata["tick"] = tick
    top_metadata["cycle"] = cycle
    top_metadata["stage"] = stage
    md = tables.metadata
    if isinstance(md, dict) and "SLiM_mutation_list" in md:
        top_metadata["SLiM_mutation_list"] = md["SLiM_mutation_list"]
    ts_metadata = set_tree_sequence_metadata(tables, **top_metadata)
    set_metadata_schemas(tables)
    _annotate_nodes_individuals(tables, age=default_ages)
    _annotate_populations(tables)
    if annotate_mutations:
        _annotate_sites_mutations(tables, ts_metadata=ts_metadata)
    if reference_sequence is not None:
        tables.reference_sequence.data = reference_sequence


def next_slim_mutation_id(ts):
    """
    Returns the next unused SLiM mutation ID for this tree sequence. This is useful
    because if you want to add more mutations to your SLiM tree sequence using
    :func:`msprime.sim_mutations`, you may need to specify the parameter
    `next_id` in your :class:`msprime.SLiMv6MutationModel` to be larger than any
    existing mutation IDs. Setting `next_id` equal to the output of this
    function will allow the mutated tree sequence to be read in by SLiM.
    To do this, recall that the "derived state" of SLiM's mutations are
    comma-separated strings of mutation IDs; this function just parses all derived
    states and returns one larger than the largest integer found. It will return an error
    if it encounters derived states that are not comma-separated strings of integers.
    """
    max_id = -1
    if ts.num_mutations > 0:
        try:
            max_id = functools.reduce(
                max,
                (x for mut in ts.mutations() for x in mut.metadata["derived_states"]),
                -1,
            )
        except TypeError:
            raise ValueError(
                "The mutation metadata in this tree sequence is not in proper format."
            )
    return max_id + 1


def _annotate_nodes_individuals(tables, age):
    """
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
    - (node_is_vacant) genomes to be non-null
    - (ind_flags) INDIVIDUAL_ALIVE

    If you have other situations, like non-alive "remembered" individuals, you
    will need to edit the tables by hand, afterwards.
    """
    if len(tables.nodes.metadata) > 0:
        warnings.warn(
            "The provided tree sequence already has some nodes with "
            "metadata; this metadata will be overwritten."
        )
    if len(tables.individuals.metadata) > 0:
        warnings.warn(
            "The provided tree sequence already has some individuals with "
            "metadata; this metadata will be overwritten."
        )
    ind_population = np.full(tables.individuals.num_rows, -1, dtype="int")
    ind_slim_id = np.full(tables.individuals.num_rows, 0, dtype="int")
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
            md = None
        node_metadata.append(md)
    nms = tables.nodes.metadata_schema
    tables.nodes.packset_metadata(
        [nms.validate_and_encode_row(x) for x in node_metadata]
    )

    slim_ind = ind_slim_id != 0
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
            md = None
        ind_metadata.append(md)
    tables.individuals.set_columns(
        flags=ind_flags,
        parents=tables.individuals.parents,
        parents_offset=tables.individuals.parents_offset,
    )
    ims = tables.individuals.metadata_schema
    tables.individuals.packset_metadata(
        [ims.validate_and_encode_row(x) for x in ind_metadata]
    )
    tables.individuals.packset_location([[0.0] * 3 if si else [] for si in slim_ind])


def _annotate_populations(tables):
    """
    Adds to a TableCollection the information about populations required for SLiM
    to load a tree sequence. This will replace anything already in the Population
    table for populations referenced by nodes alive at time zero.
    """
    alive_ind = tables.individuals.flags & INDIVIDUAL_ALIVE > 0
    do_pops = np.unique(
        tables.nodes.population[
            np.logical_and(
                tables.nodes.individual >= 0, alive_ind[tables.nodes.individual]
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


def _annotate_sites_mutations(tables, ts_metadata):
    """
    Adds to a TableCollection the information relevant to mutations required
    for SLiM to load in a tree sequence. This means adding metadata to the
    SLiM_mutation_list in top-level metadata. It will also:
    - give SLiM IDs to each mutation
    - replace ancestral states with ""
    This will replace any information already in the metadata, and the derived
    state columns of the Mutation table. We set slim_time in metadata so that
    - tick = floor(tskit time) + slim_time
    """
    if (
        isinstance(ts_metadata, dict)
        and "SLiM_mutation_list" in ts_metadata
        and len(ts_metadata["SLiM_mutation_list"]) > 0
    ):
        warnings.warn(
            "The provided tree sequence already has top-level mutation "
            "metadata; this metadata will be overwritten."
        )
    num_mutations = tables.mutations.num_rows
    default_mut = default_slim_metadata("mutation_list_entry")
    slim_time = ts_metadata["SLiM"]["tick"] - np.floor(tables.mutations.time).astype(
        "int"
    )
    mutation_list = [
        {
            "mutation_id": j,
            "mutation_type": default_mut["mutation_type"],
            "per_trait": default_mut["per_trait"],
            "subpopulation": default_mut["subpopulation"],
            "slim_time": int(st),
            "nucleotide": default_mut["nucleotide"],
            "padding": None,
        }
        for j, st in enumerate(slim_time)
    ]
    ts_metadata["SLiM_mutation_list"] = mutation_list
    tables.metadata = ts_metadata
    mut_ids = np.arange(num_mutations, dtype="int64")
    dsb, dso = tskit.pack_bytes([str(j).encode() for j in mut_ids])
    mdb, mdo = tskit.pack_bytes([np.frombuffer(j, dtype="int8") for j in mut_ids])
    tables.mutations.set_columns(
        site=tables.mutations.site,
        node=tables.mutations.node,
        time=tables.mutations.time,
        derived_state=dsb,
        derived_state_offset=dso,
        parent=tables.mutations.parent,
        metadata=mdb,
        metadata_offset=mdo,
    )
    tables.sites.set_columns(
        position=tables.sites.position,
        ancestral_state=np.array([], dtype="int8"),
        ancestral_state_offset=np.zeros(tables.sites.num_rows + 1, dtype="uint32"),
    )
