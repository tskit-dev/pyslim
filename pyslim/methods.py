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
               **kwargs):
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

    :param TreeSequence ts: The tree sequence to transform.
    :param float ancestral_Ne: If specified, then will simulate from a single
        ancestral population of this size. It is an error to specify this
        as well as ``demography``.
    :param dict kwargs: Any other arguments to :func:`msprime.sim_ancestry`.
    '''
    if ancestral_Ne is not None:
        if "demography" in kwargs:
            raise ValueError("You cannot specify both `demography` and `ancestral_Ne`.")
        demography = msprime.Demography.from_tree_sequence(ts)
        # must set pop sizes to >0 even though we merge immediately
        for pop in demography.populations:
            pop.initial_size=1.0
        ancestral_name = "ancestral"
        derived_names = [pop.name for pop in demography.populations]
        while ancestral_name in derived_names:
            ancestral_name = (ancestral_name + "_ancestral")
        ancestral_metadata = {}
        ancestral_metadata['slim_id'] = ts.num_populations
        demography.add_population(
                name=ancestral_name,
                description="ancestral population simulated by msprime",
                initial_size=ancestral_Ne,
                extra_metadata=ancestral_metadata,
        )
        # the split has to come slightly longer ago than slim_generation,
        # since that's when all the linages are at, and otherwise the event
        # won't apply to them
        demography.add_population_split(
                np.nextafter(
                    ts.metadata['SLiM']['generation'],
                    2 * ts.metadata['SLiM']['generation'],
                ),
                derived=derived_names,
                ancestral=ancestral_name,
        )
        kwargs["demography"] = demography

    recap = msprime.sim_ancestry(
                initial_state = ts,
                **kwargs)

    return SlimTreeSequence(recap)


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

    :param TreeSequence ts: The tree sequence to transform.
    """
    tables = ts.dump_tables()
    has_refseq = (
            ts.has_reference_sequence()
            and len(ts.reference_sequence.data) > 0
    )
    # unfortunately, nucleotide mutations may be stacked (e.g., substitutions
    # will appear this way) and they don't appear in any particular order;
    # so we must guess which is the most recent, by choosing the one that
    # has the largest SLiM time, doesn't appear in the parent list, or has
    # the lagest SLiM ID.
    slim_muts = np.repeat(0, ts.num_mutations)
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
        slim_muts[k] = x[-1][3]
    # gah that was terrible, ok back to it
    nuc_inds = np.array(
            [m.metadata['mutation_list'][0]['nucleotide']
                for k, m in zip(slim_muts, ts.mutations())],
            dtype='int',
    )
    if np.any(nuc_inds == -1):
        raise ValueError("All mutations must be nucleotide mutations.")
    if not has_refseq:
        raise ValueError("Tree sequence must have a valid reference sequence.")
    aa = [ts.reference_sequence.data[k] for k in tables.sites.position.astype('int')]
    da = np.array(NUCLEOTIDES)[nuc_inds]
    tables.sites.packset_ancestral_state(aa)
    tables.mutations.packset_derived_state(da)

    out = SlimTreeSequence(tables.tree_sequence())
    return out


def generate_nucleotides(ts, reference_sequence=None, keep=True, seed=None):
    """
    Returns a modified tree sequence in which mutations have been randomly assigned nucleotides
    and (optionally) a reference sequence has been randomly generated.

    If ``reference_sequence`` is a string of nucleotides (A, C, G, and T) of
    length equal to the sequence length, this is used for the reference
    sequence. Otherwise (the default), a reference sequence of independent and
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

    :param TreeSequence ts: The tree sequence to transform.
    :param bool reference_sequence: A reference sequence, or None to randomly generate one.
    :param bool keep: Whether to leave existing nucleotides in mutations that already have one.
    :param int seed: The random seed for generating new alleles.
    """
    rng = np.random.default_rng(seed=seed)
    if reference_sequence is None:
        reference_sequence = rng.choice(
                np.array([65, 67, 71, 84], dtype=np.int8),
                int(ts.sequence_length),
                replace=True,
        ).tobytes().decode('ascii')

    if len(reference_sequence) != ts.sequence_length:
        raise ValueError("Reference sequence must have length equal to sequence_length.")
    if len([x for x in reference_sequence if x not in NUCLEOTIDES]) > 0:
        raise ValueError("Reference sequence must be a string of A, C, G, and T only.")

    tables = ts.dump_tables()
    tables.reference_sequence.data = reference_sequence
    tables.mutations.clear()
    sets = [[k for k in range(4) if k != i] for i in range(4)]
    states = np.full((ts.num_mutations,), -1)
    for site in ts.sites():
        aa = NUCLEOTIDES.index(reference_sequence[int(site.position)])
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

    return SlimTreeSequence(
            tables.tree_sequence(),
    )
