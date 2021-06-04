import msprime
import tskit
import warnings

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
    ``keep_input_roots=True`` to :meth:`.simplify()`.

    If you specify an ``ancestral_Ne``, then the recapitated portion of the
    tree sequence will be simulated in a single population with this
    (diploid) size. In other words, all lineages are moved to a single
    population of this size (named "ancestral" if this name is not already
    taken), and coalescence is allowed to happen.

    You may control the ancestral demography by passing in a ``demography``
    argument: see {ref}``msprime.sim_ancestry``.
    
    In general, all defaults are whatever the defaults of
    {ref}`msprime.sim_ancestry` are; this includes recombination rate, so
    that if neither ``recombination_rate`` or a ``recombination_map`` are
    provided, there will be *no* recombination.

    :param float ancestral_Ne: If specified, then will simulate from a single
        ancestral population of this size. It is an error to specify this
        as well as ``demography``.
    :param dict kwargs: Any other arguments to :meth:`msprime.sim_ancestry`.
    '''
    demography = None
    if "demography" in kwargs:
        demography = kwargs['demography']
    if ancestral_Ne is not None:
        if demography is not None:
            raise ValueError("You cannot specify both `demography` and `ancestral_Ne`.")
        demography = msprime.Demography.from_tree_sequence(ts)
        # must set pop sizes to >0 even though we merge immediately
        for pop in demography.populations:
            pop.initial_size=1.0
        ancestral_name = "ancestral"
        derived_names = [pop.name for pop in demography.populations]
        while ancestral_name in derived_names:
            ancestral_name = (ancestral_name + "_ancestral")
        ancestral_metadata = default_slim_metadata('population')
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
                    ts.slim_generation,
                    2 * ts.slim_generation,
                ),
                derived=derived_names,
                ancestral=ancestral_name,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", msprime.IncompletePopulationMetadataWarning)
        recap = msprime.sim_ancestry(
                    initial_state = ts,
                    demography = demography,
                    **kwargs)
    return SlimTreeSequence(recap, reference_sequence=ts.reference_sequence)


