import msprime
import tskit
import json
import warnings
import numpy as np

from ._version import *
from .slim_metadata import _old_metadata_schema
from pyslim import NUCLEOTIDES

def load(*args, **kwargs):
    raise RuntimeError("This method has been removed: use tskit.load( ) instead.")


def mutation_at(ts, node, position, time=None):
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
    if position < 0 or position >= ts.sequence_length:
        raise ValueError("Position {} not valid.".format(position))
    if node < 0 or node >= ts.num_nodes:
        raise ValueError("Node {} not valid.".format(node))
    if time is None:
        time = ts.node(node).time
    tree = ts.at(position)
    out = tskit.NULL
    try:
        site = ts.site(position=position)
    except ValueError:
        pass
    else:
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
                if mut.node == n and mut.time >= time:
                    # BUG: this can fail if a mutation has two children
                    assert(out == tskit.NULL or out == mut.parent)
                    out = mut.id
    return out

def nucleotide_at(ts, node, position, time=None):
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
    if not ts.has_reference_sequence():
        raise ValueError("This tree sequence has no reference sequence.")
    mut_id = mutation_at(ts, node, position, time)
    if mut_id == tskit.NULL:
        out = NUCLEOTIDES.index(ts.reference_sequence.data[int(position)])
    else:
        mut = ts.mutation(mut_id)
        k = np.argmax([u["slim_time"] for u in mut.metadata["mutation_list"]])
        out = mut.metadata["mutation_list"][k]["nucleotide"]
    return out
