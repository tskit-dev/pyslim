import msprime
import math
import numpy as np

from .slim_tree_sequence import *

def infinite_alleles_factory(alleles):
    next_allele = str(np.random.choice(1000000000000))
    while True:
        while next_allele in alleles:
            next_allele = str(np.random.choice(1000000000000))
        alleles.add(next_allele)
        yield next_allele


def mutate(ts, mutation_rate):
    tables = ts.dump_tables()
    edge_order = tables.nodes.time[tables.edges.child].argsort()
    site_pos = {}
    alleles = set(msprime.unpack_strings(tables.sites.ancestral_state,
                                         tables.sites.ancestral_state_offset))
    alleles.update(msprime.unpack_strings(tables.mutations.derived_state,
                                          tables.mutations.derived_state_offset))
    ag = infinite_alleles_factory(alleles)
    for j, site in enumerate(ts.sites()):
        site_pos[site.position] = j


    for j in reversed(edge_order):
        e = tables.edges[j]
        left = math.ceil(e.left)
        right = math.floor(e.right)
        num_bp = 1 + right - left - (right == e.right)
        if num_bp > 0:
            edge_len = ts.node(e.parent).time - ts.node(e.child).time
            mut_prob = 1 - math.exp(-1 * edge_len * mutation_rate)
            num_muts = np.random.binomial(num_bp, mut_prob)
            mut_bp = left + np.random.choice(num_bp, size=num_muts, replace=False)
            for pos in mut_bp:
                if pos not in site_pos:
                    site = tables.sites.add_row(pos, ancestral_state='0')
                    next_allele = 1
                    site_pos[pos] = site
                else:
                    site = site_pos[pos]
                
                tables.mutations.add_row(site, e.child, derived_state=next(ag))

    tables.sort()
    tables.compute_mutation_parents()
    delabel_alleles(tables)

    return tables.tree_sequence()
