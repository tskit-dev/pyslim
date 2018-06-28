import attr
import struct
import msprime

from .slim_tables import *

def load(path, reference_time=0.0):
    '''
    Loads a standard msprime tree sequence, as produced by a SLiM simulation,
    and does the following things to it:

    - shifts time to start from `reference_time`
    - removes `sites.ancestral_state` and replaces with integers,
        retaining these as the first entries in `self.alleles`.
    - removes `mutations.derived_state` and replaces with integers,
        retaining these as the remaining entries in `self.alleles`.

    After this substitution, an allele `j` at site `k` has state
    `self.alleles[k][j]`.

    :param string path: The path to a .trees file.
    :param float reference_time: The time in tree_sequence that will correspond to zero
        in the new tree sequence (this value will be subtracted from all node times).
    '''
    ts = msprime.load(path)
    tables = ts.tables
    return load_tables(tables, reference_time=reference_time)


def load_tables(tables, reference_time=0.0):
    '''
    Loads the TreeSequence defined by the tables.
    '''
    # pull out ancestral states
    alleles = [[x] for x in msprime.unpack_bytes(tables.sites.ancestral_state,
                                                          tables.sites.ancestral_state_offset)]
    derived_state = msprime.unpack_bytes(tables.mutations.derived_state,
                                         tables.mutations.derived_state_offset)
    new_ancestral_state = [b'0' for _ in range(tables.sites.num_rows)]
    new_derived_state = [b'' for _ in derived_state]
    for j in range(tables.mutations.num_rows):
        site = tables.mutations.site[j]
        try:
            allele_index = alleles[site].index(derived_state[j])
        except ValueError:
            allele_index = len(alleles[site])
            alleles[site].append(derived_state[j])

        new_derived_state[j] = bytes(str(allele_index), encoding='utf-8')

    # reset sites and mutations
    new_ds_column, new_ds_offset = msprime.pack_bytes(new_derived_state)
    tables.mutations.set_columns(site=tables.mutations.site, node=tables.mutations.node, 
            derived_state=new_ds_column, derived_state_offset=new_ds_offset, 
            parent=tables.mutations.parent, metadata=tables.mutations.metadata, 
            metadata_offset=tables.mutations.metadata_offset)
    new_as_column, new_as_offset = msprime.pack_bytes(new_ancestral_state)
    tables.sites.set_columns(position=tables.sites.position, 
                             ancestral_state=new_as_column,
                             ancestral_state_offset=new_as_offset,
                             metadata=tables.sites.metadata,
                             metadata_offset=tables.sites.metadata_offset)

    # reset time
    tables.nodes.set_columns(flags=tables.nodes.flags, 
            time=tables.nodes.time - reference_time, 
            population=tables.nodes.population, individual=tables.nodes.individual, 
            metadata=tables.nodes.metadata, metadata_offset=tables.nodes.metadata_offset)
    tables.migrations.set_columns(left=tables.migrations.left, right=tables.migrations.right,
            node=tables.migrations.node, source=tables.migrations.source, 
            dest=tables.migrations.dest, time=tables.migrations.time - reference_time)

    ts = tables.tree_sequence()
    ts.reference_time = reference_time
    ts.alleles = alleles

    return ts


def annotate(tables):
    '''
    Takes a set of tables defining a tree sequence (as produced by msprime, for
    instance), and adds in the information necessary for SLiM to use it as an
    initial state, filling in mostly default values.
    '''
    set_nodes_individuals(tables)
    set_populations(tables)
    set_mutations(tables)
    set_provenances(tables)
    return tables
