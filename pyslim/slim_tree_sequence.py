import attr
import struct
import msprime
import json

from .slim_metadata import *

GENOME_TYPE_AUTOSOME = 0
GENOME_TYPE_X = 1
GENOME_TYPE_Y = 2
INDIVIDUAL_TYPE_HERMAPHRODITE = -1
INDIVIDUAL_TYPE_FEMALE = 0
INDIVIDUAL_TYPE_MALE = 1
INDIVIDUAL_FLAG_MIGRATED = 0x01

def load(path, slim_format):
    '''
    Loads a tree sequence, and if slim_format is True, then
    do the following things to it:

    - shifts time to start from the `generation` entry in the provenance table
    - removes `sites.ancestral_state` and `mutations.derived_state`, replacing
      these with integers, and storing translation tables in sites.metadata.

    After this substitution, the original allele `j` at site `k` can be obtained
    by doing
    ```
    alleles = extract_alleles(tables)
    alleles[k][j]
    ```

    These two operations are not necessary for working with the tree sequence,
    but are more convenient.

    :param string path: The path to a .trees file.
    :param bool slim_format: Whether the .trees file should be coverted from
        SLiM format.
    '''
    ts = msprime.load(path)
    if slim_format:
        tables = ts.tables
        ts = load_tables(tables, slim_format)
    return ts


def dump(ts, path, slim_format):
    '''
    Write out the .trees file that can be read back in by SLiM.
    '''
    if slim_format:
        tables = ts.tables
        provenance = get_provenance(tables)
        _set_slim_generation(tables, -1 * provenance.slim_generation)
        _relabel_alleles(tables)
    ts.dump(path)


def load_tables(tables, slim_format):
    '''
    Loads the TreeSequence defined by the tables: see `pyslim.load()`.
    '''
    if slim_format:
        provenance = get_provenance(tables)
        _set_slim_generation(tables, provenance.slim_generation)
        _delabel_alleles(tables)
    ts = tables.tree_sequence()
    return ts


def annotate(ts, model_type, slim_generation, remembered_node_count=0):
    '''
    Takes tree sequence (as produced by msprime, for instance), and adds in the
    information necessary for SLiM to use it as an initial state, filling in
    mostly default values.
    '''
    tables = ts.tables
    annotate_tables(tables, model_type, slim_generation, remembered_node_count)
    ts = tables.tree_sequence()
    return ts


def annotate_tables(tables, model_type, slim_generation, remembered_node_count=0):
    '''
    Does the work of `annotate()`, but modifies the tables in place.
    '''
    if (type(slim_generation) is not int) or (slim_generation < 1):
        raise ValueError("SLiM generation must be an integer and at least 1.")
    # nodes must come before populations
    _set_nodes_individuals(tables)
    _set_populations(tables)
    _set_sites_mutations(tables)
    set_provenance(tables, model_type=model_type, slim_generation=slim_generation,
                   remembered_node_count=remembered_node_count)
    ancestral_state = msprime.unpack_strings(tables.sites.ancestral_state,
                                             tables.sites.ancestral_state_offset)
    derived_state = msprime.unpack_strings(tables.mutations.derived_state,
                                           tables.mutations.derived_state_offset)
    alleles = _make_allele_map(tables)
    _store_alleles(tables, alleles)


def extract_alleles(tables, keep=True):
    '''
    Decodes the list of dictionaries giving translations for each allele
    stored in the metadata column of the site table when reading a SLiM tree sequence.
    If `keep` is False, will clear the metadata column.
    '''
    md = msprime.unpack_strings(tables.sites.metadata, tables.sites.metadata_offset)
    alleles = [json.loads(x) for x in md]
    if not keep:
        tables.sites.set_columns(
                position=tables.sites.position,
                ancestral_state=tables.sites.ancestral_state,
                ancestral_state_offset=tables.sites.ancestral_state_offset,
                metadata=None, metadata_offset=None)
    return alleles


def _store_alleles(tables, alleles):
    '''
    Stores an allele map in the metadata column of a site table, overwriting
    anything already there.
    '''
    alleles_text = [json.dumps(x) for x in alleles]
    md_val, md_off = msprime.pack_strings(alleles_text)
    tables.sites.set_columns(
            position=tables.sites.position,
            ancestral_state=tables.sites.ancestral_state,
            ancestral_state_offset=tables.sites.ancestral_state_offset,
            metadata=md_val, metadata_offset=md_off)


def _make_allele_map(tables):
    '''
    Make a list of dicts whose keys are the alleles at each site and whose
    values are a sequence of integers, starting from 0, as strings.  Used to
    make an allele map for an msprime-produced tree sequence that can later be
    used to produce a SLiM-like tree sequence.
    '''
    inv_alleles = [{x:'0'} for x in msprime.unpack_strings(tables.sites.ancestral_state,
                                                           tables.sites.ancestral_state_offset)]
    derived_state = msprime.unpack_strings(tables.mutations.derived_state,
                                           tables.mutations.derived_state_offset)
    for j in range(tables.mutations.num_rows):
        site = tables.mutations.site[j]
        if derived_state[j] not in inv_alleles[site]:
            inv_alleles[site][derived_state[j]] = str(len(inv_alleles[site]))
    return inv_alleles


def _delabel_alleles(tables):
    '''
    Replaces alleles at each site by integers, starting with '0' for the ancestral state,
    placing the resulting list of dictionaries in the metadata column of the site table
    (and removing anything already there).  Can be inverted by 
    ```
       _relabel_alleles(tables)
    ```
    Used to make SLiM-produced alleles nice to look at.
    '''
    inv_alleles = _make_allele_map(tables)
    alleles = [dict((v, k) for k, v in x.items()) for x in inv_alleles]
    _relabel_alleles(tables, inv_alleles)
    _store_alleles(tables, alleles)
    return alleles


def _relabel_alleles(tables, alleles=None):
    '''
    Replace ancestral and derived state x at site j by alleles[j][x].
    `alleles` should be a list of dicts. If alleles are not provided,
    extract them from the tables (removing them from the metadata).
    '''
    if alleles is None:
        alleles = extract_alleles(tables, keep=False)
    ancestral_state = msprime.unpack_strings(tables.sites.ancestral_state,
                                             tables.sites.ancestral_state_offset)
    derived_state = msprime.unpack_strings(tables.mutations.derived_state,
                                           tables.mutations.derived_state_offset)
    for j, x in enumerate(ancestral_state):
        if x not in alleles[j]:
            raise ValueError("Ancestral state {} not present in alleles at site {}.".format(x, j))
    for j, x in enumerate(derived_state):
        if x not in alleles[tables.mutations.site[j]]:
            raise ValueError("Derived state {} not present in alleles " \
                    + "at site {}.".format(x, tables.mutations.site[j]))
    new_ancestral_state = [alleles[j][x] for j, x in enumerate(ancestral_state)]
    new_derived_state = [alleles[j][x] for j, x in zip(tables.mutations.site, derived_state)]
    # reset sites and mutations
    new_ds_column, new_ds_offset = msprime.pack_strings(new_derived_state)
    tables.mutations.set_columns(site=tables.mutations.site, node=tables.mutations.node,
            derived_state=new_ds_column, derived_state_offset=new_ds_offset,
            parent=tables.mutations.parent, metadata=tables.mutations.metadata,
            metadata_offset=tables.mutations.metadata_offset)
    new_as_column, new_as_offset = msprime.pack_strings(new_ancestral_state)
    tables.sites.set_columns(position=tables.sites.position,
                             ancestral_state=new_as_column,
                             ancestral_state_offset=new_as_offset,
                             metadata=tables.sites.metadata,
                             metadata_offset=tables.sites.metadata_offset)


def _set_slim_generation(tables, slim_generation):
    '''
    Shifts the "time ago" entries in the tables to be measured in units of time
    *before* `slim_generation`, by adding this value to the `time` columns of
    Node and Migration tables.
    '''
    tables.nodes.set_columns(flags=tables.nodes.flags,
            time=tables.nodes.time + slim_generation,
            population=tables.nodes.population, individual=tables.nodes.individual,
            metadata=tables.nodes.metadata, metadata_offset=tables.nodes.metadata_offset)
    tables.migrations.set_columns(left=tables.migrations.left, right=tables.migrations.right,
            node=tables.migrations.node, source=tables.migrations.source,
            dest=tables.migrations.dest, time=tables.migrations.time + slim_generation)


def _set_nodes_individuals(
        tables, node_ind=None, location=(0, 0, 0), age=0, ind_id=None,
        ind_population=None, ind_sex=INDIVIDUAL_TYPE_HERMAPHRODITE,
        ind_flags=0, slim_ind_flags=0, node_id=None,
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
    - (age) individual age to 0.0
    - (ind_id) SLiM individual pedigree IDs to sequential integers starting from 0
    - (ind_population) individual populations to 0
    - (node_id) SLiM genome IDs to sequential integers starting with samples from 0
    - (node_is_null) genomes to be non-null
    - (node_type) genome type to 0 (= autosome)
    '''
    samples = list(filter(lambda j: tables.nodes.flags[j] & msprime.NODE_IS_SAMPLE,
                          range(tables.nodes.num_rows)))

    if node_ind is None:
        node_ind = [msprime.NULL_INDIVIDUAL for _ in range(tables.nodes.num_rows)]
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
        ind_population = [msprime.NULL_POPULATION for _ in range(num_individuals)]
        for j, u in enumerate(node_ind):
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

    tables.nodes.set_columns(flags=tables.nodes.flags, time=tables.nodes.time,
                             population=tables.nodes.population, individual=node_ind,
                             metadata=tables.nodes.metadata,
                             metadata_offset=tables.nodes.metadata_offset)

    tables.individuals.set_columns(
            flags=ind_flags, location=[x for y in location for x in y],
            location_offset=[0] + [len(x) for x in location])

    individual_metadata = [IndividualMetadata(*x) for x in
                           zip(age, ind_id, ind_population, ind_sex, slim_ind_flags)]
    node_metadata = [None for _ in range(num_nodes)]
    for j in samples:
        node_metadata[j] = NodeMetadata(slim_id=node_id[j], is_null=node_is_null[j],
                                        genome_type=node_type[j])

    annotate_individual_metadata(tables, individual_metadata)
    annotate_node_metadata(tables, node_metadata)


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
    for md in msprime.unpack_bytes(tables.individuals.metadata,
                                   tables.individuals.metadata_offset):
        try:
            ind_md = decode_individual(md)
        except:
            raise ValueError("Individuals do not have metadata:"
                    + "need to run set_nodes_individuals() first?")
        assert(ind_md.population < num_pops)
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
            assert(type(mr) is PopulationMigrationMetadata)

    population_metadata = [PopulationMetadata(*x) for x in
                           zip(pop_id, selfing_fraction, female_cloning_fraction,
                               male_cloning_fraction, sex_ratio, bounds_x0,
                               bounds_x1, bounds_y0, bounds_y1, bounds_z0, bounds_z1,
                               migration_records)]
    annotate_population_metadata(tables, population_metadata)


def _set_sites_mutations(
        tables, mutation_id=None, mutation_type=1, selection_coeff=0.0,
        population=msprime.NULL_POPULATION, time=None):
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

    if time is None:
        ## This may *not* make sense because we have to round:
        # time = [int(tables.nodes.time[u]) for u in tables.mutations.node]
        time = [0 for _ in range(num_mutations)]
    assert(len(time) == num_mutations)

    mutation_metadata = [[MutationMetadata(*x)] for x in
                         zip(mutation_type, selection_coeff, population, time)]
    annotate_mutation_metadata(tables, mutation_metadata)
