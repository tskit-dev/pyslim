import msprime
from .slim_metadata import *

# these aren't exported by msprime
MSP_INDIVIDUAL_MALE = 1
MSP_INDIVIDUAL_FEMALE = 2

def set_nodes_individuals(
        tables, node_ind=None, location=(0, 0, 0), age=0, ind_id=None,
        ind_flags=MSP_INDIVIDUAL_MALE | MSP_INDIVIDUAL_FEMALE, ind_population=0, 
        node_id=None, node_is_null=False, node_type=0):
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
    samples = filter(lambda j: tables.nodes.flags[j] & msprime.NODE_IS_SAMPLE,
                     range(tables.nodes.num_rows))

    if node_ind is None:
        node_ind = [msprime.NULL_INDIVIDUAL for _ in range(tables.nodes.num_rows)]
        for j, k in enumerate(samples):
            node_ind[j] = int(k/2)

    num_individuals = max(node_ind) + 1

    if type(location) is tuple:
        location = [location for _ in range(num_individuals)]
    assert(len(location) == num_individuals)

    if type(age) is int or type(age) is float:
        age = [age for _ in range(num_individuals)]
    assert(len(age) == num_individuals)

    if ind_id is None:
        ind_id = list(range(num_individuals))
    assert(len(ind_id) == num_individuals)

    if type(ind_flags) is int:
        ind_flags = [ind_flags for _ in range(num_individuals)]
    assert(len(ind_flags) == num_individuals)

    if type(ind_population) is int:
        ind_population = [ind_population for _ in range(num_individuals)]
    assert(len(ind_population) == num_individuals)

    if node_id is None:
        node_id = [-1 for _ in range(tables.nodes.num_rows)]
        for j, k in enumerate(list(samples) 
                              + sorted(list(set(range(tables.nodes.num_rows)) 
                                            - set(samples)))):
            node_id[k] = j
    assert(len(node_id) == tables.nodes.num_rows)

    if type(node_is_null) is bool:
        node_is_null = [node_is_null for _ in range(tables.nodes.num_rows)]
    assert(len(node_is_null) == tables.nodes.num_rows)

    if type(node_type) is int:
        node_type = [node_type for _ in range(tables.nodes.num_rows)]
    assert(len(node_type) == tables.nodes.num_rows)

    tables.individuals.clear()
    for j in range(num_individuals):
        metadata_obj = IndividualMetadata(age=age[j], pedigree_id=ind_id[j],
                                          population=ind_population[j])
        tables.individuals.add_row(flags=ind_flags[j], location=location[j],
                                   metadata=encode_individual(metadata_obj))

    old_nodes = tables.nodes.copy()
    tables.nodes.clear()
    for j, node in enumerate(old_nodes):
        metadata_obj = NodeMetadata(slim_id=node_id[j], is_null=node_is_null[j],
                                    genome_type=node_type[j])
        if node.flags & msprime.NODE_IS_SAMPLE:
            pop = ind_population[node_ind[j]]
        else:
            pop = msprime.NULL_POPULATION
        tables.nodes.add_row(flags=node.flags, time=node.time, 
                             population=pop, individual=node_ind[j],
                             metadata=encode_node(metadata_obj))


def set_populations(
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

    tables.populations.clear()
    for j in range(num_pops):
        metadata_obj = PopulationMetadata(
                            slim_id=pop_id[j],
                            selfing_fraction=selfing_fraction[j],
                            female_cloning_fraction=female_cloning_fraction[j],
                            male_cloning_fraction=male_cloning_fraction[j],
                            sex_ratio=sex_ratio[j], bounds_x0=bounds_x0[j],
                            bounds_x1=bounds_x1[j], bounds_y0=bounds_y0[j],
                            bounds_y1=bounds_y1[j], bounds_z0=bounds_z0[j],
                            bounds_z1=bounds_z1[j],
                            migration_records=migration_records[j])
        tables.populations.add_row(metadata=encode_population(metadata_obj))


def set_mutations(
        tables, mutation_type=1, selection_coeff=0.0,
        population=msprime.NULL_POPULATION, time=-1):
    '''
    Adds to a TableCollection the information relevant to mutations requiredj
    for SLiM to load in a tree sequence, that is found in the metadata column
    of the Mutation table.  This will replace any information already in the
    metadata column of the Mutation table.
    '''
    num_mutations = tables.mutations.num_rows

    if type(mutation_type) is int:
        mutation_type = [mutation_type for _ in range(num_mutations)]
    assert(len(mutation_type) == num_mutations)

    if type(selection_coeff) is float:
        selection_coeff = [selection_coeff for _ in range(num_mutations)]
    assert(len(selection_coeff) == num_mutations)

    if type(population) is int:
        population = [population for _ in range(num_mutations)]
    assert(len(population) == num_mutations)

    if type(time) is int:
        time = [time for _ in range(num_mutations)]
    assert(len(time) == num_mutations)

    old_mutations = tables.mutations.copy()
    tables.mutations.clear()
    for j, mut in enumerate(old_mutations):
        metadata_obj = [MutationMetadata(mutation_type=mutation_type[j],
                                         selection_coeff=selection_coeff[j],
                                         population=population[j],
                                         time=time[j])]
        tables.mutations.add_row(site=mut.site, node=mut.node, 
                                 derived_state=mut.derived_state, parent=-1, 
                                 metadata=encode_mutation(metadata_obj))


def set_provenances(tables, generation=0, remembered_node_count=0):
    '''
    Appends to the Provenance table of a TableCollection a record containing
    the information that SLiM expects to find there.
    '''
    pyslim_record = '{"program"="pyslim", "version"="{}"}'
    slim_record = '{"program"="SLiM", "version"="pyslim", "file_version"="0.1", "generation"={}, "remembered_node_count"={}}'
    tables.provenances.add_row(pyslim.record.format(0.1))
    tables.provenances.add_row(slim.record.format(generation, remembered_node_count))

