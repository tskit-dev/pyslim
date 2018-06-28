import msprime
from .slim_metadata import *


def set_nodes_individuals(
        tables, node_ind=None, location=(0, 0, 0), age=0, ind_id=None,
        ind_population=0, ind_sex=-1, ind_flags=0, slim_ind_flags=0, node_id=None,
        node_is_null=False, node_type=0):
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

    if type(ind_population) is int:
        ind_population = [ind_population for _ in range(num_individuals)]
    assert(len(ind_population) == num_individuals)

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

    node_pop = [msprime.NULL_POPULATION
                if (tables.nodes.flags[u] & msprime.NODE_IS_SAMPLE)
                else ind_population[node_ind[u]]
                for u in range(tables.nodes.num_rows)]

    tables.nodes.set_columns(flags=tables.nodes.flags, time=tables.nodes.time,
                             population=node_pop, individual=node_ind,
                             metadata=tables.nodes.metadata,
                             metadata_offset=tables.nodes.metadata_offset)

    tables.individuals.set_columns(
            flags=ind_flags, 
            location=[x for y in location for x in y],
            location_offset=[0] + [len(x) for x in location])

    individual_metadata = [IndividualMetadata(*x) for x in
                           zip(age, ind_id, ind_population, ind_sex, slim_ind_flags)]
    node_metadata = [NodeMetadata(i, n, t) for i, n, t in
                     zip(node_id, node_is_null, node_type)]

    annotate_individual_metadata(tables, individual_metadata)
    annotate_node_metadata(tables, node_metadata)


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

    population_metadata = [PopulationMetadata(*x) for x in
                           zip(pop_id, selfing_fraction, female_cloning_fraction,
                               male_cloning_fraction, sex_ratio, bounds_x0,
                               bounds_x1, bounds_y0, bounds_y1, bounds_z0, bounds_z1,
                               migration_records)]
    annotate_population_metadata(tables, population_metadata)


def set_sites_mutations(
        tables, mutation_id=None, mutation_type=1, selection_coeff=0.0, sex=-1, flags=0,
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


# The general structure of a Provenance entry is a JSON string:
# {“program”=“SLiM”, “version”=“<version>“, “file_version”=“<file_version>“,
#     “model_type”=“<model_type>“, “generation”=<generation>,
#     “remembered_node_count”=<rem_count>}
# The field values in angle brackets have the following meanings:
# <version>: The version of SLiM that generated the file, such as 3.0.
# <file_version>: The metadata format of the file; at present only 0.1 is supported.
# <model_type>: This should be either WF or nonWF, depending upon the type of
# model that generated the file.  This has some implications for the other
# metadata; in particular, some of the population metadata is required for WF
# models but unused in nonWF models, and individual ages in WF model data are
# expected to be -1.  Note that SLiM will allow WF model data to be read into a
# nonWF model, and vice versa, but since this is usually not intentional a
# warning will be issued.  Moving data between model types has not been tested,
# so be aware that issues may exist with doing so.
# <generation>: The generation number, which will be set by SLiM upon loading.
# <rem_count>: The number of remembered nodes, in addition to the sample.

def set_provenances(tables, model_type, generation=0, remembered_node_count=0):
    '''
    Appends to the Provenance table of a TableCollection a record containing
    the information that SLiM expects to find there.
    '''
    pyslim_record = '"program"="pyslim", "version"="{}"'
    slim_record = '"program"="SLiM", "version"="3.0", "file_version"="0.1", ' \
            + '"model_type"={}, "generation"={}, "remembered_node_count"={}'
    tables.provenances.add_row('{' + pyslim_record.format(0.1) + '}')
    tables.provenances.add_row('{' + slim_record.format(model_type, generation, 
                                                        remembered_node_count) + '}')

