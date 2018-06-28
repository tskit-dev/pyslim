import attr
import struct
import msprime
import json


###########
# Mutations
#  The metadata for a Mutation is a *collection* of these structs,
#  thanks to mutation stacking.
###########
# typedef struct __attribute__((__packed__)) {
#     slim_objectid_t mutation_type_id_;    // 4 bytes (int32_t): the id of the mutation type the mutation belongs to
#     slim_selcoeff_t selection_coeff_;     // 4 bytes (float): the selection coefficient
#     slim_objectid_t subpop_index_;        // 4 bytes (int32_t): the id of the subpopulation in which the mutation arose
#     slim_generation_t origin_generation_; // 4 bytes (int32_t): the generation in which the mutation arose
# } MutationMetadataRec;
#

@attr.s
class MutationMetadata(object):
    mutation_type = attr.ib()
    selection_coeff = attr.ib()
    population = attr.ib()
    time = attr.ib()

def decode_mutation(buff):
    num_muts = int(len(buff) / 16) # 4 + 4 + 4 + 4
    if len(buff) != num_muts * 16:
        raise ValueError("Mutation metadata of incorrect format.")
    struct_string = "<" + "ifii" * num_muts
    metadata = struct.unpack(struct_string, buff)
    mut_structs = []
    for k in range(num_muts):
        mutation_type = metadata[k * 4]
        selection_coeff = metadata[k * 4 + 1]
        population = metadata[k * 4 + 2]
        time = metadata[k * 4 + 3]
        mut_structs.append(MutationMetadata(mutation_type=mutation_type,
                                            selection_coeff=selection_coeff,
                                            population=population, time=time))
    return mut_structs

def encode_mutation(metadata_object):
    '''
    Encodes the list of mutation metadata objects as a bytes object,
    suitable to be put in as metadata for a mutation.
    '''
    mr_values = []
    for mr in metadata_object:
        mr_values.extend([mr.mutation_type, mr.selection_coeff, mr.population, mr.time])
    struct_string = "<" + "ifii" * len(metadata_object)
    return struct.pack(struct_string, *mr_values)


def extract_mutation_metadata(tables):
    '''
    Returns an iterator over lists of MutationMetadata objects containing
    information about the mutations in the tables.
    '''
    metadata = msprime.unpack_bytes(tables.mutations.metadata,
                                    tables.mutations.metadata_offset)
    for md in metadata:
        yield decode_mutation(md)


def annotate_mutation_metadata(tables, metadata):
    '''
    Revise the mutation table so that the metadata column is given by applying
    `encode_mutation()` to the sources given.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of (lists of MutationMetadata) or None objects
    '''
    if len(metadata) != tables.mutations.num_rows:
        raise ValueError("annotate mutations: metadata not the same length as the table.")
    orig_mutations = tables.mutations.copy()
    tables.mutations.clear()
    for mut, md in zip(orig_mutations, metadata):
        if md is None:
            metadata = b''
        else:
            metadata = encode_mutation(md)
        tables.mutations.add_row(site=mut.site, node=mut.node, derived_state=mut.derived_state,
                                 parent=mut.parent, metadata=metadata)

#######
# Nodes
#######
# typedef struct __attribute__((__packed__)) {
#     slim_genomeid_t genome_id_; // 8 bytes (int64_t): the SLiM genome ID for this genome, assigned by pedigree rec
#     uint8_t is_null_;           // 1 byte (uint8_t): true if this is a null genome (should never contain mutations)
#     GenomeType type_;           // 1 byte (uint8_t): the type of the genome (A, X, Y)
# } GenomeMetadataRec;
#

_node_struct = struct.Struct("<qBB")

@attr.s
class NodeMetadata(object):
    slim_id = attr.ib()
    is_null = attr.ib()
    genome_type = attr.ib()

def decode_node(buff):
    if len(buff) != 10: # 8 + 1 + 1
        raise ValueError("Node metadata of incorrect format.")
    slim_id, is_null, genome_type = _node_struct.unpack(buff)
    return NodeMetadata(slim_id=slim_id, is_null=is_null, genome_type=genome_type)

def encode_node(metadata_object):
    return _node_struct.pack(metadata_object.slim_id, metadata_object.is_null,
                            metadata_object.genome_type)


def extract_node_metadata(tables):
    '''
    Returns an iterator over lists of NodeMetadata objects containing
    information about the nodes in the tables.
    '''
    metadata = msprime.unpack_bytes(tables.nodes.metadata,
                                    tables.nodes.metadata_offset)
    for md in metadata:
        yield decode_node(md)


def annotate_node_metadata(tables, metadata):
    '''
    Modify the NodeTable so that the metadata
    column is given by applying `encode_node()` to the sources given.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of NodeMetadata or None objects
    '''
    if len(metadata) != tables.nodes.num_rows:
        raise ValueError("annotate nodes: metadata not the same length as the table.")
    orig_nodes = tables.nodes.copy()
    tables.nodes.clear()
    for node, md in zip(orig_nodes, metadata):
        if md is None:
            metadata = b''
        else:
            metadata = encode_node(md)
        tables.nodes.add_row(flags=node.flags, time=node.time, population=node.population,
                             individual=node.individual, metadata=metadata)


#######
# Individuals
#######
# typedef struct __attribute__((__packed__)) {
#   slim_pedigreeid_t pedigree_id_;    // 8 bytes (int64_t): the SLiM pedigree ID for this individual, assigned by pedigree rec
#   slim_age_t age_;                   // 4 bytes (int32_t): the age of the individual (-1 for WF models)
#   slim_objectid_t subpopulation_id_; // 4 bytes (int32_t): the subpopulation the individual belongs to
#   IndividualSex sex_;    // 4 bytes (int32_t): the sex of the individual, as defined by the IndividualSex enum
#   uint32_t flags_;       // 4 bytes (uint32_t): assorted flags, see below
# } IndividualMetadataRec;

_individual_struct = struct.Struct("<qiiiI")

@attr.s
class IndividualMetadata(object):
    age = attr.ib()
    pedigree_id = attr.ib()
    population = attr.ib()
    sex = attr.ib()
    flags = attr.ib()

def decode_individual(buff):
    if len(buff) != 24: # 8 + 4 + 4 + 4 + 4:
        raise ValueError("Individual metadata of incorrect format.")
    age, pedigree_id, population, sex, flags = _individual_struct.unpack(buff)
    return IndividualMetadata(
                age=age, pedigree_id=pedigree_id, population=population,
                sex=sex, flags=flags)

def encode_individual(metadata_object):
    return _individual_struct.pack(metadata_object.age, metadata_object.pedigree_id,
                                   metadata_object.population, metadata_object.sex,
                                   metadata_object.flags)


def extract_individual_metadata(tables):
    '''
    Returns an iterator over lists of IndividualMetadata objects containing
    information about the individuals in the tables.
    '''
    metadata = msprime.unpack_bytes(tables.individuals.metadata,
                                    tables.individuals.metadata_offset)
    for md in metadata:
        yield decode_individual(md)


def annotate_individual_metadata(tables, metadata):
    '''
    Modify the IndividualTable so that the metadata
    column is given by applying `encode_individual()` to the sources given.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of IndividualMetadata or None objects
    '''
    if len(metadata) != tables.individuals.num_rows:
        raise ValueError("annotate individuals: metadata not the same length as the table.")
    orig_individuals = tables.individuals.copy()
    tables.individuals.clear()
    for ind, md in zip(orig_individuals, metadata):
        if md is None:
            metadata = b''
        else:
            metadata = encode_individual(md)
        tables.individuals.add_row(flags=ind.flags, location=ind.location, metadata=metadata)

#######
# Populations
####################
# typedef struct __attribute__((__packed__)) {
#     slim_objectid_t subpopulation_id_; // 4 bytes (int32_t): the id of this subpopulation
#     double selfing_fraction_;          // 8 bytes (double): selfing fraction (unused in non-sexual models), unused in nonWF models
#     double female_clone_fraction_;     // 8 bytes (double): cloning fraction (females / hermaphrodites), unused in nonWF models
#     double male_clone_fraction_;       // 8 bytes (double): cloning fraction (males / hermaphrodites), unused in nonWF models
#     double sex_ratio_;                 // 8 bytes (double): sex ratio (M:M+F), unused in nonWF models
#     double bounds_x0_;                 // 8 bytes (double): spatial bounds, unused in non-spatial models
#     double bounds_x1_;                 // 8 bytes (double): spatial bounds, unused in non-spatial models
#     double bounds_y0_;                 // 8 bytes (double): spatial bounds, unused in non-spatial / 1D models
#     double bounds_y1_;                 // 8 bytes (double): spatial bounds, unused in non-spatial / 1D models
#     double bounds_z0_;                 // 8 bytes (double): spatial bounds, unused in non-spatial / 1D / 2D models
#     double bounds_z1_;                 // 8 bytes (double): spatial bounds, unused in non-spatial / 1D / 2D models
#     int32_t migration_rec_count_;      // 4 bytes (int32_t): the number of migration records, 0 in nonWF models
#     // followed by migration_rec_count_ instances of SubpopulationMigrationMetadataRec
# } SubpopulationMetadataRec;
#
# typedef struct __attribute__((__packed__)) {
#     slim_objectid_t source_subpop_id_; // 4 bytes (int32_t): the id of the source subpopulation, unused in nonWF models
#     double migration_rate_;            // 8 bytes (double): the migration rate from source_subpop_id_, unused in nonWF models
# } SubpopulationMigrationMetadataRec;

@attr.s
class PopulationMigrationMetadata(object):
    source_subpop = attr.ib()
    migration_rate = attr.ib()

@attr.s
class PopulationMetadata(object):
    slim_id = attr.ib()
    selfing_fraction = attr.ib()
    female_cloning_fraction = attr.ib()
    male_cloning_fraction = attr.ib()
    sex_ratio = attr.ib()
    bounds_x0 = attr.ib()
    bounds_x1 = attr.ib()
    bounds_y0 = attr.ib()
    bounds_y1 = attr.ib()
    bounds_z0 = attr.ib()
    bounds_z1 = attr.ib()
    migration_records = attr.ib()

def decode_population(buff):
    if len(buff) < 88: # 4 + 8 * 10 + 4
        raise ValueError("Population metadata of incorrect format.")
    num_migration_records = int((len(buff) - 88) / 12) # 4 + 8
    if len(buff) != 88 + 12 * num_migration_records:
        raise ValueError("Population metadata of incorrect format.")
    struct_string = "<i10di" + "id" * num_migration_records
    metadata = struct.unpack(struct_string, buff)
    if metadata[11] != num_migration_records:
        raise ValueError("Inconsistent population metadata format.")
    slim_id = metadata[0]
    selfing_fraction = metadata[1]
    female_cloning_fraction = metadata[2]
    male_cloning_fraction = metadata[3]
    sex_ratio = metadata[4]
    bounds_x0 = metadata[5]
    bounds_x1 = metadata[6]
    bounds_y0 = metadata[7]
    bounds_y1 = metadata[8]
    bounds_z0 = metadata[9]
    bounds_z1 = metadata[10]
    # num_migration_records = metadata[11]
    migration_records = []
    k = 12
    for j in range(num_migration_records):
        source_subpop = metadata[k]
        k += 1
        migration_rate = metadata[k]
        k += 1
        migration_records.append(PopulationMigrationMetadata(source_subpop=source_subpop,
                                                             migration_rate=migration_rate))

    return PopulationMetadata(slim_id=slim_id, selfing_fraction=selfing_fraction,
                              female_cloning_fraction=female_cloning_fraction,
                              male_cloning_fraction=male_cloning_fraction,
                              sex_ratio=sex_ratio, bounds_x0=bounds_x0, bounds_x1=bounds_x1,
                              bounds_y0=bounds_y0, bounds_y1=bounds_y1, bounds_z0=bounds_z0,
                              bounds_z1=bounds_z1, migration_records=migration_records)

def encode_population(metadata_object):
    num_migration_records = len(metadata_object.migration_records)
    for mr in metadata_object.migration_records:
        if len(mr) != 2:
            raise ValueError("Migration records should be tuples of (source, rate).")
    mr_values = [a for b in metadata_object.migration_records for a in b]
    struct_string = "<i10di" + "id" * num_migration_records
    metadata = struct.pack(struct_string, metadata_object.slim_id,
                           metadata_object.selfing_fraction, metadata_object.female_cloning_fraction,
                           metadata_object.male_cloning_fraction, metadata_object.sex_ratio,
                           metadata_object.bounds_x0, metadata_object.bounds_x1,
                           metadata_object.bounds_y0, metadata_object.bounds_y1,
                           metadata_object.bounds_z0, metadata_object.bounds_z1,
                           num_migration_records, *mr_values)
    return metadata


def extract_population_metadata(tables):
    '''
    Returns an iterator over lists of PopulationMetadata objects containing
    information about the populations in the tables.
    '''
    metadata = msprime.unpack_bytes(tables.populations.metadata,
                                    tables.populations.metadata_offset)
    for md in metadata:
        yield decode_population(md)


def annotate_population_metadata(tables, metadata):
    '''
    Modify the PopulationTable so that the metadata
    column is given by applying `encode_population()` to the sources given.
    This entirely removes existing information in the Population table.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of objects, each PopulationMetadata or None
    '''
    if len(metadata) != tables.populations.num_rows:
        raise ValueError("annotate populations: metadata not the same length as the table.")
    tables.populations.clear()
    for md in metadata:
        if md is None:
            metadata = b''
        else:
            metadata = encode_population(md)
        tables.populations.add_row(metadata=metadata)

#######
# Provenance
####################
# The general structure of a Provenance entry is a JSON string:
# {“program”:“SLiM”, “version”:“<version>“, “file_version”:“<file_version>“,
#     “model_type”:“<model_type>“, “generation”:<generation>,
#     “remembered_node_count”:<rem_count>}
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


@attr.s
class ProvenanceMetadata(object):
    model_type = attr.ib()
    slim_generation = attr.ib()
    remembered_node_count = attr.ib()


def get_provenance(tables):
    '''
    Extracts model type, slim generation, and remembmered node count from the last
    entry in the provenance table that is tagged with "program"="SLiM".
    '''
    prov = [json.loads(x.record) for x in tables.provenances]
    slim_prov = [u for u in prov if ('program' in u and u['program'] == "SLiM")]
    if len(slim_prov) == 0:
        raise ValueError("Tree sequence contains no SLiM provenance entries.")
    last_slim_prov = slim_prov[len(slim_prov)-1]
    return ProvenanceMetadata(last_slim_prov["model_type"], last_slim_prov["generation"],
                              last_slim_prov["remembered_node_count"])


def set_provenance(tables, model_type, slim_generation, remembered_node_count=0):
    '''
    Appends to the Provenance table of a TableCollection a record containing
    the information that SLiM expects to find there.
    '''
    pyslim_dict = {"program":"pyslim", "version":0.1}
    slim_dict = {"program":"SLiM", "version":"3.0", "file_version":"0.1",
                 "model_type":model_type, "generation":slim_generation,
                 "remembered_node_count":remembered_node_count}
    tables.provenances.add_row(json.dumps(pyslim_dict))
    tables.provenances.add_row(json.dumps(slim_dict))

