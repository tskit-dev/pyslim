import attr
import struct

## Here is how info is packed into objects:
# from slim_sim.h
#
# typedef struct __attribute__((__packed__)) {
#     slim_objectid_t mutation_type_id_;    // 4 bytes (int32_t): the id of the mutation type the mutation belongs to
#     slim_selcoeff_t selection_coeff_;     // 4 bytes (float): the selection coefficient
#     slim_objectid_t subpop_index_;        // 4 bytes (int32_t): the id of the subpopulation in which the mutation arose
#     slim_generation_t origin_generation_; // 4 bytes (int32_t): the generation in which the mutation arose
# } MutationMetadataRec;
#
# typedef struct __attribute__((__packed__)) {
#     slim_genomeid_t genome_id_;   // 8 bytes (int64_t): the SLiM genome ID for this genome, assigned by pedigree rec
#     uint8_t is_null_;             // 1 byte (uint8_t): true if this is a null genome (should never contain mutations)
#     GenomeType type_;             // 1 byte (uint8_t): the type of the genome (A, X, Y)
# } GenomeMetadataRec;
#
# typedef struct __attribute__((__packed__)) {
#     slim_pedigreeid_t pedigree_id_;     // 8 bytes (int64_t): the SLiM pedigree ID for this individual, assigned by pedigree rec
#     slim_age_t age_;                    // 4 bytes (int32_t): the age of the individual (-1 for WF models)
#     slim_objectid_t subpopulation_id_;  // 4 bytes (int32_t): the subpopulation the individual belongs to
# } IndividualMetadataRec;
#
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
class MutationMetadata(object):
    mutation_type = attr.ib()
    selection_coeff = attr.ib()
    population = attr.ib()
    time = attr.ib()

def decode_mutation(buff):
    if len(buff) != 16: # 4 + 4 + 4 + 4
        raise ValueError("Mutation metadata of incorrect format.")
    mutation_type, selection_coeff, population, time = struct.unpack("ifii", buff)
    return MutationMetadata(mutation_type=mutation_type, selection_coeff=selection_coeff,
                            population=population, time=time)


@attr.s
class NodeMetadata(object):
    slim_id = attr.ib()
    is_null = attr.ib()
    genome_type = attr.ib()

def decode_node(buff):
    if len(buff) != 10: # 8 + 1 + 1
        raise ValueError("Node metadata of incorrect format.")
    slim_id, is_null, genome_type = struct.unpack("qBB", buff)
    return NodeMetadata(slim_id=slim_id, is_null=is_null, genome_type=genome_type)


@attr.s
class IndividualMetadata(object):
    age = attr.ib()
    pedigree_id = attr.ib()
    population = attr.ib()

def decode_individual(buff):
    if len(buff) != 16: # 8 + 4 + 4:
        raise ValueError("Individual metadata of incorrect format.")
    age, pedigree_id, population = struct.unpack("qii", buff)
    return IndividualMetadata(age=age, pedigree_id=pedigree_id, population=population)


@attr.s
class MigrationMetadata(object):
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
    struct_string = "=i10di" + "id" * num_migration_records
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
        migration_records.append(MigrationMetadata(source_subpop = source_subpop,
                                                   migration_rate = migration_rate))

    return PopulationMetadata(slim_id = slim_id, selfing_fraction = selfing_fraction,
        female_cloning_fraction = female_cloning_fraction,
        male_cloning_fraction = male_cloning_fraction,
        sex_ratio = sex_ratio, bounds_x0 = bounds_x0, bounds_x1 = bounds_x1,
        bounds_y0 = bounds_y0, bounds_y1 = bounds_y1, bounds_z0 = bounds_z0,
        bounds_z1 = bounds_z1, migration_records = migration_records)

