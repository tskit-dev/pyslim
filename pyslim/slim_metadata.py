import attr
import struct
import tskit
import json

GENOME_TYPE_AUTOSOME = 0
GENOME_TYPE_X = 1
GENOME_TYPE_Y = 2
INDIVIDUAL_TYPE_HERMAPHRODITE = -1
INDIVIDUAL_TYPE_FEMALE = 0
INDIVIDUAL_TYPE_MALE = 1
INDIVIDUAL_FLAG_MIGRATED = 0x01

###########
# Mutations
#  The metadata for a Mutation is a *collection* of these structs,
#  thanks to mutation stacking.
###########
# typedef struct __attribute__((__packed__)) {
#         slim_objectid_t mutation_type_id_;    // 4 bytes (int32_t): the id of the mutation type the mutation belongs to
#         slim_selcoeff_t selection_coeff_;     // 4 bytes (float): the selection coefficient
#         slim_objectid_t subpop_index_;        // 4 bytes (int32_t): the id of the subpopulation in which the mutation arose
#         slim_generation_t origin_generation_; // 4 bytes (int32_t): the generation in which the mutation arose
#         int8_t nucleotide_;                   // 1 byte (int8_t): the nucleotide for the mutation (0='A', 1='C', 2='G', 3='T'), or -1 (added in file format v0.2)
# } MutationMetadataRec;
#

@attr.s
class MutationMetadata(object):
    mutation_type = attr.ib()
    selection_coeff = attr.ib()
    population = attr.ib()
    slim_time = attr.ib()
    nucleotide = attr.ib()

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)


def decode_mutation(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`MutationTable`.  If ``buff`` is empty, returns ``[]``.
    If `buff` is already a list of MutationMetadata objects, this will return 
    it, unchanged.

    Note that unlike other metadata decoding functions, this returns a 
    *list* of MutationMetadata objects, because (thanks to mutation stacking),
    a given "mutation" as recorded in tskit may actually represent a
    combination of several SLiM mutations.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`MutationTable`, as stored by SLiM.
    :rtype list:
    '''
    if type(buff) == type([]) and (len(buff) == 0
            or type(buff[0]) == MutationMetadata):
        mut_structs = buff
    else:
        # note that in the case that buff is of length zero
        # this returns [] instead of None, like the others do
        try:
            num_muts = int(len(buff) / 17) # 4 + 4 + 4 + 4 + 1
        except:
            raise ValueError("Metadata not a recognized format.")
        if len(buff) != num_muts * 17:
            raise ValueError("Metadata bytes of incorrect format.")
        struct_string = "<" + "ifiib" * num_muts
        try:
            metadata = struct.unpack(struct_string, buff)
        except:
            raise ValueError("Metadata bytes of incorrect format.")
        mut_structs = []
        for k in range(num_muts):
            mutation_type = metadata[k * 5]
            selection_coeff = metadata[k * 5 + 1]
            population = metadata[k * 5 + 2]
            slim_time = metadata[k * 5 + 3]
            nucleotide = metadata[k * 5 + 4]
            mut_structs.append(MutationMetadata(mutation_type = mutation_type,
                                                selection_coeff = selection_coeff,
                                                population = population,
                                                slim_time = slim_time,
                                                nucleotide = nucleotide))
    return mut_structs


def _decode_mutation_pre_nucleotides(buff):
    '''
    Decodes mutation metadata for file versions 0.1 and 0.2, before
    the 'nucleotide' was added.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`MutationTable`, as stored by SLiM.
    :rtype list:
    '''
    # note that in the case that buff is of length zero
    # this returns [] instead of None, like the others do
    num_muts = int(len(buff) / 16) # 4 + 4 + 4 + 4
    if len(buff) != num_muts * 16:
        raise ValueError("Metadata bytes of incorrect format.")
    struct_string = "<" + "ifii" * num_muts
    metadata = struct.unpack(struct_string, buff)
    mut_structs = []
    for k in range(num_muts):
        mutation_type = metadata[k * 4]
        selection_coeff = metadata[k * 4 + 1]
        population = metadata[k * 4 + 2]
        slim_time = metadata[k * 4 + 3]
        nucleotide = -1
        mut_structs.append(MutationMetadata(mutation_type=mutation_type,
                                            selection_coeff=selection_coeff,
                                            population=population,
                                            slim_time=slim_time,
                                            nucleotide = nucleotide))
    return mut_structs


def encode_mutation(metadata_object):
    '''
    Encodes the list of :class:`MutationMetadata` objects as a bytes object,
    suitable to be put in as metadata for a mutation.  Unlike other ``encode_``
    functions, this takes a *list* rather than a single value, thanks to
    stacking of SLiM mutations.

    :param MutationMetadata metadata_object: The list of
        :class:`MutationMetadata` objects to be encoded.
    :rtype bytes:
    '''
    mr_values = []
    for mr in metadata_object:
        mr_values.extend([mr.mutation_type, mr.selection_coeff,
                          mr.population, mr.slim_time, mr.nucleotide])
    struct_string = "<" + "ifiib" * len(metadata_object)
    return struct.pack(struct_string, *mr_values)


def extract_mutation_metadata(tables):
    '''
    Returns an iterator over lists of :class:`MutationMetadata` objects containing
    information about the mutations in the tables.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    metadata = tskit.unpack_bytes(tables.mutations.metadata,
                                    tables.mutations.metadata_offset)
    for md in metadata:
        yield decode_mutation(md)


def annotate_mutation_metadata(tables, metadata):
    '''
    Revise the mutation table in place so that the metadata column is given by
    applying `encode_mutation()` to the sources given.

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

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

def decode_node(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`NodeTable`.  If the buffer is empty, returns None.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`NodeTable`, as stored by SLiM.
    :rtype NodeMetadata:
    '''
    if type(buff) == NodeMetadata:
        md = buff
    else:
        try:
            buff_len = len(buff)
        except:
                raise ValueError("Metadata not a recognized format.")
        if buff_len == 0:
            md = None
        else:
            if len(buff) != 10: # 8 + 1 + 1
                raise ValueError("Metadata bytes of incorrect length.")
            try:
                slim_id, is_null, genome_type = _node_struct.unpack(buff)
            except:
                raise ValueError("Metadata not a recognized format.")
            md = NodeMetadata(slim_id=slim_id, is_null=is_null, genome_type=genome_type)
    return md

def encode_node(metadata_object):
    '''
    Encodes :class:`NodeMetadata` objects as a bytes object, suitable to be put
    in as metadata for a node. If ``metadata_object`` is ``None``, returns an
    empty bytes object.

    :param NodeMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if metadata_object is None:
        md = b''
    else:
        md = _node_struct.pack(metadata_object.slim_id, metadata_object.is_null,
                               metadata_object.genome_type)
    return md


def extract_node_metadata(tables):
    '''
    Returns an iterator over lists of :class: `NodeMetadata` objects containing
    information about the nodes in the tables.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    metadata = tskit.unpack_bytes(tables.nodes.metadata,
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
    pedigree_id = attr.ib()
    age = attr.ib()
    population = attr.ib()
    sex = attr.ib()
    flags = attr.ib()

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

def decode_individual(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`IndividualTable`. If the buffer is empty, returns None.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`IndividualTable`, as stored by SLiM.
    :rtype IndividualMetadata:
    '''
    if type(buff) == IndividualMetadata:
        md = buff
    else:
        try:
            buff_len = len(buff)
        except:
            raise ValueError("Metadata not a recognized format.")
        if buff_len == 0:
            md = None
        else:
            if buff_len != 24: # 8 + 4 + 4 + 4 + 4:
                raise ValueError("Metadata bytes of incorrect length.")
            try:
                pedigree_id, age, population, sex, flags = _individual_struct.unpack(buff)
            except:
                raise ValueError("Metadata bytes in incorrect format.")
            md = IndividualMetadata(
                        pedigree_id=pedigree_id, age=age, population=population,
                        sex=sex, flags=flags)
    return md

def encode_individual(metadata_object):
    '''
    Encodes :class:`IndividualMetadata` objects as a bytes object, suitable to be put
    in as metadata for an individual.  If ``buff`` is ``None``, returns an
    empty bytes object.

    :param IndividualMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if metadata_object is None:
        md = b''
    else:
        md = _individual_struct.pack(metadata_object.pedigree_id, metadata_object.age,
                                     metadata_object.population, metadata_object.sex,
                                     metadata_object.flags)
    return md


def extract_individual_metadata(tables):
    '''
    Returns an iterator over lists of :class:`IndividualMetadata` objects
    containing information about the individuals in the tables.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    metadata = tskit.unpack_bytes(tables.individuals.metadata,
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

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

    def __iter__(self):
        for a in (self.source_subpop, self.migration_rate):
            yield a

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

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)

def decode_population(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`PopulationTable`. If the buffer is empty, returns None.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`PopulationTable`, as stored by SLiM.
    :rtype PopulationMetadata:
    '''
    if type(buff) == PopulationMetadata:
        md = buff
    else:
        try:
            buff_len = len(buff)
        except:
            raise ValueError("Metadata not a recognized format.")
        if buff_len == 0:
            md = None
        else:
            if buff_len < 88: # 4 + 8 * 10 + 4
                raise ValueError("Metadata bytes of incorrect format.")
            num_migration_records = int((buff_len - 88) / 12) # 4 + 8
            if buff_len != 88 + 12 * num_migration_records:
                raise ValueError("Metadata bytes of incorrect format.")
            struct_string = "<i10di" + "id" * num_migration_records
            try:
                metadata = struct.unpack(struct_string, buff)
            except:
                raise ValueError("Metadata not a recognized format.")
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

            md = PopulationMetadata(slim_id=slim_id, selfing_fraction=selfing_fraction,
                                      female_cloning_fraction=female_cloning_fraction,
                                      male_cloning_fraction=male_cloning_fraction,
                                      sex_ratio=sex_ratio, bounds_x0=bounds_x0, bounds_x1=bounds_x1,
                                      bounds_y0=bounds_y0, bounds_y1=bounds_y1, bounds_z0=bounds_z0,
                                      bounds_z1=bounds_z1, migration_records=migration_records)
    return md

def encode_population(metadata_object):
    '''
    Encodes :class:`PopulationMetadata` objects as a bytes object, suitable to be put
    in as metadata for a population.  If ``buff`` is ``None``, returns an empty
    bytes object.

    :param PopulationMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if metadata_object is None:
        md = b''
    else:
        num_migration_records = len(metadata_object.migration_records)
        for mr in metadata_object.migration_records:
            if type(mr) is not PopulationMigrationMetadata:
                raise ValueError("Migration records should PopulationMigrationMetadata objects.")
        mr_values = [a for b in metadata_object.migration_records for a in b]
        struct_string = "<i10di" + "id" * num_migration_records
        md = struct.pack(struct_string, metadata_object.slim_id,
                         metadata_object.selfing_fraction, metadata_object.female_cloning_fraction,
                         metadata_object.male_cloning_fraction, metadata_object.sex_ratio,
                         metadata_object.bounds_x0, metadata_object.bounds_x1,
                         metadata_object.bounds_y0, metadata_object.bounds_y1,
                         metadata_object.bounds_z0, metadata_object.bounds_z1,
                         num_migration_records, *mr_values)
    return md


def extract_population_metadata(tables):
    '''
    Returns an iterator over lists of :class:`PopulationMetadata` objects
    containing information about the populations in the tables.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    metadata = tskit.unpack_bytes(tables.populations.metadata,
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

