import attr
import struct
import tskit
import json
import warnings

from ._version import *

GENOME_TYPE_AUTOSOME = 0
GENOME_TYPE_X = 1
GENOME_TYPE_Y = 2
INDIVIDUAL_TYPE_HERMAPHRODITE = -1
INDIVIDUAL_TYPE_FEMALE = 0
INDIVIDUAL_TYPE_MALE = 1
INDIVIDUAL_FLAG_MIGRATED = 0x01

# These are copied from slim_globals.h and modified to be python not json (eg false->False)
_raw_slim_metadata_schemas = {
        "tree_sequence" :
        {
            "$schema": "http://json-schema.org/schema#",
            "codec": "json",
            "type": "object",
            "properties": {
                "SLiM" : {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.5",
                    "type": "object",
                    "properties": {
                        "model_type" : {
                            "type": "string",
                            "enum": ["WF", "nonWF"],
                            "description": "The model type used for the last part of this simulation (WF or nonWF)."
                        },
                        "generation" : {
                            "type": "integer",
                            "description": "The 'SLiM generation' counter when this tree sequence was recorded."
                        },
                        "stage" : {
                            "type": "string",
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded."
                        },
                        "file_version" : {
                            "type": "string",
                            "description": "The SLiM 'file format version' of this tree sequence."
                        },
                        "spatial_dimensionality": {
                            "type": "string",
                            "enum": ["", "x", "xy", "xyz"],
                            "description": "The spatial dimensionality of the simulation."
                        },
                        "spatial_periodicity": {
                            "type": "string",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "description": "The spatial periodicity of the simulation."
                        },
                        "separate_sexes": {
                            "type": "boolean",
                            "description": "Whether the simulation had separate sexes."
                        },
                        "nucleotide_based": {
                            "type": "boolean",
                            "description": "Whether the simulation was nucleotide-based."
                        }
                    },
                    "required": ["model_type", "generation", "file_version", "spatial_dimensionality", "spatial_periodicity", "separate_sexes", "nucleotide_based"]
                }
            },
            "required": ["SLiM"],
            "examples": [
                {
                 "SLiM" : {
                     "model_type" : "WF",
                     "generation" : 123,
                     "file_version" : "0.5",
                     "spatial_dimensionality" : "xy",
                     "spatial_periodicity" : "x",
                     "separate_sexes" : True,
                     "nucleotide_based" : False
                     }
                }
            ]
        },
        "edge" : None,
        "site" : None,
        "mutation" : 
        {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for mutation metadata.",
            "codec": "struct",
            "type": "object",
            "properties": {
                "mutation_list": {
                    "type": "array",
                    "noLengthEncodingExhaustBuffer": True,
                    "items": {
                        "type": "object",
                        "properties": {
                            "mutation_type": {
                                "type": "integer",
                                "description": "The index of this mutation's mutationType.",
                                "binaryFormat": "i",
                                "index": 1
                            },
                            "selection_coeff": {
                                "type": "number",
                                "description": "This mutation's selection coefficient.",
                                "binaryFormat": "f",
                                "index": 2
                            },
                            "subpopulation": {
                                "type": "integer",
                                "description": "The ID of the subpopulation this mutation occurred in.",
                                "binaryFormat": "i",
                                "index": 3
                            },
                            "slim_time": {
                                "type": "integer",
                                "description": "The SLiM generation counter when this mutation occurred.",
                                "binaryFormat": "i",
                                "index": 4
                            },
                            "nucleotide": {
                                "type": "integer",
                                "description": "The nucleotide for this mutation (0=A , 1=C , 2=G, 3=T, or -1 for none)",
                                "binaryFormat": "b",
                                "index": 5
                            }
                        },
                        "required": ["mutation_type", "selection_coeff", "subpopulation",
                                     "slim_time", "nucleotide"],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["mutation_list"],
            "additionalProperties": False,
            "examples": [
                {
                 "mutation_list" : [ {
                     "mutation_type" : 1,
                     "selection_coeff" : -0.2,
                     "subpopulation" : 0,
                     "slim_time" : 243,
                     "nucleotide" : 3
                 } ]
                }
            ]
        },
        "node" : 
        {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for node metadata.",
            "codec": "struct",
            "type": ["object", "null"],
            "properties": {
                "slim_id": {
                    "type": "integer",
                    "description": "The 'pedigree ID' of this chromosome in SLiM.",
                    "binaryFormat": "q",
                    "index": 0},
                "is_null": {
                    "type": "boolean",
                    "description": "Whether this node describes a 'null' (non-existant) chromosome.",
                    "binaryFormat": "?",
                    "index": 1},
                "genome_type": {
                    "type": "integer",
                    "description": "The 'type' of this genome (0 for autosome, 1 for X, 2 for Y).",
                    "binaryFormat":
                    "B", "index": 2}
            },
            "required": ["slim_id", "is_null", "genome_type"],
            "additionalProperties": False,
            "examples": [
                {
                "slim_id" : 123,
                "is_null" : False,
                "genome_type" : 0
                }
            ]
        },
        "individual" :
        {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for individual metadata.",
            "codec": "struct",
            "type": "object",
            "properties": {
                "pedigree_id": {
                    "type": "integer",
                    "description": "The 'pedigree ID' of this individual in SLiM.",
                    "binaryFormat": "q",
                    "index": 1
                },
                "age": {
                    "type": "integer",
                    "description": "The age of this individual, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                    "binaryFormat": "i",
                    "index": 2
                },
                "subpopulation": {
                    "type": "integer",
                    "description": "The ID of the subpopulation the individual was part of, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                    "binaryFormat": "i",
                    "index": 3
                },
                "sex": {
                    "type": "integer",
                    "description": "The sex of the individual (0 for female, 1 for male, -1 for hermaphrodite).",
                    "binaryFormat": "i",
                    "index": 4
                },
                "flags": {
                    "type": "integer",
                    "description": "Other information about the individual: see 'flags'.",
                    "binaryFormat": "I",
                    "index": 5
                }
            },
            "required": ["pedigree_id", "age", "subpopulation", "sex", "flags"],
            "additionalProperties": False,
            "flags": {
                 "SLIM_INDIVIDUAL_METADATA_MIGRATED": {
                     "value": 1,
                     "description": "Whether this individual was a migrant, either in the generation when the tree sequence was written out (if the individual was alive then), or in the generation of the last time they were Remembered (if not)."
                 }
            },
            "examples": [
                {
                "pedigree_id" : 123,
                "age" : -1,
                "subpopulation" : 0,
                "sex" : 0,
                "flags" : 0
                }
            ]
        },
        "population" :
        {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for population metadata.",
            "codec": "struct",
            "type": ["object", "null"],
            "properties": {
                "slim_id": {
                    "type": "integer",
                    "description": "The ID of this population in SLiM. Note that this is called a 'subpopulation' in SLiM.",
                    "binaryFormat": "i",
                    "index": 1},
                "selfing_fraction": {
                    "type": "number",
                    "description": "The frequency with which individuals in this subpopulation self (for WF models).",
                    "binaryFormat": "d",
                    "index": 2},
                "female_cloning_fraction": {
                    "type": "number",
                    "description": "The frequency with which females in this subpopulation reproduce clonally (for WF models).",
                    "binaryFormat": "d",
                    "index": 3},
                "male_cloning_fraction": {
                    "type": "number",
                    "description": "The frequency with which males in this subpopulation reproduce clonally (for WF models).",
                    "binaryFormat": "d",
                    "index": 4},
                "sex_ratio": {
                    "type": "number",
                    "description": "This subpopulation's sex ratio (for WF models).",
                    "binaryFormat": "d",
                    "index": 5},
                "bounds_x0": {
                    "type": "number",
                    "description": "The minimum x-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 6},
                "bounds_x1": {
                    "type": "number",
                    "description": "The maximum x-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 7},
                "bounds_y0": {
                    "type": "number",
                    "description": "The minimum y-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 8},
                "bounds_y1": {
                    "type": "number",
                    "description": "The maximum y-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 9},
                "bounds_z0": {
                    "type": "number",
                    "description": "The minimum z-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 10},
                "bounds_z1": {
                    "type": "number",
                    "description": "The maximum z-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 11},
                "migration_records": {
                    "type": "array",
                    "index": 13,
                    "arrayLengthFormat": "I",
                    "items": {
                        "type": "object",
                        "properties": {
                            "source_subpop": {
                                "type": "integer",
                                "description": "The ID of the subpopulation migrants come from (in WF models).",
                                "binaryFormat": "i",
                                "index": 1},
                            "migration_rate": {
                                "type": "number",
                                "description": "The fraction of children in this subpopulation that are composed of 'migrants' from the source subpopulation (in WF models).",
                                "binaryFormat": "d",
                                "index": 2}
                        },
                        "required": ["source_subpop", "migration_rate"],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["slim_id", "selfing_fraction", "female_cloning_fraction",
                         "male_cloning_fraction", "sex_ratio", "bounds_x0", "bounds_x1",
                         "bounds_y0", "bounds_y1", "bounds_z0", "bounds_z1", "migration_records"],
            "additionalProperties": False,
            "examples": [
                {
                   "slim_id": 2,
                   "selfing_fraction": 0.5,
                   "female_cloning_fraction": 0.25,
                   "male_cloning_fraction": 0.0,
                   "sex_ratio": 0.5,
                   "bounds_x0": 0.0,
                   "bounds_x1": 100.0,
                   "bounds_y0": 0.0,
                   "bounds_y1": 100.0,
                   "bounds_z0": 0.0,
                   "bounds_z1": 100.0,
                   "migration_records": [
                       {"source_subpop": 1,
                        "migration_rate": 0.9 },
                       {"source_subpop": 2,
                        "migration_rate": 0.1 }
                   ]
                }
            ]
        },
    }


slim_metadata_schemas = {k: tskit.MetadataSchema(_raw_slim_metadata_schemas[k])
        for k in _raw_slim_metadata_schemas}


def default_slim_metadata(name):
    """
    Returns default metadata of type ``name``, where ``name`` is one of
    "tree_sequence", "edge", "site", "mutation", "node", "individual", or
    "population".

    :param str name: The type of metadata requested.
    :rtype dict:
    """
    if name == "tree_sequence":
        out = {
             "SLiM" : {
                 "model_type" : "nonWF",
                 "generation" : 1,
                 "file_version" : slim_file_version,
                 "spatial_dimensionality" : "",
                 "spatial_periodicity" : "",
                 "separate_sexes" : False,
                 "nucleotide_based" : False,
                 "stage" : "late"
             }
         }
    elif name == "edge":
        out = None
    elif name == "site":
        out = None
    elif name == "mutation":
        out = {
            "mutation_list": []
        }
    elif name == "node":
        out = {
            "slim_id": tskit.NULL,
            "is_null": False,
            "genome_type": 0,
        }
    elif name == "individual":
        out = {
            "pedigree_id": tskit.NULL,
            "age": -1,
            "subpopulation": tskit.NULL,
            "sex": -1,
            "flags": 0,
        }
    elif name == "population":
        out = {
            "slim_id": tskit.NULL,
            "selfing_fraction": 0.0,
            "female_cloning_fraction": 0.0,
            "male_cloning_fraction": 0.0,
            "sex_ratio": 0.0,
            "bounds_x0": 0.0,
            "bounds_x1": 0.0,
            "bounds_y0": 0.0,
            "bounds_y1": 0.0,
            "bounds_z0": 0.0,
            "bounds_z1": 0.0,
            "migration_records": []
        }
    else:
        raise ValueError(
            "Unknown metadata request: name should be one of 'tree_sequence', "
            "'edge', 'site', 'mutation', 'node', 'individual', or 'population'.")
    return out


###########
# Top-level, a.k.a., tree sequence metadata
###########

def set_tree_sequence_metadata(tables,
        model_type,
        generation,
        spatial_dimensionality='',
        spatial_periodicity='',
        separate_sexes=False,
        nucleotide_based=False,
        stage='late',
        file_version=None):
    if file_version is None:
        file_version = slim_file_version
    if isinstance(tables.metadata, bytes):
        if len(tables.metadata) > 0:
            raise ValueError("Tree sequence has top-level metadata but no schema: this is a problem "
                             "since pyslim is trying to add to the metadata.")
        schema_dict = slim_metadata_schemas['tree_sequence'].schema
        metadata_dict = {}
    else:
        # we need to keep other keys in the metadata (and schema) if there are any
        schema_dict = tables.metadata_schema.schema
        metadata_dict = tables.metadata
    assert(schema_dict['codec'] == 'json')
    assert(schema_dict['type'] == 'object')
    if "properties" not in schema_dict:
        schema_dict["properties"] = {}
    schema_dict['properties']['SLiM'] = slim_metadata_schemas['tree_sequence'].schema['properties']['SLiM']
    tables.metadata_schema = tskit.MetadataSchema(schema_dict)
    metadata_dict['SLiM'] = {
            "model_type": model_type,
            "generation": generation,
            "file_version": file_version,
            "spatial_dimensionality": spatial_dimensionality,
            "spatial_periodicity": spatial_periodicity,
            "separate_sexes": separate_sexes,
            "nucleotide_based": nucleotide_based,
            "stage": stage,
            }
    tables.metadata = metadata_dict
    _set_metadata_schemas(tables)


def _set_metadata_schemas(tables):
    tables.edges.metadata_schema = slim_metadata_schemas['edge']
    tables.sites.metadata_schema = slim_metadata_schemas['site']
    tables.mutations.metadata_schema = slim_metadata_schemas['mutation']
    tables.nodes.metadata_schema = slim_metadata_schemas['node']
    tables.individuals.metadata_schema = slim_metadata_schemas['individual']
    tables.populations.metadata_schema = slim_metadata_schemas['population']



################################
# Old-style metadata:

def _deprecation_warning(name):
    warnings.warn(f"This method ({name}) will dissappear at some point, "
                  "along with all other old-style metadata tools: "
                  "see `the documentation <https://pyslim.readthedocs.io/en/latest/metadata.html#legacy-metadata>`_"
                  "for more details.", FutureWarning)


def _legacy_error(name):
    raise ValueError(
            f"{name} recieved a dict instead of bytes:"
            "it looks like you're trying to use pyslim legacy metadata tools: "
            "see `the documentation <https://pyslim.readthedocs.io/en/latest/metadata.html#legacy-metadata>`_"
            "for how to update your script.")


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

    def asdict(self):
        out = {
                "mutation_type" : int(self.mutation_type),
                "selection_coeff" : float(self.selection_coeff),
                "subpopulation" : int(self.population),
                "slim_time" : int(self.slim_time),
                "nucleotide" : int(self.nucleotide),
              }
        return out

    @classmethod
    def fromdict(cls, md):
        return cls(mutation_type=md['mutation_type'],
                   selection_coeff=md['selection_coeff'],
                   population=md['subpopulation'],
                   slim_time=md['slim_time'],
                   nucleotide=md['nucleotide'])


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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`MutationTable`, as stored by SLiM.
    :rtype list:
    '''
    if isinstance(buff, dict):
        _legacy_error("decode_mutation")
    _deprecation_warning("decode_mutation")
    if isinstance(buff, list) and (len(buff) == 0
            or isinstance(buff[0], MutationMetadata)):
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
    
    **Only used in upgrading old tables.**

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
    mut_schema = slim_metadata_schemas['mutation']
    recoded = mut_schema.encode_row(
        {'mutation_list': [mm.asdict() for mm in mut_structs]})
    return recoded


def encode_mutation(metadata_object):
    '''
    Encodes the list of :class:`MutationMetadata` objects as a bytes object,
    suitable to be put in as metadata for a mutation.  Unlike other ``encode_``
    functions, this takes a *list* rather than a single value, thanks to
    stacking of SLiM mutations.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param MutationMetadata metadata_object: The list of
        :class:`MutationMetadata` objects to be encoded.
    :rtype bytes:
    '''
    if isinstance(metadata_object, dict):
        _legacy_error("encode_mutation")
    _deprecation_warning("encode_mutation")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    _deprecation_warning("extract_mutation_metadata")
    metadata = tskit.unpack_bytes(tables.mutations.metadata,
                                    tables.mutations.metadata_offset)
    for mut in tables.mutations:
        yield [MutationMetadata.fromdict(mm) for mm in mut.metadata['mutation_list']]


def annotate_mutation_metadata(tables, metadata):
    '''
    Revise the mutation table in place so that the metadata column is given by
    applying `encode_mutation()` to the sources given.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of (lists of MutationMetadata) or None objects
    '''
    _deprecation_warning("annotate_mutation_metadata")
    if len(metadata) != tables.mutations.num_rows:
        raise ValueError("annotate mutations: metadata not the same length as the table.")
    orig_mutations = tables.mutations.copy()
    tables.mutations.clear()
    for mut, md in zip(orig_mutations, metadata):
        if md is not None:
            md = {'mutation_list': [mm.asdict() for mm in md]}
        tables.mutations.add_row(site=mut.site, node=mut.node, derived_state=mut.derived_state,
                                 parent=mut.parent, time=mut.time, metadata=md)


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
        return self.asdict() == other.asdict()

    def __ne__(self, other):
        return not self.__eq__(other)

    def asdict(self):
        out = {
                "slim_id" : int(self.slim_id),
                "is_null" : bool(self.is_null),
                "genome_type" : int(self.genome_type),
              }
        return out

    @classmethod
    def fromdict(cls, md):
        if md is None:
            out = None
        else:
            out = cls(**md)
        return out


def decode_node(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`NodeTable`.  If the buffer is empty, returns None.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`NodeTable`, as stored by SLiM.
    :rtype NodeMetadata:
    '''
    if isinstance(buff, dict):
        _legacy_error("decode_node")
    _deprecation_warning("decode_node")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param NodeMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if isinstance(metadata_object, dict):
        _legacy_error("encode_node")
    _deprecation_warning("encode_node")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    _deprecation_warning("extract_node_metadata")
    for n in tables.nodes:
        yield NodeMetadata.fromdict(n.metadata)


def annotate_node_metadata(tables, metadata):
    '''
    Modify the NodeTable so that the metadata
    column is given by applying `encode_node()` to the sources given.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of NodeMetadata or None objects
    '''
    _deprecation_warning("annotate_node_metadata")
    if len(metadata) != tables.nodes.num_rows:
        raise ValueError("annotate nodes: metadata not the same length as the table.")
    orig_nodes = tables.nodes.copy()
    tables.nodes.clear()
    for node, md in zip(orig_nodes, metadata):
        tables.nodes.add_row(flags=node.flags, time=node.time, population=node.population,
                             individual=node.individual, metadata=md.asdict())


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
        return self.asdict() == other.asdict()

    def __ne__(self, other):
        return not self.__eq__(other)

    def asdict(self):
        out = {
                "pedigree_id": int(self.pedigree_id),
                "age": int(self.age),
                "subpopulation": int(self.population),
                "sex": int(self.sex),
                "flags": int(self.flags),
                }
        return out

    @classmethod
    def fromdict(cls, md):
        return cls(pedigree_id=md['pedigree_id'],
                   age=md['age'],
                   population=md['subpopulation'],
                   sex=md['sex'],
                   flags=md['flags'])


def decode_individual(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`IndividualTable`. If the buffer is empty, returns None.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`IndividualTable`, as stored by SLiM.
    :rtype IndividualMetadata:
    '''
    if isinstance(buff, dict):
        _legacy_error("decode_individual")
    _deprecation_warning("decode_individual")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param IndividualMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if isinstance(metadata_object, dict):
        _legacy_error("encode_individual")
    _deprecation_warning("encode_individual")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    _deprecation_warning("extract_individual_metadata")
    for ind in tables.individuals:
        yield IndividualMetadata.fromdict(ind.metadata)


def annotate_individual_metadata(tables, metadata):
    '''
    Modify the IndividualTable so that the metadata
    column is given by applying `encode_individual()` to the sources given.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of IndividualMetadata or None objects
    '''
    _deprecation_warning("annotate_individual_metadata")
    if len(metadata) != tables.individuals.num_rows:
        raise ValueError("annotate individuals: metadata not the same length as the table.")
    orig_individuals = tables.individuals.copy()
    tables.individuals.clear()
    for ind, md in zip(orig_individuals, metadata):
        if md is None:
            md = default_slim_metadata('individual')
        tables.individuals.add_row(flags=ind.flags, location=ind.location, metadata=md.asdict())

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
        return self.asdict() == other.asdict()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __iter__(self):
        for a in (self.source_subpop, self.migration_rate):
            yield a

    def asdict(self):
        out = {
                "source_subpop" : int(self.source_subpop),
                "migration_rate" : float(self.migration_rate),
              }
        return out

    @classmethod
    def fromdict(cls, md):
        return cls(**md)


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
        return self.asdict() == other.asdict()

    def __ne__(self, other):
        return not self.__eq__(other)

    def asdict(self):
        out = {
                "slim_id" : int(self.slim_id),
                "selfing_fraction" : float(self.selfing_fraction),
                "female_cloning_fraction" : float(self.female_cloning_fraction),
                "male_cloning_fraction" : float(self.male_cloning_fraction),
                "sex_ratio" : float(self.sex_ratio),
                "bounds_x0" : float(self.bounds_x0),
                "bounds_x1" : float(self.bounds_x1),
                "bounds_y0" : float(self.bounds_y0),
                "bounds_y1" : float(self.bounds_y1),
                "bounds_z0" : float(self.bounds_z0),
                "bounds_z1" : float(self.bounds_z1),
                "migration_records" : [u.asdict() for u in self.migration_records],
              }
        return out

    @classmethod
    def fromdict(cls, md):
        if md is None:
            return None
        else:
            md['migration_records'] = [
                    PopulationMigrationMetadata.fromdict(u) for u in md['migration_records']
                    ]
            return cls(**md)


def decode_population(buff):
    '''
    Extracts the information stored in binary by SLiM in the ``metadata``
    column of a :class:`PopulationTable`. If the buffer is empty, returns None.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param bytes buff: The ``metadata`` entry of a row of a
        :class:`PopulationTable`, as stored by SLiM.
    :rtype PopulationMetadata:
    '''
    if isinstance(buff, dict):
        _legacy_error("decode_population")
    _deprecation_warning("decode_population")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param PopulationMetadata metadata_object: The object to be encoded.
    :rtype bytes:
    '''
    if isinstance(metadata_object, dict):
        _legacy_error("encode_population")
    _deprecation_warning("encode_population")
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

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: The tables, as produced by SLiM.
    '''
    _deprecation_warning("extract_population_metadata")
    for pop in tables.populations:
        yield PopulationMetadata.fromdict(pop.metadata)


def annotate_population_metadata(tables, metadata):
    '''
    Modify the PopulationTable so that the metadata
    column is given by applying `encode_population()` to the sources given.
    This entirely removes existing information in the Population table.

    .. warning::

        This method is deprecated, since metadata handling has been taken over
        by tskit. It will dissappear at some point in the future.

    :param TableCollection tables: a table collection to be modified
    :param iterable metadata: a list of objects, each PopulationMetadata or None
    '''
    _deprecation_warning("annotate_population_metadata")
    if len(metadata) != tables.populations.num_rows:
        raise ValueError("annotate populations: metadata not the same length as the table.")
    tables.populations.clear()
    for md in metadata:
        if md is not None:
            md = md.asdict()
        tables.populations.add_row(metadata=md)

