import copy
import json
import warnings

import numpy as np
import tskit

from ._version import *  # noqa F403
from .provenance import get_environment, slim_provenance_version

_INDIVIDUAL_FLAG_MIGRATED_DEPRECATED = np.uint32(1 << 0)


def is_vacant_num_bytes(num_chromosomes):
    """
    TODO document
    """
    # see _is_chrom_vacant
    return int((num_chromosomes + 7) / 8)


# def _isvacant_values(vacancy):
#     """
#     Returns the correct is_vacant metadata for the given vacancy pattern.
#
#     :param list vacancy: A list of booleans.
#     """
#     out = [0 for _ in is_vacant_num_bytes(len(vacancy))]
#     powers = [1, 2, 4, 8, 16, 32, 64, 128]
#     for k, v in enumerate(vacancy):
#         if v:
#             n = k % 8
#             m = int(k / 8)
#             out[m] += powers[n]
#     return out


# These are copied from slim_globals.cpp, modified to be python not json
# by changing false->False and true->True, then printed with
# json.dump(..., indent=True), then lightly edited.

_raw_slim_metadata_schemas = {
    "tree_sequence": {
        "$schema": "http://json-schema.org/schema#",
        "codec": "json+struct",
        "json": {
            "codec": "json",
            "description": "SLiM schema for JSON top-level metadata.",
            "examples": [
                {
                    "SLiM": {
                        "chromosomes": [
                            {"id": 1, "name": "autosome_1", "symbol": "1", "type": "A"},
                            {"id": 35, "name": "mtDNA", "symbol": "MT", "type": "HF"},
                        ],
                        "cycle": 123,
                        "description": "foxes on Catalina island",
                        "file_version": "1.0",
                        "model_type": "WF",
                        "name": "fox",
                        "nucleotide_based": False,
                        "separate_sexes": True,
                        "spatial_dimensionality": "xy",
                        "spatial_periodicity": "x",
                        "this_chromosome": {
                            "id": 1,
                            "index": 0,
                            "name": "autosome_1",
                            "symbol": "1",
                            "type": "A",
                        },
                        "tick": 123,
                        "traits": [
                            {"index": 0, "name": "simT", "type": "multiplicative"}
                        ],
                    }
                }
            ],
            "properties": {
                "SLiM": {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 1.0",
                    "properties": {
                        "chromosomes": {
                            "description": "The chromosomes represented by the collection of tree sequences, of which this tree sequence is one member.",
                            "items": {
                                "properties": {
                                    "id": {
                                        "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
                                        "type": "integer",
                                    },
                                    "name": {
                                        "description": "A user-specified name for the chromosome, such as an accession identifier.",
                                        "type": "string",
                                    },
                                    "symbol": {
                                        "description": 'A short string symbol for the chromosome, unique within this set of tree sequences, such as "1" or "MT".',
                                        "type": "string",
                                    },
                                    "type": {
                                        "description": "The type of chromosome, as specified by SLiM.",
                                        "type": "string",
                                    },
                                },
                                "required": ["id", "symbol", "type"],
                                "type": "object",
                            },
                            "type": "array",
                        },
                        "cycle": {
                            "description": "The 'SLiM cycle' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "description": {
                            "description": "A user-configurable description of the species represented by this tree sequence.",
                            "type": "string",
                        },
                        "file_version": {
                            "description": "The SLiM 'file format version' of this tree sequence.",
                            "type": "string",
                        },
                        "model_type": {
                            "description": "The model type used for the last part of this simulation (WF or nonWF).",
                            "enum": ["WF", "nonWF"],
                            "type": "string",
                        },
                        "name": {
                            "description": "The SLiM species name represented by this tree sequence.",
                            "type": "string",
                        },
                        "nucleotide_based": {
                            "description": "Whether the simulation was nucleotide-based.",
                            "type": "boolean",
                        },
                        "separate_sexes": {
                            "description": "Whether the simulation had separate sexes.",
                            "type": "boolean",
                        },
                        "spatial_dimensionality": {
                            "description": "The spatial dimensionality of the simulation.",
                            "enum": ["", "x", "xy", "xyz"],
                            "type": "string",
                        },
                        "spatial_periodicity": {
                            "description": "The spatial periodicity of the simulation.",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "type": "string",
                        },
                        "stage": {
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
                            "type": "string",
                        },
                        "this_chromosome": {
                            "description": "The chromosome represented by the tree sequence in this file.",
                            "properties": {
                                "id": {
                                    "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
                                    "type": "integer",
                                },
                                "index": {
                                    "description": "The (zero-based) index of this chromosome in the chromosomes metadata array (if present), which should match the information given here.",
                                    "type": "integer",
                                },
                                "name": {
                                    "description": "A user-specified name for the chromosome, such as an accession identifier.",
                                    "type": "string",
                                },
                                "symbol": {
                                    "description": 'A short string symbol for the chromosome, unique within this set of tree sequences, such as "1" or "MT".',
                                    "type": "string",
                                },
                                "type": {
                                    "description": "The type of chromosome, as specified by SLiM.",
                                    "type": "string",
                                },
                            },
                            "required": ["id", "index", "symbol", "type"],
                            "type": "object",
                        },
                        "tick": {
                            "description": "The 'SLiM tick' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "traits": {
                            "description": "The traits defined for this tree sequence; each mutation and individual will have per-trait metadata.",
                            "items": {
                                "properties": {
                                    "baselineAccumulation": {
                                        "description": "Whether the baseline offset includes accumulated effects from fixed (substituted) mutations.",
                                        "type": "boolean",
                                    },
                                    "baselineOffsetFromUser": {
                                        "type": "number",
                                        "description": "The from-user component of the baseline offset of the trait.",
                                    },
                                    "baselineOffsetFromSubstitutions": {
                                        "type": "number",
                                        "description": "The from-substitutions component of the baseline offset of the trait.",
                                    },
                                    "directFitnessEffect": {
                                        "description": "Whether the trait's effects are used directly as fitness effects.",
                                        "type": "boolean",
                                    },
                                    "index": {
                                        "description": "The integer index for the trait; indices must be sequential starting from zero.",
                                        "type": "integer",
                                    },
                                    "individualOffsetMean": {
                                        "description": "The mean of the trait's individual offset distribution (which might or might not be used).",
                                        "type": "number",
                                    },
                                    "individualOffsetSD": {
                                        "description": "The standard deviation of the trait's individual offset distribution (which might or might not be used).",
                                        "type": "number",
                                    },
                                    "name": {
                                        "description": "The string name for the trait.",
                                        "type": "string",
                                    },
                                    "type": {
                                        "description": "The type of the trait; this must be 'additive', 'multiplicative', or 'logistic'.",
                                        "enum": [
                                            "additive",
                                            "multiplicative",
                                            "logistic",
                                        ],
                                        "type": "string",
                                    },
                                },
                                "required": ["index", "name", "type"],
                                "type": "object",
                            },
                            "type": "array",
                        },
                    },
                    "required": [
                        "model_type",
                        "tick",
                        "file_version",
                        "spatial_dimensionality",
                        "spatial_periodicity",
                        "this_chromosome",
                        "separate_sexes",
                        "nucleotide_based",
                        "traits",
                    ],
                    "type": "object",
                }
            },
            "required": ["SLiM"],
            "type": "object",
        },
        "struct": {
            "codec": "struct",
            "description": "SLiM schema for binary top-level metadata.",
            "properties": {
                "SLiM_mutation_list": {
                    "arrayLengthFormat": "Q",
                    "items": {
                        "additionalProperties": False,
                        "properties": {
                            "mutation_id": {
                                "binaryFormat": "q",
                                "description": "The SLiM mutation ID for this mutation.",
                                "index": 1,
                                "type": "integer",
                            },
                            "mutation_type": {
                                "binaryFormat": "i",
                                "description": "The id of this mutation's mutationType.",
                                "index": 2,
                                "type": "integer",
                            },
                            "nucleotide": {
                                "binaryFormat": "b",
                                "description": "The nucleotide for this mutation (0=A , 1=C , 2=G, 3=T, or -1 for none)",
                                "index": 5,
                                "type": "integer",
                            },
                            "padding": {
                                "binaryFormat": "3x",
                                "description": "Padding bytes for alignment",
                                "index": 6,
                                "type": "null",
                            },
                            "per_trait": {
                                "index": 7,
                                "items": {
                                    "additionalProperties": False,
                                    "properties": {
                                        "dominance": {
                                            "binaryFormat": "f",
                                            "description": "The dominance coefficient for this trait.",
                                            "index": 2,
                                            "type": "number",
                                        },
                                        "effect_size": {
                                            "binaryFormat": "f",
                                            "description": "The effect size for this trait.",
                                            "index": 1,
                                            "type": "number",
                                        },
                                        "hemizygous_dominance": {
                                            "binaryFormat": "f",
                                            "description": "The hemizygous dominance coefficient for this trait.",
                                            "index": 3,
                                            "type": "number",
                                        },
                                    },
                                    "required": [
                                        "dominance",
                                        "effect_size",
                                        "hemizygous_dominance",
                                    ],
                                    "type": "object",
                                },
                                "length": 1,  # NOTE this may need to be changed to match the number of traits!
                                "type": "array",
                            },
                            "slim_time": {
                                "binaryFormat": "i",
                                "description": "The SLiM tick counter when this mutation occurred.",
                                "index": 4,
                                "type": "integer",
                            },
                            "subpopulation": {
                                "binaryFormat": "i",
                                "description": "The ID of the subpopulation this mutation occurred in.",
                                "index": 3,
                                "type": "integer",
                            },
                        },
                        "required": [
                            "mutation_id",
                            "mutation_type",
                            "slim_time",
                            "subpopulation",
                            "nucleotide",
                            "per_trait",
                        ],
                        "type": "object",
                    },
                    "type": "array",
                }
            },
            "required": ["SLiM_mutation_list"],
            "type": "object",
        },
    },
    "edge": None,
    "site": None,
    "mutation": {
        "$schema": "http://json-schema.org/schema#",
        "additionalProperties": False,
        "codec": "struct",
        "type": "object",
        "description": "SLiM schema for representing binary derived state data in mutation metadata (the metadata for each unique SLiM mutation is stored in top-level metadata).",
        "examples": [{"derived_states": [0, 1, 17]}],
        "properties": {
            "derived_states": {
                "index": 1,
                "type": "array",
                "noLengthEncodingExhaustBuffer": True,
                "description": "An array of SLiM mutation IDs (int64t), representing the (stacked) mutations contained by the derived state for the mutation.",
                "items": {"binaryFormat": "q", "type": "number"},
            }
        },
        "required": ["derived_states"],
    },
    "node": {
        "$schema": "http://json-schema.org/schema#",
        "additionalProperties": False,
        "codec": "struct",
        "description": "SLiM schema for node metadata.",
        "examples": [{"slim_id": 123, "is_vacant": 0}],
        "properties": {
            "slim_id": {
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of the haplosomes associated with this node in SLiM.",
                "index": 1,
                "type": "integer",
            },
            "is_vacant": {
                "description": "A vector of byte (uint8_t) values, with each bit representing whether the node represents a vacant position, either unused or a null haplosome (1), or a non-null haplosome (0), in the corresponding chromosome. This field encodes vacancy for all of the chromosomes in the model, not just the chromosome represented in this file (so that the node table is identical across all chromosomes for a multi-chromosome model). Each chromosome receives one bit here; there are two node table entries per individual, used for the two haplosomes of every chromosome, so only one bit is needed in each entry (making two bits total per chromosome, across the two node table entries). The least significant bit of the first byte is used first (for one haplosome of the first chromosome); the most significant bit of the last byte is used last. The number of bytes present in this field is indicated by this schema's 'binaryFormat' field, which is variable (!), and can also be deduced from the number of chromosomes in the model as given in the top-level 'chromosomes' metadata key, which should always be present if this metadata is present.",
                "index": 2,
                "type": "array",
                "length": 1,  # MAY NEED TO BE CHANGED (in SLiM code is "%d")
                "items": {"type": "number", "binaryFormat": "B"},
            },
        },
        "required": ["slim_id", "is_vacant"],
        "type": ["object", "null"],
    },
    "individual": {
        "$schema": "http://json-schema.org/schema#",
        "codec": "struct",
        "type": "object",
        "description": "SLiM schema for individual metadata.",
        "examples": [
            {
                "age": -1,
                "pedigree_id": 123,
                "pedigree_p1": 12,
                "pedigree_p2": 23,
                "per_trait": [{"offset": 1.0, "phenotype": 1.1}],
                "sex": 0,
                "subpopulation": 0,
                "tag": 1,
                "tagF": 5.5,
                "tagL0_set": True,
                "tagL0": True,
                "tagL1_set": True,
                "tagL1": False,
                "tagL2_set": False,
                "tagL2": False,
                "tagL3_set": False,
                "tagL3": False,
                "tagL4_set": False,
                "tagL4": False,
            }
        ],
        "properties": {
            "pedigree_id": {
                "index": 1,
                "type": "integer",
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual in SLiM.",
            },
            "pedigree_p1": {
                "index": 2,
                "type": "integer",
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual's first parent in SLiM.",
            },
            "pedigree_p2": {
                "index": 3,
                "type": "integer",
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual's second parent in SLiM.",
            },
            "tag": {
                "index": 4,
                "type": "integer",
                "binaryFormat": "q",
                "description": "The `tag` property of this individual; INT64_MIN if unset.",
            },
            "tagF": {
                "index": 5,
                "type": "number",
                "binaryFormat": "d",
                "description": "The `tagF` property of this individual; -DBL_MAX if unset.",
            },
            "age": {
                "index": 6,
                "type": "integer",
                "binaryFormat": "i",
                "description": "The age of this individual, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
            },
            "subpopulation": {
                "index": 7,
                "type": "integer",
                "binaryFormat": "i",
                "description": "The ID of the subpopulation the individual was part of, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
            },
            "sex": {
                "index": 8,
                "type": "integer",
                "binaryFormat": "h",
                "description": "The sex of the individual (0 for female, 1 for male, -1 for hermaphrodite).",
            },
            "tagL0_set": {
                "index": 9,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "A flag indicating whether the `tagL0` property is set; if false, accessing `tagL0` is invalid.",
            },
            "tagL0": {
                "index": 10,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "The `tagL0` property of this individual; only valid if `tagL0_set` is true.",
            },
            "tagL1_set": {
                "index": 11,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "A flag indicating whether the `tagL1` property is set; if false, accessing `tagL1` is invalid.",
            },
            "tagL1": {
                "index": 12,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "The `tagL1` property of this individual; only valid if `tagL1_set` is true.",
            },
            "tagL2_set": {
                "index": 13,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "A flag indicating whether the `tagL2` property is set; if false, accessing `tagL2` is invalid.",
            },
            "tagL2": {
                "index": 14,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "The `tagL2` property of this individual; only valid if `tagL2_set` is true.",
            },
            "tagL3_set": {
                "index": 15,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "A flag indicating whether the `tagL3` property is set; if false, accessing `tagL3` is invalid.",
            },
            "tagL3": {
                "index": 16,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "The `tagL3` property of this individual; only valid if `tagL3_set` is true.",
            },
            "tagL4_set": {
                "index": 17,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "A flag indicating whether the `tagL4` property is set; if false, accessing `tagL4` is invalid.",
            },
            "tagL4": {
                "index": 18,
                "type": "boolean",
                "binaryFormat": "?",
                "description": "The `tagL4` property of this individual; only valid if `tagL4_set` is true.",
            },
            "per_trait": {
                "index": 19,
                "type": "array",
                "length": 1,  # MAY NEED TO BE CHANGED (in SLiM code is "%d")
                "items": {
                    "additionalProperties": False,
                    "properties": {
                        "phenotype": {
                            "index": 1,
                            "type": "number",
                            "binaryFormat": "d",
                            "description": "The phenotype for this trait.",
                        },
                        "offset": {
                            "index": 2,
                            "type": "number",
                            "binaryFormat": "d",
                            "description": "The individual offset for this trait.",
                        },
                    },
                    "required": ["offset", "phenotype"],
                    "type": "object",
                },
            },
        },
        "additionalProperties": False,
        "required": [
            "pedigree_id",
            "pedigree_p1",
            "pedigree_p2",
            "tag",
            "tagF",
            "age",
            "subpopulation",
            "sex",
            "tagL0_set",
            "tagL0",
            "tagL1_set",
            "tagL1",
            "tagL2_set",
            "tagL2",
            "tagL3_set",
            "tagL3",
            "tagL4_set",
            "tagL4",
            "per_trait",
        ],
    },
    "population": {
        "$schema": "http://json-schema.org/schema#",
        "additionalProperties": True,
        "codec": "json",
        "description": "SLiM schema for population metadata.",
        "examples": [
            {
                "bounds_x0": 0.0,
                "bounds_x1": 100.0,
                "bounds_y0": 0.0,
                "bounds_y1": 100.0,
                "female_cloning_fraction": 0.25,
                "male_cloning_fraction": 0.0,
                "migration_records": [
                    {"migration_rate": 0.9, "source_subpop": 1},
                    {"migration_rate": 0.1, "source_subpop": 2},
                ],
                "selfing_fraction": 0.5,
                "sex_ratio": 0.5,
                "slim_id": 2,
                "name": "p2",
            }
        ],
        "properties": {
            "bounds_x0": {
                "description": "The minimum x-coordinate in this subpopulation.",
                "type": "number",
            },
            "bounds_x1": {
                "description": "The maximum x-coordinate in this subpopulation.",
                "type": "number",
            },
            "bounds_y0": {
                "description": "The minimum y-coordinate in this subpopulation.",
                "type": "number",
            },
            "bounds_y1": {
                "description": "The maximum y-coordinate in this subpopulation.",
                "type": "number",
            },
            "bounds_z0": {
                "description": "The minimum z-coordinate in this subpopulation.",
                "type": "number",
            },
            "bounds_z1": {
                "description": "The maximum z-coordinate in this subpopulation.",
                "type": "number",
            },
            "description": {
                "description": "A description of this subpopulation.",
                "type": "string",
            },
            "female_cloning_fraction": {
                "description": "The frequency with which females in this subpopulation reproduce clonally (for WF models).",
                "type": "number",
            },
            "male_cloning_fraction": {
                "description": "The frequency with which males in this subpopulation reproduce clonally (for WF models).",
                "type": "number",
            },
            "migration_records": {
                "items": {
                    "properties": {
                        "migration_rate": {
                            "description": "The fraction of children in this subpopulation that are composed of 'migrants' from the source subpopulation (in WF models).",
                            "type": "number",
                        },
                        "source_subpop": {
                            "description": "The ID of the subpopulation migrants come from (in WF models).",
                            "type": "integer",
                        },
                    },
                    "required": ["source_subpop", "migration_rate"],
                    "type": "object",
                },
                "type": "array",
            },
            "name": {
                "description": "A human-readable name for this subpopulation.",
                "type": "string",
            },
            "selfing_fraction": {
                "description": "The frequency with which individuals in this subpopulation self (for WF models).",
                "type": "number",
            },
            "sex_ratio": {
                "description": "This subpopulation's sex ratio (for WF models).",
                "type": "number",
            },
            "slim_id": {
                "description": "The ID of this population in SLiM. Note that this is called a 'subpopulation' in SLiM.",
                "type": "integer",
            },
        },
        "required": [],
        "type": ["object", "null"],
    },
}


def slim_tree_sequence_metadata_schema(num_traits=1):
    """
    The top-level metadata schema depends on the number of traits, and
    {data}`.slim_metadata_schemas`
    returns the schema for a single-trait simulation. This function
    returns the correct schema for a simulation with arbitrary number of
    traits. (The resulting schemas only differ in the
    "length" of the "per_trait" property of
    ``schema["properties"]["SLiM_mutation_list"]["items"]``).

    :param int num_traits: The number of traits in the model.
    :return tskit.MetadataSchema: The metadata schema to be used
        in the node table.
    """
    schema = _raw_slim_metadata_schemas["tree_sequence"]
    schema["struct"]["properties"]["SLiM_mutation_list"]["items"]["properties"][
        "per_trait"
    ]["length"] = num_traits
    return tskit.MetadataSchema(schema)


def slim_individual_metadata_schema(num_traits=1):
    """
    The individual metadata schema depends on the number of traits, and
    {data}`.slim_metadata_schemas`
    returns the schema for a single-trait simulation. This function
    returns the correct schema for a simulation with arbitrary number of
    traits. (The resulting schemas only differ in the
    "length" of the "per_trait" property.)

    :param int num_traits: The number of traits in the model.
    :return tskit.MetadataSchema: The metadata schema to be used
        in the node table.
    """
    schema = _raw_slim_metadata_schemas["individual"]
    schema["properties"]["per_trait"]["length"] = num_traits
    return tskit.MetadataSchema(schema)


def slim_node_metadata_schema(num_chromosomes=1):
    """
    The node metadata schema depends on the number of
    chromosomes in a multichromosome simulation, and
    {data}`.slim_metadata_schemas`
    returns the schema for a single-chromosome simulation. This function
    returns the correct schema for a simulation with arbitrary number of
    chromosomes. (The resulting schemas only differ in the
    ``schema["properties"]["is_vacant"]["length"]`` property,
    which is set to ``floor(num_chromosomes+7)/8``, as described in
    the SLiM manual.)

    :param int num_chromosomes: The number of chromosomes in the model.
    :return tskit.MetadataSchema: The metadata schema to be used
        in the node table.
    """
    # From the SLiM manual on is_vacant:
    # M bytes (uint8_t): a series of bytes comprising a bitfield of is_vacant
    # values, true (1) if this node represents a vacant haplosome for a given
    # chromosome, false (0) otherwise. For chromosomes with indices 0...N−1, the
    # chromosome with index k has its is_vacant bit in bit k%8 of byte k/8, where
    # byte 0 is the first byte in the series of bytes provided, and bit 0 is the
    # least-significant bit, the one with value 0x01 (hexadecimal 1). The number
    # of bytes present, M, is equal to (N+7)/8, the minimum number of bytes
    # necessary. The operators / and % here are integer divide (rounding down)
    # and integer modulo, respectively.
    num_bytes = is_vacant_num_bytes(num_chromosomes)
    schema = _raw_slim_metadata_schemas["node"]
    schema["properties"]["is_vacant"]["length"] = num_bytes
    return tskit.MetadataSchema(schema)


slim_metadata_schemas = {
    k: tskit.MetadataSchema(_raw_slim_metadata_schemas[k])
    for k in _raw_slim_metadata_schemas
    if k != "node"
}
slim_metadata_schemas["node"] = slim_node_metadata_schema()
"""
A dictionary containing the metadata schemas used by SLiM for each of the tables
and for top-level metadata. **Warning:** node metadata schema depends on the
number of chromosomes, and so `slim_metadata_schemas["node"]` is not valid for
a SLiM simulation with more than 8 chromosomes.
"""


def default_slim_metadata(name, num_chromosomes=1, num_traits=1, **kwargs):
    """
    Returns default metadata of type ``name``, where ``name`` is one of
    "tree_sequence", "edge", "site", "mutation", "mutation_list_entry",
    "node", "individual", or "population".

    Additional keyword arguments may be passed in to modify the returned
    metadata dictionary: for each such key=value pair the returned dictionary
    has `value` (re)assigned to `key`.

    :param str name: The type of metadata requested.
    :param int num_chromosomes: The number of chromosomes (only relevant for "node").
    :param int num_traits: The number of traits
        (only relevant for "individual" and "mutation_list_entry").
    :rtype dict:
    """
    if name == "tree_sequence":
        if num_traits == 1:
            trait_names = [
                "simT",
            ]
        else:
            trait_names = [f"trait{j}" for j in range(num_traits)]
        out = {
            "SLiM": {
                "model_type": "nonWF",
                "cycle": 1,
                "tick": 1,
                "file_version": slim_file_version,
                "spatial_dimensionality": "",
                "spatial_periodicity": "",
                "separate_sexes": False,
                "nucleotide_based": False,
                "stage": "late",
                "name": "sim",
                "description": "",
                "this_chromosome": {
                    "id": 1,
                    "index": 0,
                    "symbol": "A",
                    "type": "A",
                },
                "chromosomes": [{"id": 1, "index": 0, "symbol": "A", "type": "A"}],
                "traits": [
                    {"index": j, "name": n, "type": "multiplicative"}
                    for j, n in enumerate(trait_names)
                ],
            },
            "SLiM_mutation_list": [],
        }
    elif name == "edge":
        out = None
    elif name == "site":
        out = None
    elif name == "mutation":
        out = {"derived_states": []}
    elif name == "mutation_list_entry":
        out = {
            "mutation_id": 0,
            "mutation_type": 0,
            "subpopulation": tskit.NULL,
            "slim_time": 0,
            "nucleotide": -1,
            "per_trait": num_traits
            * [{"effect_size": 0.0, "dominance": 0.5, "hemizygous_dominance": 1.0}],
            "padding": None,
        }
    elif name == "node":
        out = {
            "slim_id": tskit.NULL,
            "is_vacant": [0 for _ in range(is_vacant_num_bytes(num_chromosomes))],
        }
    elif name == "individual":
        out = {
            "pedigree_id": tskit.NULL,
            "age": -1,
            "subpopulation": tskit.NULL,
            "sex": -1,
            "pedigree_p1": tskit.NULL,
            "pedigree_p2": tskit.NULL,
            "tag": np.iinfo(np.int64).min,
            "tagF": np.finfo(np.float64).min,
            "tagL0_set": False,
            "tagL0": False,
            "tagL1_set": False,
            "tagL1": False,
            "tagL2_set": False,
            "tagL2": False,
            "tagL3_set": False,
            "tagL3": False,
            "tagL4_set": False,
            "tagL4": False,
            "per_trait": num_traits * [{"phenotype": np.nan, "offset": 1.0}],
        }
    elif name == "population":
        out = {
            "slim_id": tskit.NULL,
            "name": "default",
            "description": "",
            "selfing_fraction": 0.0,
            "female_cloning_fraction": 0.0,
            "male_cloning_fraction": 0.0,
            "sex_ratio": 0.0,
            "bounds_x0": 0.0,
            "bounds_x1": 1.0,
            "bounds_y0": 0.0,
            "bounds_y1": 1.0,
            "bounds_z0": 0.0,
            "bounds_z1": 1.0,
            "migration_records": [],
        }
    else:
        raise ValueError(
            "Unknown metadata request: name should be one of 'tree_sequence', "
            "'edge', 'site', 'mutation', 'mutation_list_entry', 'node', "
            "'individual', or 'population'."
        )
    if out is not None:
        out.update(kwargs)
    return out


###########
# Top-level, a.k.a., tree sequence metadata
###########


def set_tree_sequence_metadata(
    tables,
    model_type,
    tick,
    *,
    cycle=None,
    spatial_dimensionality="",
    spatial_periodicity="",
    separate_sexes=False,
    nucleotide_based=False,
    stage="late",
    name="",
    description="",
    this_chromosome=None,
    chromosomes=None,
    file_version=None,
    set_table_schemas=True,
    traits=None,
    SLiM_mutation_list=None,
):
    if file_version is None:
        file_version = slim_file_version
    if traits is None:
        traits = [{"index": 0, "name": "simT", "type": "multiplicative"}]
    num_traits = len(traits)
    schema_dict = slim_tree_sequence_metadata_schema(num_traits).schema
    old_schema_dict = tables.metadata_schema.schema
    old_json_schema_dict = {}
    old_struct_schema_dict = {}
    tmd = tables.metadata
    if isinstance(tmd, bytes):
        if len(tmd) > 0:
            raise ValueError(
                "Tree sequence has top-level metadata but no schema: this is a problem "
                "since pyslim is trying to add to the metadata."
            )
        metadata_dict = {}
    else:
        # we need to keep other keys in the metadata (and schema) if there are any
        metadata_dict = tables.metadata
        if old_schema_dict["codec"] == "json":
            old_json_schema_dict = tables.metadata_schema.schema
        else:
            assert old_schema_dict["codec"] == "json+struct", (
                "You are using an unexpected codec; "
                "please raise an issue on pyslim if "
                "you need this functionality."
            )
            old_json_schema_dict = tables.metadata_schema.schema["json"]
            old_struct_schema_dict = tables.metadata_schema.schema["struct"]
    if cycle is None:
        cycle = tick
    if chromosomes is None:
        num_chromosomes = 1
    else:
        num_chromosomes = len(chromosomes)
    defaults = default_slim_metadata(
        "tree_sequence", num_chromosomes=num_chromosomes, num_traits=num_traits
    )
    if this_chromosome is None:
        this_chromosome = defaults["SLiM"]["this_chromosome"]
    if chromosomes is None:
        chromosomes = defaults["SLiM"]["chromosomes"]
    if "properties" in old_json_schema_dict:
        schema_dict["json"]["properties"].update(old_json_schema_dict["properties"])
    if "properties" in old_struct_schema_dict:
        schema_dict["struct"]["properties"].update(old_struct_schema_dict["properties"])
    if SLiM_mutation_list is None:
        SLiM_mutation_list = []
    tables.metadata_schema = tskit.MetadataSchema(schema_dict)
    metadata_dict["SLiM"] = {
        "model_type": model_type,
        "tick": tick,
        "cycle": cycle,
        "file_version": file_version,
        "spatial_dimensionality": spatial_dimensionality,
        "spatial_periodicity": spatial_periodicity,
        "separate_sexes": separate_sexes,
        "nucleotide_based": nucleotide_based,
        "stage": stage,
        "name": name,
        "description": description,
        "this_chromosome": this_chromosome,
        "chromosomes": chromosomes,
        "traits": traits,
    }
    metadata_dict["SLiM_mutation_list"] = SLiM_mutation_list
    tables.metadata = metadata_dict
    return metadata_dict


def set_metadata_schemas(tables, num_chromosomes=1, num_traits=1):
    tables.edges.metadata_schema = slim_metadata_schemas["edge"]
    tables.sites.metadata_schema = slim_metadata_schemas["site"]
    tables.mutations.metadata_schema = slim_metadata_schemas["mutation"]
    tables.nodes.metadata_schema = slim_node_metadata_schema(num_chromosomes)
    tables.individuals.metadata_schema = slim_individual_metadata_schema(num_traits)
    tables.populations.metadata_schema = slim_metadata_schemas["population"]


################################
# Previous versions of metadata schema:


def _old_metadata_schema(name, file_version):
    # Returns a metadata schema *if the format has changed*,
    # and None otherwise.
    ms = None
    if name == "tree_sequence" and file_version == "0.9":
        pre_1_0_tree_sequence = {
            "$schema": "http://json-schema.org/schema#",
            "codec": "json",
            "examples": [
                {
                    "SLiM": {
                        "file_version": "0.9",
                        "name": "fox",
                        "description": "foxes on Catalina island",
                        "cycle": 123,
                        "tick": 123,
                        "model_type": "WF",
                        "this_chromosome": {
                            "id": 1,
                            "index": 0,
                            "symbol": "1",
                            "name": "autosome_1",
                            "type": "A",
                        },
                        "chromosomes": [
                            {"id": 1, "symbol": "1", "name": "autosome_1", "type": "A"},
                            {"id": 35, "symbol": "MT", "name": "mtDNA", "type": "HF"},
                        ],
                        "nucleotide_based": False,
                        "separate_sexes": True,
                        "spatial_dimensionality": "xy",
                        "spatial_periodicity": "x",
                    }
                }
            ],
            "properties": {
                "SLiM": {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.9",
                    "properties": {
                        "file_version": {
                            "description": "The SLiM 'file format version' of this tree sequence.",
                            "type": "string",
                        },
                        "name": {
                            "description": "The SLiM species name represented by this tree sequence.",
                            "type": "string",
                        },
                        "description": {
                            "description": "A user-configurable description of the species represented by this tree sequence.",
                            "type": "string",
                        },
                        "cycle": {
                            "description": "The 'SLiM cycle' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "tick": {
                            "description": "The 'SLiM tick' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "model_type": {
                            "description": "The model type used for the last part of this simulation (WF or nonWF).",
                            "enum": ["WF", "nonWF"],
                            "type": "string",
                        },
                        "this_chromosome": {
                            "description": "The chromosome represented by the tree sequence in this file.",
                            "properties": {
                                "id": {
                                    "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
                                    "type": "integer",
                                },
                                "index": {
                                    "description": "The (zero-based) index of this chromosome in the chromosomes metadata array (if present), which should match the information given here.",
                                    "type": "integer",
                                },
                                "symbol": {
                                    "description": 'A short string symbol for the chromosome, unique within this set of tree sequences, such as "1" or "MT".',
                                    "type": "string",
                                },
                                "name": {
                                    "description": "A user-specified name for the chromosome, such as an accession identifier.",
                                    "type": "string",
                                },
                                "type": {
                                    "description": "The type of chromosome, as specified by SLiM.",
                                    "type": "string",
                                },
                            },
                            "required": ["id", "index", "symbol", "type"],
                            "type": "object",
                        },
                        "chromosomes": {
                            "description": "The chromosomes represented by the collection of tree sequences, of which this tree sequence is one member.",
                            "items": {
                                "properties": {
                                    "id": {
                                        "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
                                        "type": "integer",
                                    },
                                    "symbol": {
                                        "description": 'A short string symbol for the chromosome, unique within this set of tree sequences, such as "1" or "MT".',
                                        "type": "string",
                                    },
                                    "name": {
                                        "description": "A user-specified name for the chromosome, such as an accession identifier.",
                                        "type": "string",
                                    },
                                    "type": {
                                        "description": "The type of chromosome, as specified by SLiM.",
                                        "type": "string",
                                    },
                                },
                                "required": ["id", "symbol", "type"],
                                "type": "object",
                            },
                            "type": "array",
                        },
                        "nucleotide_based": {
                            "description": "Whether the simulation was nucleotide-based.",
                            "type": "boolean",
                        },
                        "separate_sexes": {
                            "description": "Whether the simulation had separate sexes.",
                            "type": "boolean",
                        },
                        "spatial_dimensionality": {
                            "description": "The spatial dimensionality of the simulation.",
                            "enum": ["", "x", "xy", "xyz"],
                            "type": "string",
                        },
                        "spatial_periodicity": {
                            "description": "The spatial periodicity of the simulation.",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "type": "string",
                        },
                        "stage": {
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
                            "type": "string",
                        },
                    },
                    "required": [
                        "model_type",
                        "tick",
                        "file_version",
                        "spatial_dimensionality",
                        "spatial_periodicity",
                        "this_chromosome",
                        "separate_sexes",
                        "nucleotide_based",
                    ],
                    "type": "object",
                }
            },
            "required": ["SLiM"],
            "type": "object",
        }

        ms = pre_1_0_tree_sequence

    if name == "tree_sequence" and file_version == "0.8":
        pre_0_9_tree_sequence = {
            "$schema": "http://json-schema.org/schema#",
            "codec": "json",
            "examples": [
                {
                    "SLiM": {
                        "file_version": "0.8",
                        "name": "fox",
                        "description": "foxes on Catalina island",
                        "cycle": 123,
                        "tick": 123,
                        "model_type": "WF",
                        "nucleotide_based": False,
                        "separate_sexes": True,
                        "spatial_dimensionality": "xy",
                        "spatial_periodicity": "x",
                    }
                }
            ],
            "properties": {
                "SLiM": {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.8",
                    "properties": {
                        "file_version": {
                            "description": "The SLiM 'file format version' of this tree sequence.",
                            "type": "string",
                        },
                        "name": {
                            "description": "The SLiM species name represented by this tree sequence.",
                            "type": "string",
                        },
                        "description": {
                            "description": "A user-configurable description of the species represented by this tree sequence.",
                            "type": "string",
                        },
                        "cycle": {
                            "description": "The 'SLiM cycle' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "tick": {
                            "description": "The 'SLiM tick' counter when this tree sequence was recorded.",
                            "type": "integer",
                        },
                        "model_type": {
                            "description": "The model type used for the last part of this simulation (WF or nonWF).",
                            "enum": ["WF", "nonWF"],
                            "type": "string",
                        },
                        "nucleotide_based": {
                            "description": "Whether the simulation was nucleotide-based.",
                            "type": "boolean",
                        },
                        "separate_sexes": {
                            "description": "Whether the simulation had separate sexes.",
                            "type": "boolean",
                        },
                        "spatial_dimensionality": {
                            "description": "The spatial dimensionality of the simulation.",
                            "enum": ["", "x", "xy", "xyz"],
                            "type": "string",
                        },
                        "spatial_periodicity": {
                            "description": "The spatial periodicity of the simulation.",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "type": "string",
                        },
                        "stage": {
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
                            "type": "string",
                        },
                    },
                    "required": [
                        "model_type",
                        "tick",
                        "file_version",
                        "spatial_dimensionality",
                        "spatial_periodicity",
                        "separate_sexes",
                        "nucleotide_based",
                    ],
                    "type": "object",
                }
            },
            "required": ["SLiM"],
            "type": "object",
        }

        ms = pre_0_9_tree_sequence

    if name == "tree_sequence" and file_version in [
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
    ]:
        pre_0_8_tree_sequence = {
            "$schema": "http://json-schema.org/schema#",
            "codec": "json",
            "type": "object",
            "properties": {
                "SLiM": {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.7",
                    "type": "object",
                    "properties": {
                        "model_type": {
                            "type": "string",
                            "enum": ["WF", "nonWF"],
                            "description": "The model type used for the last part of this simulation (WF or nonWF).",
                        },
                        "generation": {
                            "type": "integer",
                            "description": "The 'SLiM generation' counter when this tree sequence was recorded.",
                        },
                        "stage": {
                            "type": "string",
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
                        },
                        "file_version": {
                            "type": "string",
                            "description": "The SLiM 'file format version' of this tree sequence.",
                        },
                        "spatial_dimensionality": {
                            "type": "string",
                            "enum": ["", "x", "xy", "xyz"],
                            "description": "The spatial dimensionality of the simulation.",
                        },
                        "spatial_periodicity": {
                            "type": "string",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "description": "The spatial periodicity of the simulation.",
                        },
                        "separate_sexes": {
                            "type": "boolean",
                            "description": "Whether the simulation had separate sexes.",
                        },
                        "nucleotide_based": {
                            "type": "boolean",
                            "description": "Whether the simulation was nucleotide-based.",
                        },
                    },
                    "required": [
                        "model_type",
                        "generation",
                        "file_version",
                        "spatial_dimensionality",
                        "spatial_periodicity",
                        "separate_sexes",
                        "nucleotide_based",
                    ],
                }
            },
            "required": ["SLiM"],
            "examples": [
                {
                    "SLiM": {
                        "model_type": "WF",
                        "generation": 123,
                        "file_version": "0.7",
                        "spatial_dimensionality": "xy",
                        "spatial_periodicity": "x",
                        "separate_sexes": True,
                        "nucleotide_based": False,
                    }
                }
            ],
        }
        ms = pre_0_8_tree_sequence

    if name == "population" and file_version in [
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
    ]:
        pre_0_7_population = {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for population metadata.",
            "codec": "struct",
            "type": ["object", "null"],
            "properties": {
                "slim_id": {
                    "type": "integer",
                    "description": "The ID of this population in SLiM. Note that this is called a 'subpopulation' in SLiM.",
                    "binaryFormat": "i",
                    "index": 1,
                },
                "selfing_fraction": {
                    "type": "number",
                    "description": "The frequency with which individuals in this subpopulation self (for WF models).",
                    "binaryFormat": "d",
                    "index": 2,
                },
                "female_cloning_fraction": {
                    "type": "number",
                    "description": "The frequency with which females in this subpopulation reproduce clonally (for WF models).",
                    "binaryFormat": "d",
                    "index": 3,
                },
                "male_cloning_fraction": {
                    "type": "number",
                    "description": "The frequency with which males in this subpopulation reproduce clonally (for WF models).",
                    "binaryFormat": "d",
                    "index": 4,
                },
                "sex_ratio": {
                    "type": "number",
                    "description": "This subpopulation's sex ratio (for WF models).",
                    "binaryFormat": "d",
                    "index": 5,
                },
                "bounds_x0": {
                    "type": "number",
                    "description": "The minimum x-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 6,
                },
                "bounds_x1": {
                    "type": "number",
                    "description": "The maximum x-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 7,
                },
                "bounds_y0": {
                    "type": "number",
                    "description": "The minimum y-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 8,
                },
                "bounds_y1": {
                    "type": "number",
                    "description": "The maximum y-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 9,
                },
                "bounds_z0": {
                    "type": "number",
                    "description": "The minimum z-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 10,
                },
                "bounds_z1": {
                    "type": "number",
                    "description": "The maximum z-coordinate in this subpopulation.",
                    "binaryFormat": "d",
                    "index": 11,
                },
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
                                "index": 1,
                            },
                            "migration_rate": {
                                "type": "number",
                                "description": "The fraction of children in this subpopulation that are composed of 'migrants' from the source subpopulation (in WF models).",
                                "binaryFormat": "d",
                                "index": 2,
                            },
                        },
                        "required": ["source_subpop", "migration_rate"],
                        "additionalProperties": False,
                    },
                },
            },
        }
        ms = pre_0_7_population

    if name == "individual" and file_version in ["0.7", "0.8", "0.9"]:
        pre_1_0_individual = {
            "$schema": "http://json-schema.org/schema#",
            "additionalProperties": False,
            "codec": "struct",
            "description": "SLiM schema for individual metadata.",
            "examples": [
                {
                    "age": -1,
                    "flags": 0,
                    "pedigree_id": 123,
                    "pedigree_p1": 12,
                    "pedigree_p2": 23,
                    "sex": 0,
                    "subpopulation": 0,
                }
            ],
            "flags": {
                "SLIM_INDIVIDUAL_METADATA_MIGRATED": {
                    "description": "Whether this individual was a migrant, either in the tick when the tree sequence "
                    "was written out (if the individual was alive then), or in the tick of the last time "
                    "they were Remembered (if not).",
                    "value": 1,
                }
            },
            "properties": {
                "age": {
                    "binaryFormat": "i",
                    "description": "The age of this individual, either when the tree sequence was written out "
                    "(if the individual was alive then), or the last time they were Remembered (if not).",
                    "index": 4,
                    "type": "integer",
                },
                "flags": {
                    "binaryFormat": "I",
                    "description": "Other information about the individual: see 'flags'.",
                    "index": 7,
                    "type": "integer",
                },
                "pedigree_id": {
                    "binaryFormat": "q",
                    "description": "The 'pedigree ID' of this individual in SLiM.",
                    "index": 1,
                    "type": "integer",
                },
                "pedigree_p1": {
                    "binaryFormat": "q",
                    "description": "The 'pedigree ID' of this individual's first parent in SLiM.",
                    "index": 2,
                    "type": "integer",
                },
                "pedigree_p2": {
                    "binaryFormat": "q",
                    "description": "The 'pedigree ID' of this individual's second parent in SLiM.",
                    "index": 3,
                    "type": "integer",
                },
                "sex": {
                    "binaryFormat": "i",
                    "description": "The sex of the individual (0 for female, 1 for male, -1 for hermaphrodite).",
                    "index": 6,
                    "type": "integer",
                },
                "subpopulation": {
                    "binaryFormat": "i",
                    "description": "The ID of the subpopulation the individual was part of, either when the tree sequence "
                    "was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                    "index": 5,
                    "type": "integer",
                },
            },
            "required": [
                "pedigree_id",
                "pedigree_p1",
                "pedigree_p2",
                "age",
                "subpopulation",
                "sex",
                "flags",
            ],
            "type": "object",
        }
        ms = pre_1_0_individual

    if name == "individual" and file_version in [
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
    ]:
        pre_0_7_individual = {
            "$schema": "http://json-schema.org/schema#",
            "description": "SLiM schema for individual metadata.",
            "codec": "struct",
            "type": "object",
            "properties": {
                "pedigree_id": {
                    "type": "integer",
                    "description": "The 'pedigree ID' of this individual in SLiM.",
                    "binaryFormat": "q",
                    "index": 1,
                },
                "age": {
                    "type": "integer",
                    "description": "The age of this individual, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                    "binaryFormat": "i",
                    "index": 2,
                },
                "subpopulation": {
                    "type": "integer",
                    "description": "The ID of the subpopulation the individual was part of, either when the tree sequence was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                    "binaryFormat": "i",
                    "index": 3,
                },
                "sex": {
                    "type": "integer",
                    "description": "The sex of the individual (0 for female, 1 for male, -1 for hermaphrodite).",
                    "binaryFormat": "i",
                    "index": 4,
                },
                "flags": {
                    "type": "integer",
                    "description": "Other information about the individual: see 'flags'.",
                    "binaryFormat": "I",
                    "index": 5,
                },
            },
            "required": ["pedigree_id", "age", "subpopulation", "sex", "flags"],
            "additionalProperties": False,
            "flags": {
                "SLIM_INDIVIDUAL_METADATA_MIGRATED": {
                    "value": 1,
                    "description": "Whether this individual was a migrant, either in the generation when the tree sequence was written out (if the individual was alive then), or in the generation of the last time they were Remembered (if not).",
                }
            },
        }
        ms = pre_0_7_individual

    if name == "mutation" and file_version in [
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
        "0.8",
        "0.9",
    ]:
        mutation_pre_1_0 = {
            "$schema": "http://json-schema.org/schema#",
            "additionalProperties": False,
            "codec": "struct",
            "description": "SLiM schema for mutation metadata.",
            "examples": [
                {
                    "mutation_list": [
                        {
                            "mutation_type": 1,
                            "nucleotide": 3,
                            "selection_coeff": -0.2,
                            "slim_time": 243,
                            "subpopulation": 0,
                        }
                    ]
                }
            ],
            "properties": {
                "mutation_list": {
                    "items": {
                        "additionalProperties": False,
                        "properties": {
                            "mutation_type": {
                                "binaryFormat": "i",
                                "description": "The index of this mutation's mutationType.",
                                "index": 1,
                                "type": "integer",
                            },
                            "nucleotide": {
                                "binaryFormat": "b",
                                "description": "The nucleotide for this mutation (0=A , 1=C , 2=G, 3=T, or -1 for none)",
                                "index": 5,
                                "type": "integer",
                            },
                            "selection_coeff": {
                                "binaryFormat": "f",
                                "description": "This mutation's selection coefficient.",
                                "index": 2,
                                "type": "number",
                            },
                            "slim_time": {
                                "binaryFormat": "i",
                                "description": "The SLiM tick counter when this mutation occurred.",
                                "index": 4,
                                "type": "integer",
                            },
                            "subpopulation": {
                                "binaryFormat": "i",
                                "description": "The ID of the subpopulation this mutation occurred in.",
                                "index": 3,
                                "type": "integer",
                            },
                        },
                        "required": [
                            "mutation_type",
                            "selection_coeff",
                            "subpopulation",
                            "slim_time",
                            "nucleotide",
                        ],
                        "type": "object",
                    },
                    "noLengthEncodingExhaustBuffer": True,
                    "type": "array",
                }
            },
            "required": ["mutation_list"],
            "type": "object",
        }
        ms = mutation_pre_1_0

    if name == "mutation" and file_version in ["0.1", "0.2"]:
        mutation_pre_0_3 = {
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
                                "index": 1,
                            },
                            "selection_coeff": {
                                "type": "number",
                                "description": "This mutation's selection coefficient.",
                                "binaryFormat": "f",
                                "index": 2,
                            },
                            "subpopulation": {
                                "type": "integer",
                                "description": "The ID of the subpopulation this mutation occurred in.",
                                "binaryFormat": "i",
                                "index": 3,
                            },
                            "slim_time": {
                                "type": "integer",
                                "description": "The SLiM generation counter when this mutation occurred.",
                                "binaryFormat": "i",
                                "index": 4,
                            },
                        },
                        "required": [
                            "mutation_type",
                            "selection_coeff",
                            "subpopulation",
                            "slim_time",
                        ],
                        "additionalProperties": False,
                    },
                }
            },
            "required": ["mutation_list"],
            "additionalProperties": False,
        }
        ms = mutation_pre_0_3

    if name == "node" and file_version == "0.9":
        node_0_9 = {
            "$schema": "http://json-schema.org/schema#",
            "additionalProperties": False,
            "codec": "struct",
            "description": "SLiM schema for node metadata.",
            "examples": [{"slim_id": 123, "is_vacant": 0}],
            "properties": {
                "slim_id": {
                    "binaryFormat": "q",
                    "description": "The 'pedigree ID' of the haplosomes associated with this node in SLiM.",
                    "index": 0,
                    "type": "integer",
                },
                "is_vacant": {
                    "description": "A vector of byte (uint8_t) values, with each bit representing whether the node represents a vacant position, either unused or a null haplosome (1), or a non-null haplosome (0), in the corresponding chromosome. This field encodes vacancy for all of the chromosomes in the model, not just the chromosome represented in this file (so that the node table is identical across all chromosomes for a multi-chromosome model). Each chromosome receives one bit here; there are two node table entries per individual, used for the two haplosomes of every chromosome, so only one bit is needed in each entry (making two bits total per chromosome, across the two node table entries). The least significant bit of the first byte is used first (for one haplosome of the first chromosome); the most significant bit of the last byte is used last. The number of bytes present in this field is indicated by this schema's 'binaryFormat' field, which is variable (!), and can also be deduced from the number of chromosomes in the model as given in the top-level 'chromosomes' metadata key, which should always be present if this metadata is present.",
                    "index": 1,
                    "type": "array",
                    "length": 1,  # MAY NEED TO BE CHANGED (in SLiM code is "%d")
                    "items": {"type": "number", "binaryFormat": "B"},
                },
            },
            "required": ["slim_id", "is_vacant"],
            "type": ["object", "null"],
        }
        ms = node_0_9

    if name == "node" and file_version in [
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
        "0.8",
    ]:
        node_pre_0_9 = {
            "$schema": "http://json-schema.org/schema#",
            "additionalProperties": False,
            "codec": "struct",
            "description": "SLiM schema for node metadata.",
            "examples": [{"genome_type": 0, "is_null": False, "slim_id": 123}],
            "properties": {
                "genome_type": {
                    "binaryFormat": "B",
                    "description": "The 'type' of this genome (0 for autosome, 1 for X, 2 for Y).",
                    "index": 2,
                    "type": "integer",
                },
                "is_null": {
                    "binaryFormat": "?",
                    "description": "Whether this node describes a 'null' (non-existant) chromosome.",
                    "index": 1,
                    "type": "boolean",
                },
                "slim_id": {
                    "binaryFormat": "q",
                    "description": "The 'pedigree ID' of this chromosome in SLiM.",
                    "index": 0,
                    "type": "integer",
                },
            },
            "required": ["slim_id", "is_null", "genome_type"],
            "type": ["object", "null"],
        }
        ms = node_pre_0_9

    # everything else's format has remained unchanged
    if ms is not None:
        ms = tskit.MetadataSchema(ms)
    return ms


def _make_mutation_list(mutations, file_version):
    # Prior to 1.0, mutation metadata was a list of entries like:
    # {'mutation_type': 1, 'selection_coeff': -0.1,
    # 'subpopulation': 1, 'slim_time': 5, 'nucleotide': -1}
    # with mutation id stored in the derived state.
    #
    # As of 1.0, this lives in the top-level ts.metadata['SLiM_mutation_list'],
    # with entries like
    # {'mutation_id': 87, 'mutation_type': 1, 'subpopulation': 1, 'slim_time': 5,
    # 'nucleotide': -1, 'padding': None,
    # 'per_trait': [{'effect_size': -0.1, 'dominance': 0.5,
    # 'hemizygous_dominance': 1.0}]}
    if mutations.metadata_schema == tskit.MetadataSchema(None):
        mutations.metadata_schema = _old_metadata_schema("mutation", file_version)
    mutation_list = []
    for mut in mutations:
        for sid, md in zip(mut.derived_state.split(","), mut.metadata["mutation_list"]):
            if "nucleotide" not in md:
                md["nucleotide"] = -1
            md["mutation_id"] = int(sid)
            md["per_trait"] = [
                {
                    "effect_size": md["selection_coeff"],
                    "dominance": 0.5,  # WE DON'T KNOW THIS
                    "hemizygous_dominance": 1.0,  # OR THIS
                },
            ]
            del md["selection_coeff"]
            md["padding"] = None
            mutation_list.append(md)
    return mutation_list


def is_current_version(ts, _warn=False):
    """
    Tests whether the metadata provided is the current SLiM file format or not.
    If not, use `pyslim.update( )` to bring it up to date.

    This method may be provided either a TreeSequence or TableCollection directly,
    or the metadata from one of these. The latter is useful because
    accessing top-level metadata can be a costly operation.

    :param ts: Either the top-level metadata of a tree sequence,
        or a TreeSequence or TableCollection that carries this metadata.
    :return bool: Whether the tree sequence is the current version.
    """
    if (
        isinstance(ts, tskit.TreeSequence)
        or isinstance(ts, tskit.TableCollection)
        or isinstance(ts, tskit.ImmutableTableCollection)
    ):
        ts = ts.metadata
    out = (
        isinstance(ts, dict)
        and ("SLiM" in ts)
        and (ts["SLiM"]["file_version"] == slim_file_version)
    )
    if _warn and not out:
        warnings.warn(
            "This tree sequence is not the current SLiM format. "
            "Use `pyslim.update( )` to update the tree sequence."
        )
    return out


def update(ts):
    """
    Update a tree sequence produced by a previous version of SLiM
    to the current file version.

    :return TreeSequence: The updated tree sequence.
    """
    tables = ts.dump_tables()
    update_tables(tables)
    return tables.tree_sequence()


def update_tables(tables):
    """
    Update tables produced by a previous version of SLiM to the current file version.
    Modifies the tables in place.
    """
    # First we ensure we can find the file format version number
    # in top-level metadata. Then we proceed to fix up the tables as necessary.
    md = tables.metadata
    if not (isinstance(md, dict) and "SLiM" in md):
        # Old versions kept information in provenance, not top-level metadata.
        # Note this uses defaults on keys not present in provenance,
        # which prior to 0.5 was everything but generation and model_type.
        # Recovering from provenance has also been useful for operations
        # that discard metadata (eg as msprime did prior to 0.7.5).
        values = default_slim_metadata("tree_sequence")["SLiM"]
        prov = None
        file_version = "unknown"
        # use only the last SLiM provenance
        for p in tables.provenances:
            is_slim, this_file_version = slim_provenance_version(p)
            if is_slim:
                prov = p
                file_version = this_file_version
        values["file_version"] = file_version
        try:
            record = json.loads(prov.record)
            if file_version == "0.1":
                values["model_type"] = record["model_type"]
                values["tick"] = record["generation"]
                values["cycle"] = record["generation"]
            else:
                if "generation" in record["slim"]:
                    values["tick"] = record["slim"]["generation"]
                    values["cycle"] = record["slim"]["generation"]
                for k in values:
                    if k in record["parameters"]:
                        values[k] = record["parameters"][k]
                    if k in record["slim"]:
                        values[k] = record["slim"][k]
        except:
            raise ValueError("Failed to obtain metadata from provenance.")
        md = set_tree_sequence_metadata(tables, **values)

    file_version = md["SLiM"]["file_version"]
    if file_version != slim_file_version:
        warnings.warn(
            f"This is a version {file_version} SLiM tree sequence. "
            "If you write this out to a file, "
            f"it will be converted to version {slim_file_version}."
        )

        old_schema = _old_metadata_schema("tree_sequence", file_version)
        if old_schema is not None:
            assert (
                "struct" not in old_schema.schema
                or "SLiM_mutation_list" not in old_schema.schema["struct"]["properties"]
            )
            # we should get the number of traits from the metadata,
            # but for old file versions, this won't be present
            assert (
                "json" not in old_schema.schema
                or "traits" not in old_schema.schema["json"]["properties"]
            )
            md["SLiM_mutation_list"] = _make_mutation_list(
                tables.mutations, file_version
            )
            num_traits = 1
            new_schema = slim_tree_sequence_metadata_schema(num_traits=num_traits)
            new_properties = new_schema.asdict()["json"]["properties"]["SLiM"][
                "required"
            ]
            tables.metadata_schema = new_schema
            defaults = default_slim_metadata("tree_sequence", num_traits=num_traits)
            for k in new_properties:
                if k not in md["SLiM"]:
                    if k == "tick":
                        md["SLiM"]["tick"] = md["SLiM"]["generation"]
                        md["SLiM"]["cycle"] = md["SLiM"]["generation"]
                    else:
                        md["SLiM"][k] = defaults["SLiM"][k]
            tables.metadata = md

        old_schema = _old_metadata_schema("node", file_version)
        if old_schema is not None:
            nodes = tables.nodes.copy()
            tables.nodes.clear()
            if nodes.metadata_schema == tskit.MetadataSchema(None):
                nodes.metadata_schema = old_schema
            if "chromosomes" not in md["SLiM"]:
                num_chroms = 1
            else:
                num_chroms = len(md["SLiM"]["chromosomes"])
            new_schema = slim_node_metadata_schema(num_chroms)
            tables.nodes.metadata_schema = new_schema
            new_node_schema = new_schema
            if file_version in ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"):
                # 0.8->0.9 switched from is_null to is_vacant,
                # and moved chromosome type from node metadata to top-level
                not_vacant = [0]  # single chromosome
                yes_vacant = [1]
                gt = None
                for n in nodes:
                    md = n.metadata
                    if len(md) > 0:
                        md["is_vacant"] = yes_vacant if md["is_null"] else not_vacant
                        if not md["is_null"]:
                            if gt is None:
                                gt = md["genome_type"]
                            else:
                                assert md["genome_type"] == gt, (
                                    "Inconsistent tables: "
                                    f"mismatching genome types {gt} and "
                                    f"{md['genome_type']} in node metadata."
                                )
                        del md["is_null"]
                        del md["genome_type"]
                    tables.nodes.append(n.replace(metadata=md))
                # flags for node genome type pre-0.9:
                # confusingly and sub-optimally, these were redundant:
                # all non-null nodes in the same sim would have the same genome type
                GENOME_TYPE_AUTOSOME = 0
                GENOME_TYPE_X = 1
                GENOME_TYPE_Y = 2
                top_md = tables.metadata
                i = top_md["SLiM"]["this_chromosome"]["index"]
                # 'chromosomes' is not required so won't be inserted,
                # so we won't update it also
                assert "chromosomes" not in top_md["SLiM"]["this_chromosome"], (
                    ""
                    "This is an unexpected result: if you hit this, "
                    "please file a bug at "
                    "https://github.com/tskit-dev/pyslim."
                )
                if gt == GENOME_TYPE_X:
                    top_md["SLiM"]["this_chromosome"]["type"] = "X"
                    top_md["SLiM"]["this_chromosome"]["symbol"] = "X"
                    # if "chromosomes" in tables.metadata['SLiM']:
                    #     top_md['SLiM']['chromosomes'][i]['type'] = "X"
                    #     top_md['SLiM']['chromosomes'][i]['symbol'] = "X"
                elif gt == GENOME_TYPE_Y:
                    top_md["SLiM"]["this_chromosome"]["type"] = "-Y"
                    top_md["SLiM"]["this_chromosome"]["symbol"] = "Y"
                    #     top_md['SLiM']['chromosomes'][i]['type'] = "Y"
                    #     top_md['SLiM']['chromosomes'][i]['symbol'] = "Y"
                else:
                    assert gt == GENOME_TYPE_AUTOSOME
                tables.metadata = top_md
            else:
                assert file_version == "0.9"
                # just needs recoding (and doesn't really need that, we just
                # changed the index of some entries in the schema)
                for n in nodes:
                    tables.nodes.append(n)

        old_schema = _old_metadata_schema("population", file_version)
        if old_schema is not None:
            pops = tables.populations.copy()
            tables.populations.clear()
            if pops.metadata_schema == tskit.MetadataSchema(None):
                pops.metadata_schema = old_schema
            new_schema = slim_metadata_schemas["population"]
            tables.populations.metadata_schema = new_schema
            # just needs recoding
            for pop in pops:
                tables.populations.append(pop)

        old_schema = _old_metadata_schema("individual", file_version)
        if old_schema is not None:
            inds = tables.individuals.copy()
            tables.individuals.clear()
            if inds.metadata_schema == tskit.MetadataSchema(None):
                inds.metadata_schema = old_schema
            num_traits = len(tables.metadata["SLiM"]["traits"])
            new_schema = slim_individual_metadata_schema(num_traits=num_traits)
            tables.individuals.metadata_schema = new_schema
            # new(er) additions are pedigree_pX in 0.7
            # and per_trait in 1.0
            # also flags got removed in 1.0
            defaults = default_slim_metadata("individual", num_traits=num_traits)
            for ind in inds:
                md = ind.metadata
                flags = ind.flags
                for k in [
                    "pedigree_p1",
                    "pedigree_p2",
                    "tag",
                    "tagF",
                    "tagL0",
                    "tagL0_set",
                    "tagL1",
                    "tagL1_set",
                    "tagL2",
                    "tagL2_set",
                    "tagL3",
                    "tagL3_set",
                    "tagL4",
                    "tagL4_set",
                ]:
                    md.setdefault(k, defaults[k])
                assert "flags" in md
                # MIGRATED got moved over in 1.0, and previous flag was =1
                flags |= (
                    pyslim.INDIVIDUAL_MIGRATED
                    if (md["flags"] & _INDIVIDUAL_FLAG_MIGRATED_DEPRECATED > 0)
                    else 0
                )
                del md["flags"]
                if "per_trait" not in md:
                    md["per_trait"] = copy.deepcopy(defaults["per_trait"])
                tables.individuals.append(ind.replace(metadata=md, flags=flags))

        old_schema = _old_metadata_schema("mutation", file_version)
        if old_schema is not None:
            muts = tables.mutations.copy()
            tables.mutations.clear()
            if muts.metadata_schema == tskit.MetadataSchema(None):
                muts.metadata_schema = old_schema
            tables.mutations.metadata_schema = slim_metadata_schemas["mutation"]
            for mut in muts:
                md = {"derived_states": [int(j) for j in mut.derived_state.split(",")]}
                tables.mutations.append(mut.replace(metadata=md))

        if file_version == "0.1":
            # shift times
            slim_generation = tables.metadata["SLiM"]["tick"]
            node_times = tables.nodes.time + slim_generation
            tables.nodes.set_columns(
                flags=tables.nodes.flags,
                time=node_times,
                population=tables.nodes.population,
                individual=tables.nodes.individual,
                metadata=tables.nodes.metadata,
                metadata_offset=tables.nodes.metadata_offset,
            )
            migration_times = tables.migrations.time + slim_generation
            tables.migrations.set_columns(
                left=tables.migrations.left,
                right=tables.migrations.right,
                node=tables.migrations.node,
                source=tables.migrations.source,
                dest=tables.migrations.dest,
                time=migration_times,
            )

        new_record = {
            "schema_version": "1.0.0",
            "software": {
                "name": "pyslim",
                "version": pyslim_version,
            },
            "parameters": {
                "command": ["updrade_tables"],
                "old_file_version": file_version,
                "new_file_version": slim_file_version,
            },
            "environment": get_environment(),
        }
        tskit.validate_provenance(new_record)
        tables.provenances.add_row(json.dumps(new_record))

        md = tables.metadata
        md["SLiM"]["file_version"] = slim_file_version
        tables.metadata = md
