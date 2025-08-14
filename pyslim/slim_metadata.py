import tskit
import json
import warnings
import numpy as np

from ._version import *
from .provenance import slim_provenance_version, get_environment

from pyslim import INDIVIDUAL_ALIVE
from pyslim import GENOME_TYPE_AUTOSOME
from pyslim import GENOME_TYPE_X
from pyslim import GENOME_TYPE_Y

def is_vacant_num_bytes(num_chromosomes):
    """
    TODO document
    """
    # see _is_chrom_vacant
    return int((num_chromosomes + 7)/8)


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
    "tree_sequence" : {
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
             "type": "A"
            },
            "chromosomes": [
             {
              "id": 1,
              "symbol": "1",
              "name": "autosome_1",
              "type": "A"
             },
             {
              "id": 35,
              "symbol": "MT",
              "name": "mtDNA",
              "type": "HF"
             }
            ],
            "nucleotide_based": False,
            "separate_sexes": True,
            "spatial_dimensionality": "xy",
            "spatial_periodicity": "x"
           }
          }
         ],
         "properties": {
          "SLiM": {
           "description": "Top-level metadata for a SLiM tree sequence, file format version 0.9",
           "properties": {
            "file_version": {
             "description": "The SLiM 'file format version' of this tree sequence.",
             "type": "string"
            },
            "name": {
             "description": "The SLiM species name represented by this tree sequence.",
             "type": "string"
            },
            "description": {
             "description": "A user-configurable description of the species represented by this tree sequence.",
             "type": "string"
            },
            "cycle": {
             "description": "The 'SLiM cycle' counter when this tree sequence was recorded.",
             "type": "integer"
            },
            "tick": {
             "description": "The 'SLiM tick' counter when this tree sequence was recorded.",
             "type": "integer"
            },
            "model_type": {
             "description": "The model type used for the last part of this simulation (WF or nonWF).",
             "enum": [
              "WF",
              "nonWF"
             ],
             "type": "string"
            },
            "this_chromosome": {
             "description": "The chromosome represented by the tree sequence in this file.",
             "properties": {
              "id": {
               "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
               "type": "integer"
              },
              "index": {
               "description": "The (zero-based) index of this chromosome in the chromosomes metadata array (if present), which should match the information given here.",
               "type": "integer"
              },
              "symbol": {
               "description": "A short string symbol for the chromosome, unique within this set of tree sequences, such as \"1\" or \"MT\".",
               "type": "string"
              },
              "name": {
               "description": "A user-specified name for the chromosome, such as an accession identifier.",
               "type": "string"
              },
              "type": {
               "description": "The type of chromosome, as specified by SLiM.",
               "type": "string"
              }
             },
             "required": [
              "id",
              "index",
              "symbol",
              "type"
             ],
             "type": "object"
            },
            "chromosomes": {
             "description": "The chromosomes represented by the collection of tree sequences, of which this tree sequence is one member.",
             "items": {
              "properties": {
               "id": {
                "description": "An integer identifier for the chromosome, unique within this set of tree sequences; often the chromosome number in the organism being represented, such as 1.",
                "type": "integer"
               },
               "symbol": {
                "description": "A short string symbol for the chromosome, unique within this set of tree sequences, such as \"1\" or \"MT\".",
                "type": "string"
               },
               "name": {
                "description": "A user-specified name for the chromosome, such as an accession identifier.",
                "type": "string"
               },
               "type": {
                "description": "The type of chromosome, as specified by SLiM.",
                "type": "string"
               }
              },
              "required": [
               "id",
               "symbol",
               "type"
              ],
              "type": "object"
             },
             "type": "array"
            },
            "nucleotide_based": {
             "description": "Whether the simulation was nucleotide-based.",
             "type": "boolean"
            },
            "separate_sexes": {
             "description": "Whether the simulation had separate sexes.",
             "type": "boolean"
            },
            "spatial_dimensionality": {
             "description": "The spatial dimensionality of the simulation.",
             "enum": [
              "",
              "x",
              "xy",
              "xyz"
             ],
             "type": "string"
            },
            "spatial_periodicity": {
             "description": "The spatial periodicity of the simulation.",
             "enum": [
              "",
              "x",
              "y",
              "z",
              "xy",
              "xz",
              "yz",
              "xyz"
             ],
             "type": "string"
            },
            "stage": {
             "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
             "type": "string"
            }
           },
           "required": [
            "model_type",
            "tick",
            "file_version",
            "spatial_dimensionality",
            "spatial_periodicity",
            "this_chromosome",
            "separate_sexes",
            "nucleotide_based"
           ],
           "type": "object"
          }
         },
         "required": [
          "SLiM"
         ],
         "type": "object"
    }
    ,
    "edge" : None,
    "site" : None,
    "mutation" : 
    {
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
                        "subpopulation": 0
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
                            "type": "integer"
                        },
                        "nucleotide": {
                            "binaryFormat": "b",
                            "description": "The nucleotide for this mutation (0=A , 1=C , 2=G, 3=T, or -1 for none)",
                            "index": 5,
                            "type": "integer"
                        },
                        "selection_coeff": {
                            "binaryFormat": "f",
                            "description": "This mutation's selection coefficient.",
                            "index": 2,
                            "type": "number"
                        },
                        "slim_time": {
                            "binaryFormat": "i",
                            "description": "The SLiM tick counter when this mutation occurred.",
                            "index": 4,
                            "type": "integer"
                        },
                        "subpopulation": {
                            "binaryFormat": "i",
                            "description": "The ID of the subpopulation this mutation occurred in.",
                            "index": 3,
                            "type": "integer"
                        }
                    },
                    "required": [
                        "mutation_type",
                        "selection_coeff",
                        "subpopulation",
                        "slim_time",
                        "nucleotide"
                    ],
                    "type": "object"
                },
                "noLengthEncodingExhaustBuffer": True,
                "type": "array"
            }
        },
        "required": [
                "mutation_list"
        ],
        "type": "object"
    },
    "node" : 
    {
     "$schema": "http://json-schema.org/schema#",
     "additionalProperties": False,
     "codec": "struct",
     "description": "SLiM schema for node metadata.",
     "examples": [
      {
       "slim_id": 123,
       "is_vacant": 0
      }
     ],
     "properties": {
      "slim_id": {
       "binaryFormat": "q",
       "description": "The 'pedigree ID' of the haplosomes associated with this node in SLiM.",
       "index": 0,
       "type": "integer"
      },
      "is_vacant": {
       "description": "A vector of byte (uint8_t) values, with each bit representing whether the node represents a vacant position, either unused or a null haplosome (1), or a non-null haplosome (0), in the corresponding chromosome. This field encodes vacancy for all of the chromosomes in the model, not just the chromosome represented in this file (so that the node table is identical across all chromosomes for a multi-chromosome model). Each chromosome receives one bit here; there are two node table entries per individual, used for the two haplosomes of every chromosome, so only one bit is needed in each entry (making two bits total per chromosome, across the two node table entries). The least significant bit of the first byte is used first (for one haplosome of the first chromosome); the most significant bit of the last byte is used last. The number of bytes present in this field is indicated by this schema's 'binaryFormat' field, which is variable (!), and can also be deduced from the number of chromosomes in the model as given in the top-level 'chromosomes' metadata key, which should always be present if this metadata is present.",
       "index": 1,
       "type": "array",
       "length": 1,  # MAY NEED TO BE CHANGED (in SLiM code is "%d")
       "items": {
        "type": "number",
        "binaryFormat": "B"
       }
      }
     },
     "required": [
      "slim_id",
      "is_vacant"
     ],
     "type": [
      "object",
      "null"
     ]
    },
    "individual" :
    {
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
                "subpopulation": 0
            }
        ],
        "flags": {
            "SLIM_INDIVIDUAL_METADATA_MIGRATED": {
                "description": "Whether this individual was a migrant, either in the tick when the tree sequence "
                               "was written out (if the individual was alive then), or in the tick of the last time "
                               "they were Remembered (if not).",
                "value": 1
            }
        },
        "properties": {
            "age": {
                "binaryFormat": "i",
                "description": "The age of this individual, either when the tree sequence was written out "
                               "(if the individual was alive then), or the last time they were Remembered (if not).",
                               "index": 4,
                               "type": "integer"
            },
            "flags": {
                "binaryFormat": "I",
                "description": "Other information about the individual: see 'flags'.",
                "index": 7,
                "type": "integer"
            },
            "pedigree_id": {
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual in SLiM.",
                "index": 1,
                "type": "integer"
            },
            "pedigree_p1": {
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual's first parent in SLiM.",
                "index": 2,
                "type": "integer"
            },
            "pedigree_p2": {
                "binaryFormat": "q",
                "description": "The 'pedigree ID' of this individual's second parent in SLiM.",
                "index": 3,
                "type": "integer"
            },
            "sex": {
                "binaryFormat": "i",
                "description": "The sex of the individual (0 for female, 1 for male, -1 for hermaphrodite).",
                "index": 6,
                "type": "integer"
            },
            "subpopulation": {
                "binaryFormat": "i",
                "description": "The ID of the subpopulation the individual was part of, either when the tree sequence "
                               "was written out (if the individual was alive then), or the last time they were Remembered (if not).",
                               "index": 5,
                               "type": "integer"
            }
        },
        "required": [
                "pedigree_id",
                "pedigree_p1",
                "pedigree_p2",
                "age",
                "subpopulation",
                "sex",
                "flags"
            ],
        "type": "object"
    },
    "population": 
    {
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
                        {
                            "migration_rate": 0.9,
                            "source_subpop": 1
                        },
                        {
                            "migration_rate": 0.1,
                            "source_subpop": 2
                        }
                    ],
                    "selfing_fraction": 0.5,
                    "sex_ratio": 0.5,
                    "slim_id": 2,
                    "name": "p2"
                }
            ],
            "properties": {
                "bounds_x0": {
                    "description": "The minimum x-coordinate in this subpopulation.",
                    "type": "number"
                },
                "bounds_x1": {
                    "description": "The maximum x-coordinate in this subpopulation.",
                    "type": "number"
                },
                "bounds_y0": {
                    "description": "The minimum y-coordinate in this subpopulation.",
                    "type": "number"
                },
                "bounds_y1": {
                    "description": "The maximum y-coordinate in this subpopulation.",
                    "type": "number"
                },
                "bounds_z0": {
                    "description": "The minimum z-coordinate in this subpopulation.",
                    "type": "number"
                },
                "bounds_z1": {
                    "description": "The maximum z-coordinate in this subpopulation.",
                    "type": "number"
                },
                "description": {
                    "description": "A description of this subpopulation.",
                    "type": "string"
                },
                "female_cloning_fraction": {
                    "description": "The frequency with which females in this subpopulation reproduce clonally (for WF models).",
                    "type": "number"
                },
                "male_cloning_fraction": {
                    "description": "The frequency with which males in this subpopulation reproduce clonally (for WF models).",
                    "type": "number"
                },
                "migration_records": {
                    "items": {
                        "properties": {
                            "migration_rate": {
                                "description": "The fraction of children in this subpopulation that are composed of 'migrants' from the source subpopulation (in WF models).",
                                "type": "number"
                            },
                            "source_subpop": {
                                "description": "The ID of the subpopulation migrants come from (in WF models).",
                                "type": "integer"
                            }
                        },
                        "required": [
                            "source_subpop",
                            "migration_rate"
                        ],
                        "type": "object"
                    },
                    "type": "array"
                },
                "name": {
                    "description": "A human-readable name for this subpopulation.",
                    "type": "string"
                },
                "selfing_fraction": {
                    "description": "The frequency with which individuals in this subpopulation self (for WF models).",
                    "type": "number"
                },
                "sex_ratio": {
                    "description": "This subpopulation's sex ratio (for WF models).",
                    "type": "number"
                },
                "slim_id": {
                        "description": "The ID of this population in SLiM. Note that this is called a 'subpopulation' in SLiM.",
                        "type": "integer"
                    }
        },
        "required": [],
        "type": [
                "object",
                "null"
        ]
    }
}


def slim_node_metadata_schema(num_chromosomes=1):
    """
    Unlike other schema, the node metadata schema depends on the number of
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


def default_slim_metadata(name, num_chromosomes=1):
    """
    Returns default metadata of type ``name``, where ``name`` is one of
    "tree_sequence", "edge", "site", "mutation", "mutation_list_entry",
    "node", "individual", or "population".

    :param str name: The type of metadata requested.
    :rtype dict:
    """
    if name == "tree_sequence":
        out = {
             "SLiM" : {
                 "model_type" : "nonWF",
                 "cycle" : 1,
                 "tick" : 1,
                 "file_version" : slim_file_version,
                 "spatial_dimensionality" : "",
                 "spatial_periodicity" : "",
                 "separate_sexes" : False,
                 "nucleotide_based" : False,
                 "stage" : "late",
                 "name" : "sim",
                 "description" : "",
                 "this_chromosome" : {
                     "id": 1,
                     "index": 0,
                     "symbol": "A",
                     "type": "A",
                 },
                 "chromosomes": [{
                     'id': 1,
                     "index": 0,
                     'symbol': 'A',
                     'type': 'A'
                 }],
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
    elif name == "mutation_list_entry":
        out = {
            "mutation_type": 0,
            "selection_coeff": 0.0,
            "subpopulation": tskit.NULL,
            "slim_time": 0,
            "nucleotide": -1,
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
            "flags": 0,
            "pedigree_p1": tskit.NULL,
            "pedigree_p2": tskit.NULL,
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
            "migration_records": []
        }
    else:
        raise ValueError(
            "Unknown metadata request: name should be one of 'tree_sequence', "
            "'edge', 'site', 'mutation', 'mutation_list_entry', 'node', "
            "'individual', or 'population'.")
    return out


###########
# Top-level, a.k.a., tree sequence metadata
###########

def set_tree_sequence_metadata(tables,
        model_type,
        tick,
        cycle=None,
        spatial_dimensionality='',
        spatial_periodicity='',
        separate_sexes=False,
        nucleotide_based=False,
        stage='late',
        name='',
        description='',
        this_chromosome=None,
        chromosomes=None,
        file_version=None,
        set_table_schemas=True):
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
    if cycle is None:
        cycle = tick
    defaults = default_slim_metadata("tree_sequence")
    if this_chromosome is None:
        this_chromosome = defaults['SLiM']['this_chromosome']
    if chromosomes is None:
        chromosomes = defaults['SLiM']['chromosomes']
    assert(schema_dict['codec'] == 'json')
    assert(schema_dict['type'] == 'object')
    if "properties" not in schema_dict:
        schema_dict["properties"] = {}
    schema_dict['properties']['SLiM'] = slim_metadata_schemas['tree_sequence'].schema['properties']['SLiM']
    tables.metadata_schema = tskit.MetadataSchema(schema_dict)
    metadata_dict['SLiM'] = {
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
            }
    tables.metadata = metadata_dict


def set_metadata_schemas(tables, num_chromosomes=1):
    tables.edges.metadata_schema = slim_metadata_schemas['edge']
    tables.sites.metadata_schema = slim_metadata_schemas['site']
    tables.mutations.metadata_schema = slim_metadata_schemas['mutation']
    tables.nodes.metadata_schema = slim_node_metadata_schema(num_chromosomes)
    tables.individuals.metadata_schema = slim_metadata_schemas['individual']
    tables.populations.metadata_schema = slim_metadata_schemas['population']


################################
# Previous versions of metadata schema:

def _old_metadata_schema(name, file_version):
    # Returns a metadata schema *if the format has changed*,
    # and None otherwise.
    ms = None
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
                        "spatial_periodicity": "x"
                    }
                }
            ],
            "properties": {
                "SLiM": {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.8",
                    "properties": {
                        "file_version": {
                            "description": "The SLiM 'file format version' of this tree sequence.",
                            "type": "string"
                        },
                        "name": {
                            "description": "The SLiM species name represented by this tree sequence.",
                            "type": "string"
                        },
                        "description": {
                            "description": "A user-configurable description of the species represented by this tree sequence.",
                            "type": "string"
                        },
                        "cycle": {
                            "description": "The 'SLiM cycle' counter when this tree sequence was recorded.",
                            "type": "integer"
                        },
                        "tick": {
                            "description": "The 'SLiM tick' counter when this tree sequence was recorded.",
                            "type": "integer"
                        },
                        "model_type": {
                            "description": "The model type used for the last part of this simulation (WF or nonWF).",
                            "enum": [ "WF", "nonWF" ],
                            "type": "string"
                        },
                        "nucleotide_based": {
                            "description": "Whether the simulation was nucleotide-based.",
                            "type": "boolean"
                        },
                        "separate_sexes": {
                            "description": "Whether the simulation had separate sexes.",
                            "type": "boolean"
                        },
                        "spatial_dimensionality": {
                            "description": "The spatial dimensionality of the simulation.",
                            "enum": ["", "x", "xy", "xyz"],
                            "type": "string"
                        },
                        "spatial_periodicity": {
                            "description": "The spatial periodicity of the simulation.",
                            "enum": ["", "x", "y", "z", "xy", "xz", "yz", "xyz"],
                            "type": "string"
                        },
                        "stage": {
                            "description": "The stage of the SLiM life cycle when this tree sequence was recorded.",
                            "type": "string"
                        }
                    },
                    "required": [
                            "model_type",
                            "tick",
                            "file_version",
                            "spatial_dimensionality",
                            "spatial_periodicity",
                            "separate_sexes",
                            "nucleotide_based"
                    ],
                    "type": "object"
                }
            },
            "required": ["SLiM"],
            "type": "object"
        }

        ms = pre_0_9_tree_sequence

    if (name == "tree_sequence"
        and file_version in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7"]):
        pre_0_8_tree_sequence = {
            "$schema": "http://json-schema.org/schema#",
            "codec": "json",
            "type": "object",
            "properties": {
                "SLiM" : {
                    "description": "Top-level metadata for a SLiM tree sequence, file format version 0.7",
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
                     "file_version" : "0.7",
                     "spatial_dimensionality" : "xy",
                     "spatial_periodicity" : "x",
                     "separate_sexes" : True,
                     "nucleotide_based" : False
                     }
                }
            ]
        }
        ms = pre_0_8_tree_sequence 

    if (name == "population"
        and file_version in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6"]):
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
            }
        }
        ms = pre_0_7_population 

    if (name == "individual"
        and file_version in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6"]):
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
        }
        ms = pre_0_7_individual 

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
                        },
                        "required": ["mutation_type", "selection_coeff", "subpopulation",
                                     "slim_time"],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["mutation_list"],
            "additionalProperties": False,
        }
        ms = mutation_pre_0_3

    if (name == "node"
        and file_version in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"]):
        node_pre_0_9 = {
           "$schema": "http://json-schema.org/schema#",
           "additionalProperties": False,
           "codec": "struct",
           "description": "SLiM schema for node metadata.",
           "examples": [
               {
                   "genome_type": 0,
                   "is_null": False,
                   "slim_id": 123
                }
            ],
           "properties": {
               "genome_type": {
                   "binaryFormat": "B",
                   "description": "The 'type' of this genome (0 for autosome, 1 for X, 2 for Y).",
                   "index": 2,
                   "type": "integer"
                },
               "is_null": {
                   "binaryFormat": "?",
                   "description": "Whether this node describes a 'null' (non-existant) chromosome.",
                   "index": 1,
                   "type": "boolean"
                },
               "slim_id": {
                   "binaryFormat": "q",
                   "description": "The 'pedigree ID' of this chromosome in SLiM.",
                   "index": 0,
                   "type": "integer"
                }
            },
           "required": [
               "slim_id",
               "is_null",
               "genome_type"
            ],
           "type": [
               "object",
               "null"
            ]
        }
        ms = node_pre_0_9

    # everything else's format has remained unchanged
    if ms is not None:
        ms = tskit.MetadataSchema(ms)
    return ms


def is_current_version(ts, _warn=False):
    """
    Tests whether the tree sequence or table collection provided is the current
    SLiM file format or not. If not, use `pyslim.update( )` to bring it up to
    date.

    :param TreeSequence ts: The tree sequence or table collection.
    :return bool: Whether the tree sequence is the current version.
    """
    out = (
        isinstance(ts.metadata, dict)
        and
        ('SLiM' in ts.metadata)
        and
        (ts.metadata['SLiM']['file_version'] == slim_file_version)
    )
    if _warn and not out:
        warnings.warn(
                "This tree sequence is not the current SLiM format, "
                "so some operations may not work. "
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
    if not (isinstance(tables.metadata, dict) and 'SLiM' in tables.metadata):
        # Old versions kept information in provenance, not top-level metadata.
        # Note this uses defaults on keys not present in provenance,
        # which prior to 0.5 was everything but generation and model_type.
        values = default_slim_metadata('tree_sequence')['SLiM']
        prov = None
        file_version = 'unknown'
        # use only the last SLiM provenance
        for p in tables.provenances:
            is_slim, this_file_version = slim_provenance_version(p) 
            if is_slim:
                prov = p
                file_version = this_file_version
        values['file_version'] = file_version
        try:
            record = json.loads(prov.record)
            if file_version == "0.1":
                values['model_type'] = record['model_type']
                values['tick'] = record['generation']
                values['cycle'] = record['generation']
            else:
                if 'generation' in record['slim']:
                    values['tick'] = record['slim']['generation']
                    values['cycle'] = record['slim']['generation']
                for k in values:
                    if k in record['parameters']:
                        values[k] = record['parameters'][k]
                    if k in record['slim']:
                        values[k] = record['slim'][k]
        except:
            raise ValueError("Failed to obtain metadata from provenance.")
        set_tree_sequence_metadata(tables, **values)

    file_version = tables.metadata['SLiM']['file_version']
    if file_version != slim_file_version:
        warnings.warn("This is a version {} SLiM tree sequence.".format(file_version) +
                      " If you write this out to a file, " +
                      "it will be converted to version {}.".format(slim_file_version))

        # the only tables to have metadata schema changed thus far
        # are nodes, populations, individuals, mutations, and top-level:
        old_schema = _old_metadata_schema("tree_sequence", file_version)
        if old_schema is not None:
            md = tables.metadata
            new_schema = slim_metadata_schemas["tree_sequence"]
            new_properties = new_schema.asdict()['properties']['SLiM']['required']
            tables.metadata_schema = new_schema
            defaults = default_slim_metadata("tree_sequence")
            for k in new_properties:
                if k not in md['SLiM']:
                    if k == "tick":
                        md['SLiM']['tick'] = md['SLiM']['generation']
                        md['SLiM']['cycle'] = md['SLiM']['generation']
                    else:
                        md['SLiM'][k] = defaults['SLiM'][k]
            tables.metadata = md

        old_schema = _old_metadata_schema("node", file_version)
        if old_schema is not None:
            nodes = tables.nodes.copy()
            tables.nodes.clear()
            if nodes.metadata_schema == tskit.MetadataSchema(None):
                nodes.metadata_schema = old_schema
            assert "chromosomes" not in tables.metadata
            # if chromosomes was in metadata
            # we should use its length to get num_chroms,
            # but old file versions did not have this.
            num_chroms = 1
            new_schema = slim_node_metadata_schema(num_chroms)
            tables.nodes.metadata_schema = new_schema
            not_vacant = [0]  # single chromosome
            yes_vacant = [1]
            gt = None
            for n in nodes:
                md = n.metadata
                if len(md) > 0:
                    md['is_vacant'] = yes_vacant if md['is_null'] else not_vacant
                    if not md['is_null']:
                        if gt is None:
                            gt = md['genome_type']
                        else:
                            assert md['genome_type'] == gt, ("Inconsistent tables: "
                                            f"mismatching genome types {gt} and "
                                            f"{md['genome_type']} in node metadata.")
                    del md['is_null']
                    del md['genome_type']
                tables.nodes.append(n.replace(metadata=md))
            # flags for node genome type pre-0.9:
            # confusingly and sub-optimally, these were redundant:
            # all non-null nodes in the same sim would have the same genome type
            assert gt is not None
            GENOME_TYPE_AUTOSOME = 0
            GENOME_TYPE_X = 1
            GENOME_TYPE_Y = 2
            top_md = tables.metadata
            i = top_md['SLiM']['this_chromosome']['index']
            # 'chromosomes' is not required so won't be inserted,
            # so we won't update it also
            assert "chromosomes" not in top_md['SLiM']['this_chromosome'], (""
                    "This is an unexpected result: if you hit this, "
                    "please file a bug at "
                    "https://github.com/tskit-dev/pyslim.")
            if gt == GENOME_TYPE_X:
                top_md['SLiM']['this_chromosome']['type'] = "X"
                top_md['SLiM']['this_chromosome']['symbol'] = "X"
                # if "chromosomes" in tables.metadata['SLiM']:
                #     top_md['SLiM']['chromosomes'][i]['type'] = "X"
                #     top_md['SLiM']['chromosomes'][i]['symbol'] = "X"
            elif gt == GENOME_TYPE_Y:
                top_md['SLiM']['this_chromosome']['type'] = "-Y"
                top_md['SLiM']['this_chromosome']['symbol'] = "Y"
                #     top_md['SLiM']['chromosomes'][i]['type'] = "Y"
                #     top_md['SLiM']['chromosomes'][i]['symbol'] = "Y"
            else:
                assert gt == GENOME_TYPE_AUTOSOME
            tables.metadata = top_md

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
            new_schema = slim_metadata_schemas["individual"]
            tables.individuals.metadata_schema = new_schema
            defaults = default_slim_metadata("individual")
            d = {}
            for k in ["pedigree_p1", "pedigree_p2"]:
                d[k] = defaults[k]
            for ind in inds:
                md = ind.metadata
                md.update(d)
                tables.individuals.append(ind.replace(metadata=md))

        old_schema = _old_metadata_schema("mutation", file_version)
        if old_schema is not None:
            muts = tables.mutations.copy()
            tables.mutations.clear()
            if muts.metadata_schema == tskit.MetadataSchema(None):
                muts.metadata_schema = old_schema
            tables.mutations.metadata_schema = slim_metadata_schemas["mutation"]
            for mut in muts:
                md = mut.metadata
                for ml in md['mutation_list']:
                    ml['nucleotide'] = -1
                tables.mutations.append(mut.replace(metadata=md))

        if file_version == "0.1":
            # shift times
            slim_generation = tables.metadata['SLiM']['tick']
            node_times = tables.nodes.time + slim_generation
            tables.nodes.set_columns(
                    flags=tables.nodes.flags,
                    time=node_times,
                    population=tables.nodes.population,
                    individual=tables.nodes.individual,
                    metadata=tables.nodes.metadata,
                    metadata_offset=tables.nodes.metadata_offset)
            migration_times = tables.migrations.time + slim_generation
            tables.migrations.set_columns(
                    left=tables.migrations.left,
                    right=tables.migrations.right,
                    node=tables.migrations.node,
                    source=tables.migrations.source,
                    dest=tables.migrations.dest,
                    time=migration_times)

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

        set_metadata_schemas(tables)
        md = tables.metadata
        md['SLiM']['file_version'] = slim_file_version
        tables.metadata = md
