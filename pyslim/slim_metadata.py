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
                    },
                ],
                "flags": {
                    "SLIM_INDIVIDUAL_METADATA_MIGRATED": {
                        "description": "Whether this individual was a migrant, either in the generation when the tree sequence was written out "
                                       "(if the individual was alive then), or in the generation of the last time they were Remembered (if not).",
                        "value": 1
                    }
                },
                "properties": {
                    "age": {
                        "binaryFormat": "i",
                        "description": "The age of this individual, either when the tree sequence "
                                       "was written out (if the individual was alive then), or the "
                                       "last time they were Remembered (if not).",
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
                        "description": "The sex of the individual (0 for female, 1 for male, "
                                       "-1 for hermaphrodite).",
                        "index": 6,
                        "type": "integer"
                    },
                    "subpopulation": {
                        "binaryFormat": "i",
                        "description": "The ID of the subpopulation the individual was part of, "
                                       "either when the tree sequence was written out (if the "
                                       "individual was alive then), or the last time they were "
                                       "Remembered (if not).",
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
                "type": "object",
        },
        "population": 
            {"$schema":"http://json-schema.org/schema#","additionalProperties":True,"codec":"json","description":"SLiM schema for population metadata.","examples":[{"bounds_x0":0.0,"bounds_x1":100.0,"bounds_y0":0.0,"bounds_y1":100.0,"female_cloning_fraction":0.25,"male_cloning_fraction":0.0,"migration_records":[{"migration_rate":0.9,"source_subpop":1},{"migration_rate":0.1,"source_subpop":2}],"selfing_fraction":0.5,"sex_ratio":0.5,"slim_id":2,"name":"p2"}],"properties":{"bounds_x0":{"description":"The minimum x-coordinate in this subpopulation.","type":"number"},"bounds_x1":{"description":"The maximum x-coordinate in this subpopulation.","type":"number"},"bounds_y0":{"description":"The minimum y-coordinate in this subpopulation.","type":"number"},"bounds_y1":{"description":"The maximum y-coordinate in this subpopulation.","type":"number"},"bounds_z0":{"description":"The minimum z-coordinate in this subpopulation.","type":"number"},"bounds_z1":{"description":"The maximum z-coordinate in this subpopulation.","type":"number"},"description":{"description":"A description of this subpopulation.","type":"string"},"female_cloning_fraction":{"description":"The frequency with which females in this subpopulation reproduce clonally (for WF models).","type":"number"},"male_cloning_fraction":{"description":"The frequency with which males in this subpopulation reproduce clonally (for WF models).","type":"number"},"migration_records":{"items":{"properties":{"migration_rate":{"description":"The fraction of children in this subpopulation that are composed of 'migrants' from the source subpopulation (in WF models).","type":"number"},"source_subpop":{"description":"The ID of the subpopulation migrants come from (in WF models).","type":"integer"}},"required":["source_subpop","migration_rate"],"type":"object"},"type":"array"},"name":{"description":"A human-readable name for this subpopulation.","type":"string"},"selfing_fraction":{"description":"The frequency with which individuals in this subpopulation self (for WF models).","type":"number"},"sex_ratio":{"description":"This subpopulation's sex ratio (for WF models).","type":"number"},"slim_id":{"description":"The ID of this population in SLiM. Note that this is called a 'subpopulation' in SLiM.","type":"integer"}},"required":["slim_id"],"type":["object","null"]}
    }


slim_metadata_schemas = {k: tskit.MetadataSchema(_raw_slim_metadata_schemas[k])
        for k in _raw_slim_metadata_schemas}


def default_slim_metadata(name):
    '''
    Returns default metadata of type ``name``, where ``name`` is one of
    "tree_sequence", "edge", "site", "mutation", "mutation_list_entry",
    "node", "individual", or "population".

    :param str name: The type of metadata requested.
    :rtype dict:
    '''
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
        generation,
        spatial_dimensionality='',
        spatial_periodicity='',
        separate_sexes=False,
        nucleotide_based=False,
        stage='late',
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


def set_metadata_schemas(tables):
    tables.edges.metadata_schema = slim_metadata_schemas['edge']
    tables.sites.metadata_schema = slim_metadata_schemas['site']
    tables.mutations.metadata_schema = slim_metadata_schemas['mutation']
    tables.nodes.metadata_schema = slim_metadata_schemas['node']
    tables.individuals.metadata_schema = slim_metadata_schemas['individual']
    tables.populations.metadata_schema = slim_metadata_schemas['population']


################################
# Previous versions of metadata schema:

def _old_metadata_schema(name, file_version):
    # Returns a metadata schema *if the format has changed*,
    # and None otherwise.
    ms = None
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

    # everything else's format has remained unchanged
    if ms is not None:
        ms = tskit.MetadataSchema(ms)
    return ms
