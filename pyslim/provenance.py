from __future__ import print_function

import platform
import warnings
import attr
import json
import msprime
import tskit

from . import _version

__version__ = _version.pyslim_version


@attr.s
class ProvenanceMetadata(object):
    model_type = attr.ib()
    slim_generation = attr.ib()
    file_version = attr.ib()


def get_provenance(ts):
    '''
    Extracts model type, slim generation, and remembmered node count from the last
    entry in the provenance table that is tagged with "program"="SLiM".

    :param SlimTreeSequence: The tree sequence.
    :rtype ProvenanceMetadata:
    '''
    prov = [json.loads(x.record) for x in ts.tables.provenances]
    prov_info = [(_slim_provenance_version(u), u) for u in prov]
    slim_prov = [x for x in prov_info if x[0][0]]
    if len(slim_prov) == 0:
        raise ValueError("Tree sequence contains no SLiM provenance entries.")
    info, record = slim_prov[len(slim_prov)-1]
    file_version = info[1]
    if file_version == "0.1":
        out = ProvenanceMetadata(record['model_type'],
                                 record['generation'],
                                 file_version)
    else: # >= 0.2
        out = ProvenanceMetadata(record['parameters']['model_type'],
                                 record['slim']["generation"],
                                 file_version)
    return out


def upgrade_slim_provenance(tables):
    """
    Copies the last provenance entry from a previous SLiM file version to that
    required by the current file version.

    :param TableCollection tables: the table collection
    """
    provlist = [json.loads(x.record) for x in tables.provenances]
    prov_info = [(_slim_provenance_version(u), u) for u in provlist]
    slim_prov = [x for x in prov_info if x[0][0]]
    if len(slim_prov) == 0:
        raise ValueError("Tree sequence contains no SLiM provenance entries.")
    info, record = slim_prov[len(slim_prov)-1]
    file_version = info[1]
    if not (file_version == "0.1" or file_version == "0.2"):
        warnings.warn("File version is not v0.1 or v0.2; not doing anything.")
    is_slim, version = _slim_provenance_version(record)
    if not is_slim:
        raise ValueError("Not a SLiM provenance entry.")
    if file_version == "0.1":
        new_record = make_slim_provenance_dict(
                        record['model_type'],
                        record['generation'])
        new_record['parameters']['command'] = ['pyslim', 'convert']
    else:
        new_record = make_slim_provenance_dict(
                        record['parameters']['model_type'],
                        record['slim']['generation'])
        new_record['parameters']['command'] = ['pyslim', 'convert']
    tskit.validate_provenance(new_record)
    tables.provenances.add_row(json.dumps(new_record))


def _slim_provenance_version(record):
    software_name = "unknown"
    file_version = "unknown"
    # >= SLiM 3.1 // file version >= 0.2
    try:
        software_name = record["software"]["name"]
    except:
        software_name = "unknown"

    if software_name == "SLiM":
        try:
            file_version = record["slim"]["file_version"]
        except:
            pass
    else:
        # SLiM 3.0 // file version 0.1
        try:
            software_name = record["program"]
        except:
            pass
        try:
            file_version = record["file_version"]
        except:
            pass
    is_slim = (software_name == "SLiM") and (file_version in ["0.1", "0.2", "0.3"])
    return is_slim, file_version


def get_environment():
    """
    Returns a dictionary describing the environment in which msprime
    is currently running.
    """
    env = {
        "libraries": {
        },
        "parameters" : {
            "command" : []
        },
        "os": {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version_tuple(),
        }
    }
    return env


def make_pyslim_provenance_dict():
    """
    Returns a dictionary encoding the information about this version of pyslim.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name" : "pyslim",
            "version": __version__,
            },
        "parameters": {
            "command": {}
            },
        "environment": get_environment()
    }
    return document

def make_slim_provenance_dict(model_type, slim_generation):
    """
    Returns a dictionary encoding necessary provenance information for a SLiM tree sequence.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name" : "SLiM",
            "version": "3.3"
            },
        "parameters": {
            "command": ['pyslim'],
            "model_type": model_type,
            },
        "environment": {},
        "metadata": {
            "individuals": {
                "flags": {
                    "16": {
                        "name" : "SLIM_TSK_INDIVIDUAL_ALIVE",
                        "description" : "the individual was alive "
                              + "at the time the file was written",
                          },
                    "17": {
                        "name" : "SLIM_TSK_INDIVIDUAL_REMEMBERED",
                        "description" : "the individual was requested "
                              + "by the user to be remembered",
                          },
                    "18": {
                        "name" : "SLIM_TSK_INDIVIDUAL_FIRST_GEN",
                        "description" : "the individual was in the first "
                              + "generation of a new population"
                          }
                }
            }
        },
        "slim": {
            "file_version": "0.3",
            "generation": slim_generation,
            "model": ""
            }
    }
    return document

_slim_v3_1_example = '''
{
    "environment": {
        "os": {
            "machine": "x86_64",
            "node": "d93-172.uoregon.edu",
            "release": "17.6.0",
            "system": "Darwin",
            "version": "Darwin Kernel Version 17.6.0: Tue May  8 15:22:16 PDT 2018; root:xnu-4570.61.1~1/RELEASE_X86_64"
        }
    },
    "metadata": {
        "individuals": {
            "flags": {
                "16": {
                    "description": "the individual was alive at the time the file was written",
                    "name": "SLIM_TSK_INDIVIDUAL_ALIVE"
                },
                "17": {
                    "description": "the individual was requested by the user to be remembered",
                    "name": "SLIM_TSK_INDIVIDUAL_REMEMBERED"
                },
                "18": {
                    "description": "the individual was in the first generation of a new population",
                    "name": "SLIM_TSK_INDIVIDUAL_FIRST_GEN"
                }
            }
        }
    },
    "parameters": {
        "command": [],
        "model": "initialize() {\n\tinitializeTreeSeq();\n\tinitializeMutationRate(1e-7);\n\tinitializeMutationType(\"m1\", 0.5, \"f\", 0.0);\n\tinitializeGenomicElementType(\"g1\", m1, 1.0);\n\tinitializeGenomicElement(g1, 0, 99999);\n\tinitializeRecombinationRate(1e-8);\n}\n1 {\n\tsim.addSubpop(\"p1\", 500);\n}\n2000 late() { sim.treeSeqOutput(\"~/Desktop/junk.trees\"); }\n",
        "model_type": "WF",
        "seed": 1783301962445
    },
    "schema_version": "1.0.0",
    "slim": {
        "file_version": "0.2",
        "generation": 2000
    },
    "software": {
        "name": "SLiM",
        "version": "3.1"
    }
}
'''

_slim_v3_0_example = '''
{'id': 0, 'timestamp': '2018-08-25T14:59:13', 'record': '{"program": "SLiM", "version": "3.0", "file_version": "0.1", "model_type": "WF", "generation": 10, "remembered_node_count": 0}'}
'''
