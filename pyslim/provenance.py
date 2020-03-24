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
    has_slim_provenance = False
    for p in ts.tables.provenances:
        this_record = json.loads(p.record)
        is_slim, this_file_version = _slim_provenance_version(this_record) 
        if is_slim:
            record = this_record
            file_version = this_file_version
            has_slim_provenance = True

    if not has_slim_provenance:
        raise ValueError("Tree sequence contains no SLiM provenance entries"
                          "(or your pyslim is out of date).")

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
    if not (float(file_version) < 0.4):
        warnings.warn("File version is not older than 0.4; not doing anything.")
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
    is_slim = (software_name == "SLiM") and (file_version in ["0.1", "0.2", "0.3", "0.4"])
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
            "version": "3.3.2"
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
            "file_version": "0.4",
            "generation": slim_generation,
            "model": ""
            }
    }
    return document

