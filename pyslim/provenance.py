from __future__ import print_function

import platform
import warnings
import attr
import json
import msprime
import tskit

from . import _version

__version__ = _version.pyslim_version


def slim_provenance_version(provenance):
    """
    Parses a provenance record, returning whether the record is a SLiM
    provenance entry, and version is the file format version, or "unknown" if
    it is not a SLiM entry.

    :param Provenance provenance: The provenance entry, as for instance obtained
        from ts.provenance(0).
    :return: A (bool, string) tuple (is_slim, version).
    """
    record = json.loads(provenance.record)
    software_name = "unknown"
    file_version = "unknown"
    # >= SLiM 3.1 // file version >= 0.2
    try:
        software_name = record["software"]["name"]
    except:
        pass

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
    is_slim = (software_name == "SLiM")
    return is_slim, file_version


def get_environment():
    """
    Returns a dictionary describing the environment in which we are
    currently running.
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

def make_slim_provenance_dict(model_type, slim_generation,
        stage='late', spatial_dimensionality='', spatial_periodicity='',
        separate_sexes=False, nucleotide_based=False):
    """
    Returns a dictionary encoding provenance information for a SLiM tree sequence.
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
            "stage": stage,
            "spatial_dimensionality": spatial_dimensionality,
            "spatial_periodicity": spatial_periodicity,
            "separate_sexes": separate_sexes,
            "nucleotide_based": nucleotide_based,
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
                              + "by the user to be permanently remembered",
                          },
                    "18": {
                        "name" : "SLIM_TSK_INDIVIDUAL_RETAINED",
                        "description" : "the individual was requested "
                              + "by the user to be retained only if its "
                              + "nodes continue to exist in the tree sequence",
                          },
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

