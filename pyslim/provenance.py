from __future__ import print_function

import platform
import warnings
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
