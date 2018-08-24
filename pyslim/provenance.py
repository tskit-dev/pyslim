
from __future__ import print_function

import platform

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass

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


def get_provenance_dict():
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

def make_slim_dict(model_type, slim_generation):
    """
    Returns a dictionary encoding necessary provenance information for a SLiM tree sequence.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name" : "SLiM",
            "version": "3.0"
            },
        "parameters": {
            "command": ['pyslim']
            },
        "environment": {},
        "slim": {
            "file_version": "0.2",
            "generation": str(slim_generation),
            "model": "",
            "model_type": model_type
            }
    }
    return document

## AN EXAMPLE:
# {
#     "environment": {
#         "os": {
#             "machine": "x86_64",
#             "node": "darwin-447.local",
#             "release": "17.6.0",
#             "system": "Darwin",
#             "version": "Darwin Kernel Version 17.6.0: Tue May  8 15:22:16 PDT 2018; root:xnu-4570.61.1~1/RELEASE_X86_64"
#         }
#     },
#     "parameters": {
#         "command": []
#     },
#     "schema_version": "1.0.0",
#     "slim": {
#         "file_version": "0.2",
#         "generation": 2000,
#         "model": "initialize() {\n\tinitializeTreeSeq();\n\tinitializeMutationRate(1e-7);\n\tinitializeMutationType(\"m1\", 0.5, \"f\", 0.0);\n\tinitializeGenomicElementType(\"g1\", m1, 1.0);\n\tinitializeGenomicElement(g1, 0, 99999);\n\tinitializeRecombinationRate(1e-8);\n}\n1 {\n\tsim.addSubpop(\"p1\", 500);\n}\n2000 late() { sim.treeSeqOutput(\"~/Desktop/junk.trees\"); }\n",
#         "model_type": "WF",
#         "remembered_node_count": 0,
#         "seed": 1722162964723
#     },
#     "software": {
#         "name": "SLiM",
#         "version": "3.0"
#     }
# }
