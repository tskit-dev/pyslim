__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass

def get_provenance_dict():
    """
    Returns a dictionary encoding the version of pyslim used.
    """
    document = {
        "software": "pyslim",
        "version": __version__
    }
    return document


