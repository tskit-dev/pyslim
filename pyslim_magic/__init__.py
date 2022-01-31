"""pyslim magic"""
__version__ = '0.0.1'

from .pyslim_magic import PySlimMagic

def load_ipython_extension(ipython):
    ipython.register_magics(PySlimMagic)