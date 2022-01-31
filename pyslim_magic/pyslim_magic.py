from IPython.core.magic import (Magics, magics_class, line_magic, cell_magic)
from io import StringIO
import os
import subprocess
import pandas as pd


@magics_class
class PySlimMagic(Magics):

    @cell_magic
    def slim_stats(self, name, cell):
        script = cell
        logfile = name + ".log"
        os.system("echo '" + script + "' | slim > " + logfile)
        # deal with ugliness of slim output and cut off
        # the top few header lines
        pattern = "generation,"
        count = 0
        with open(logfile, "r") as f:
            for line in f:
                if pattern in line:
                    break
                else:
                    count += 1
        return pd.read_csv(logfile, skiprows=count)


