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
        df = pd.read_csv(logfile, skiprows=count)
        df = df.set_index('generation')
        return df

    @cell_magic
    def slim_stats_reps(self, reps, cell):
        script = cell
        n = int(reps)
        aList = []
        for i in range(n):
            logfile = "tmp.log"
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
            df = pd.read_csv(logfile, skiprows=count)
            df = df.set_index('generation')
            # TODO generalize me to more statistics/summaries
            df.rename(columns={'H':'H'+str(i)}, inplace=True)
            aList.append(df)
        dff = pd.concat(aList, axis=1, join="inner")
        return dff
        
    @cell_magic
    def slim_stats_reps_stack(self, reps, cell):
        """
        slim_stats_reps_stack returns a pandas df in which
        reps number of replicate slim simulations have been
        run, and their output captured, and then stacked, merging
        along columns. 
        
        output from simulation is expected to be a comma-delimited
        list of summaries printed from slim to stdout. by convention
        the header row begins with 'generation' e.g.,
        
        generation,stat1,stat2,...,statn
        """
        script = cell
        n = int(reps)
        aList = []
        for i in range(n):
            logfile = "tmp.log"
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
            df = pd.read_csv(logfile, skiprows=count)
            df = df.set_index('generation')
            aList.append(df)
        dff = pd.concat(aList)
        return dff
        

