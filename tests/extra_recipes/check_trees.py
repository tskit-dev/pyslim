#!/usr/bin/env python3

import sys
import tskit

treefile = sys.argv[1] 

ts = tskit.load(treefile)
print(f"The file '{treefile}' has {ts.num_trees} trees relating {ts.num_individuals} individuals "
      f"with {ts.num_mutations} mutations at {ts.num_sites} sites.")
