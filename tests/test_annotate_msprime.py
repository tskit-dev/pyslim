import msprime
import pyslim

ts = msprime.simulate(10, mutation_rate = 1.0, recombination_rate = 1.0)
slim_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
slim_ts.dump("new_slim.trees")

slim_ts2 = pyslim.load("new_slim.trees", slim_format=True)
