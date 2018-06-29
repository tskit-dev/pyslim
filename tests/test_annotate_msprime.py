import msprime
import pyslim

ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
slim_ts = pyslim.annotate(ts, model_type="nonWF", slim_generation=1)
pyslim.dump(slim_ts, "new_slim.trees", slim_format=True)

slim_ts2 = pyslim.load("new_slim.trees", slim_format=True)
