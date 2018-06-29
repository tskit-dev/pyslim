import msprime
import pyslim
import random

ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
ts.dump("msprime.trees")

ts = pyslim.load("msprime.trees", slim_format=False)
# now we will be editing it so we need to switch to tables
tables = ts.tables
pyslim.annotate_tables(tables, model_type="nonWF", slim_generation=1)
individual_metadata = list(pyslim.extract_individual_metadata(tables))
for j in range(len(individual_metadata)):
    individual_metadata[j].sex = random.choice([0, 1])

pyslim.annotate_individual_metadata(tables, individual_metadata)
slim_ts = pyslim.load_tables(tables, slim_format=True)
slim_ts.dump("new_msprime.trees")

slim_ts2 = pyslim.load("new_msprime.trees", slim_format=True)
