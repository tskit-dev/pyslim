import pyslim
import msprime
import random

# Default annotation
ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
new_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)

# Annotate, then modify
ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
new_ts = pyslim.annotate_defaults(ts, model_type="nonWF", slim_generation=1)
tables = new_ts.tables

individual_metadata = list(pyslim.extract_individual_metadata(tables))

for j in range(len(individual_metadata)):
    individual_metadata[j].sex = random.choice([0, 1])

pyslim.annotate_individual_metadata(tables, individual_metadata)

slim_ts = pyslim.load_tables(tables, slim_format=True)
