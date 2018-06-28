import pyslim
import msprime
import random

ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
tables = ts.tables
pyslim.annotate(tables, model_type="nonWF")

individual_metadata = list(pyslim.extract_individual_metadata(tables))

for j in range(len(individual_metadata)):
    individual_metadata[j].sex = random.choice([0, 1])

pyslim.annotate_individual_metadata(tables, individual_metadata)

slim_ts = pyslim.load_tables(tables)
