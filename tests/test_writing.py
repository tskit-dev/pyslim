import pyslim
import msprime

ts = msprime.simulate(10, mutation_rate = 0.0, recombination_rate = 1.0)
tables = ts.tables

print(tables)
print("---Individuals---\n")
print(tables.individuals)
print("---Populations---\n")
print(tables.populations)

pyslim.annotate(tables, model_type="nonWF")

print(tables)
print("---Individuals---\n")
print(tables.individuals)
print("---Populations---\n")
print(tables.populations)

new_ts = pyslim.load_tables(tables)

for t in new_ts.trees():
    print(t)
