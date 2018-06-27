import pyslim
import msprime

ts = msprime.simulate(10, mutation_rate = 1.0, recombination_rate = 1.0)
tables = ts.tables

print(tables)
print("---Individuals---\n")
print(tables.individuals)
print("---Populations---\n")
print(tables.populations)

pyslim.set_nodes_individuals(tables)

pyslim.set_populations(tables)

print("---Populations---\n")
print(tables.populations)

pyslim.set_mutations(tables)

print(tables)

new_ts = pyslim.SlimTreeSequence.load_tables(tables)

for t in new_ts.trees():
    print(t)
