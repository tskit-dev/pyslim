import pyslim
import msprime

ts = msprime.load("tests/simple.trees")

m = ts.mutation(0)
pyslim.decode_mutation(m.metadata)

n = ts.node(0)
pyslim.decode_node(n.metadata)

i = ts.individual(0)
pyslim.decode_individual(i.metadata)

p = ts.population(1)
pyslim.decode_population(p.metadata)
