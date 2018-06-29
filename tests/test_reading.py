import pyslim
import msprime

ts = pyslim.load("simple.trees", slim_format=True)
tables = ts.tables

# mutations

mut_metadata = []
for md in msprime.unpack_bytes(tables.mutations.metadata, 
                               tables.mutations.metadata_offset):
    dm = pyslim.decode_mutation(md)
    edm = pyslim.encode_mutation(dm)
    assert(md == edm)
    mut_metadata.append(dm)

pyslim.annotate_mutation_metadata(tables, mut_metadata)

# nodes

node_metadata = []
for md in msprime.unpack_bytes(tables.nodes.metadata, 
                               tables.nodes.metadata_offset):
    dn = pyslim.decode_node(md)
    edn = pyslim.encode_node(dn)
    assert(md == edn)
    node_metadata.append(dn)

pyslim.annotate_node_metadata(tables, node_metadata)

# individuals

individual_metadata = []
for md in msprime.unpack_bytes(tables.individuals.metadata,
                               tables.individuals.metadata_offset):
    di = pyslim.decode_individual(md)
    edi = pyslim.encode_individual(di)
    assert(md == edi)
    individual_metadata.append(di)

pyslim.annotate_individual_metadata(tables, individual_metadata)

# populations

population_metadata = []
for md in msprime.unpack_bytes(tables.populations.metadata, 
                               tables.populations.metadata_offset):
    if len(md) > 0:
        dp = pyslim.decode_population(md)
        edp = pyslim.encode_population(dp)
        assert(md == edp)
    else:
        dp = None
    population_metadata.append(dp)

pyslim.annotate_population_metadata(tables, population_metadata)

