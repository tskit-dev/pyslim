import tskit, pyslim

"""
Takes an old tree sequence and update the metadata *without* properly updating
the top-level metadata.
"""

ts = tskit.load("recipe_WF.v3.5.trees")
tables = ts.dump_tables()

tables.populations.clear()
tables.populations.metadata_schema = pyslim.slim_metadata_schemas['population']
for p in ts.populations():
    tables.populations.append(p)

tables.individuals.clear()
tables.individuals.metadata_schema = pyslim.slim_metadata_schemas['individual']
d = pyslim.default_slim_metadata("individual")
for i in ts.individuals():
    d.update(i.metadata)
    ii = i.replace(metadata=d)
    tables.individuals.append(ii)

tables.mutations.clear()
tables.mutations.metadata_schema = pyslim.slim_metadata_schemas['mutation']
d = pyslim.default_slim_metadata("mutation")
for m in ts.mutations():
    tables.mutations.append(m)

ts = tables.tree_sequence()
ts.dump("recipe_WF.v3.5_and_v3.6.trees")
