#!/bin/bash


SETUP="
import tskit, pyslim
ts = tskit.load('benchmark/metadata.trees')
pts = pyslim.load('benchmark/metadata.trees')
tables = ts.tables
tables.individuals.metadata_schema = tskit.MetadataSchema(None)
tables.nodes.metadata_schema = tskit.MetadataSchema(None)
tables.mutations.metadata_schema = tskit.MetadataSchema(None)
nts = tables.tree_sequence()

def tskit_fn(ts):
    pid = [ind.metadata['pedigree_id'] for ind in ts.individuals()]
    sid = [n.metadata['slim_id'] for n in ts.nodes()]
    s = [m.metadata['mutation_list'][0]['selection_coeff'] for m in ts.mutations()]
    return 0

def pyslim_fn(ts):
    pid = [ind.metadata.pedigree_id for ind in pts.individuals()]
    sid = [n.metadata.slim_id for n in pts.nodes()]
    s = [m.metadata[0].selection_coeff for m in pts.mutations()]
    return 0

def decode_fn(ts):
    pid = [pyslim.decode_individual(ind.metadata).pedigree_id for ind in nts.individuals()]
    sid = [pyslim.decode_node(n.metadata).slim_id for n in nts.nodes()]
    s = [pyslim.decode_mutation(m.metadata)[0].selection_coeff for m in nts.mutations()]
    return 0
"

echo "tskit:"
python3 -m timeit -s "$SETUP" "tskit_fn(ts)"

echo "pyslim:"
python3 -m timeit -s "$SETUP" "pyslim_fn(ts)"

echo "pyslim decode:"
python3 -m timeit -s "$SETUP" "decode_fn(nts)"
