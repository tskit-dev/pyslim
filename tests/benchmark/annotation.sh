#!/bin/bash

SETUP_SETUP="
import msprime

def setup():
    ts = msprime.simulate(
            2000,
            Ne=1000,
            length=1e8,
            discrete_genome=True, 
            recombination_rate=1e-8)
    ts.dump('to_annotate.trees')
    return 0
"

SETUP="
import pyslim, tskit
ts = tskit.load('to_annotate.trees')

def pyslim_fn(ts):
    tables = ts.tables
    pyslim.annotate_defaults_tables(tables, model_type='WF', slim_generation=1)
    return 0

def tskit_fn(ts):
    tables = ts.tables
    dmd = {}
    for k in ('node', 'individual', 'population'):
        md = pyslim.default_slim_metadata(k)
        ms = pyslim.slim_metadata_schemas[k]
        dmd[k] = ms.validate_and_encode_row(md)

    pop_md = [b'' for _ in tables.populations]
    for p in set(tables.nodes.population):
        pop_md[p] = dmd['population']
    tables.populations.packset_metadata(pop_md)

    node_md = [b'' if (n.flags & 1) else dmd['node'] for n in tables.nodes]
    tables.nodes.packset_metadata(node_md)

    ind_md = [dmd['individual']] * tables.individuals.num_rows
    tables.individuals.packset_metadata(ind_md)
    return 0
    
"
echo "setup:"
python3 -m timeit -s "$SETUP_SETUP" "setup()"

echo "tskit:"
python3 -m timeit -s "$SETUP" "tskit_fn(ts)"

echo "pyslim:"
python3 -m timeit -s "$SETUP" "pyslim_fn(ts)"
