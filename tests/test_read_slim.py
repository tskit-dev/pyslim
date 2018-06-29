import msprime
import pyslim

ts = pyslim.load("simple.trees", slim_format=True)

# selection coefficients and locations of all selected mutations
sel_loci = []
for mut in ts.mutations():
    md = pyslim.decode_mutation(mut.metadata)
    sel = [x.selection_coeff for x in md]
    if md is not None and any([s != 0 for s in sel]):
        sel_loci.append((mut.position, sel))

