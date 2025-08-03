"""specify which tests to run on each of the recipes in ./test_recipes"""

import os

# possible attributes to simulation scripts are
#  (the script name itself, to specify a particular one)
#  WF, nonWF
#  adds_mutations
#  nucleotides, non-nucleotides: has the respective sorts of mutations
#  mutation_spectrum: writes out info file on mutation spectrum
#  everyone: records everyone ever
#  pedigree: writes out accompanying info file containing the pedigree
#  remembered_early: remembering and saving the ts happens during early
#  retained: has retained individuals
#  multipop: has more than one population
#  multichrom: has more than one chromosome
#  long: kinda big
#  (chromosome type)
# All files are of the form `tests/test_recipes/{key}`
recipe_specs = {
    "recipe_nonWF.slim":                       {"nonWF": True, "pedigree": True},
    "recipe_nonWF_X.slim":                     {"nonWF": True, "pedigree": True, "X": True},
    "recipe_nonWF_Y.slim":                     {"nonWF": True, "pedigree": True, "Y": True},
    "recipe_nonWF_H.slim":                     {"nonWF": True, "pedigree": True, "H": True},
    "recipe_WF_X.slim":                        {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "X": True},
    "recipe_WF_Y.slim":                        {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "Y": True},
    "recipe_WF_H.slim":                        {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "H": True},
    "recipe_WF_Z.slim":                        {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "Z": True},
    "recipe_WF_W.slim":                        {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "W": True},
    "recipe_WF_HF.slim":                       {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "HF": True},
    "recipe_WF_HM.slim":                       {"WF": True, "pedigree": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True, "HM": True},
    "recipe_long_nonWF.slim":                  {"nonWF": True, "long": True},
    "recipe_old_nonWF.slim":                   {"nonWF": True, "remembered_first": True},
    "recipe_WF.slim":                          {"WF": True, "pedigree": True},
    "recipe_long_WF.slim":                     {"WF": True, "long": True},
    "recipe_WF_migration.slim":                {"WF": True, "pedigree": True, "multipop": True},
    "recipe_nucleotides_WF.slim":              {"WF": True, "pedigree": True, "nucleotides": True},
    "recipe_nucleotides_nonWF.slim":           {"nonWF": True, "pedigree": True, "nucleotides": True},
    "recipe_nucleotides_plus_others.slim":     {"WF": True, "pedigree": True, "nucleotides": True, "non-nucleotides": True, "adds_mutations": True},
    "recipe_long_nucleotides.slim":            {"WF": True, "nucleotides": True, "long": True, "nucleotides": True},
    "recipe_mutation_spectrum.slim":           {"WF": True, "nucleotides": True, "mutation_spectrum": True, "nucleotides": True},
    "recipe_roots.slim":                       {"WF": True, "pedigree": True, "long": True},
    "recipe_starts_later.slim":                {"nonWF": True, "starts_later": True},
    "recipe_nonWF_selfing.slim":               {"nonWF": True, "pedigree": True},
    "recipe_nonWF_nonstacked.slim":            {"nonWF": True, "nonstacked": True},
    "recipe_init_mutated_WF.slim":             {"WF": True, "init_mutated": True},
    "recipe_init_mutated_nonWF.slim":          {"nonWF": True, "init_mutated": True},
    "recipe_with_metadata.slim":               {"user_metadata": True},
    "recipe_resettable_WF.slim":               {"WF": True, "resettable": True, "user_metadata": True},
    "recipe_resettable_nonWF.slim":            {"nonWF": True, "resettable": True, "user_metadata": True},
    "recipe_retain_everyone_WF_late.slim":     {"WF": True, "pedigree": True, "retained": True},
    "recipe_retain_sometimes_WF_late.slim":    {"WF": True, "pedigree": True, "retained": True},
    "recipe_retain_unary_WF_late.slim":        {"WF": True, "pedigree": True, "retained": True, "retainCoalescentOnly": False},
    "recipe_retain_everyone_nonWF_late.slim":  {"nonWF": True, "pedigree": True, "retained": True},
    "recipe_retain_sometimes_nonWF_late.slim": {"nonWF": True, "pedigree": True, "retained": True},
    "recipe_retain_unary_nonWF_late.slim":     {"nonWF": True, "pedigree": True, "retained": True, "retainCoalescentOnly": False},
    "recipe_remember_and_retain.slim":         {"nonWF": True, "pedigree": True, "retained": True},
    "recipe_all_the_chromosome_types.slim":    {"nonWF": True, "pedigree": True, "multichrom": True, "nucleotides": True, "X": True, "Y": True, "H": True, "Z": True, "W": True, "HF": True, "FL": True, "HM": True, "ML": True, "Y-": True},
    "recipe_chromosomes_adds_muts.slim":       {"WF": True, "pedigree": True, "multichrom": True, "adds_mutations": True, "nucleotides": True, "non-nucleotides": True},
    "recipe_many_chromosomes.slim":            {"nonWF": True, "pedigree": True, "multichrom": True},
    "recipe_H-_chromosome.slim":               {"nonWF": True, "pedigree": True, "multichrom": True, "H-": True},
}

for x in ("first", "early", "late"):
    for y in ("first", "early", "late"):
        for t in ("WF", "nonWF"):
            d = {"pedigree": True}
            d[t] = True
            if y != "late":
                # "remembered_early" or "remembered_first"
                d[f"remembered_{y}"] = True
            if x != "early":
                # "begun_late" or "begun_first"
                d[f"begun_{x}"] = True
            recipe_specs[f"recipe_{t}_{x}_{y}.slim"] = d

for y in ("first", "early", "late"):
    for t in ("WF", "nonWF"):
        d = {"pedigree": True, "everyone": True}
        d[t] = True
        if y != "late":
            # "remembered_early" or "remembered_first"
            d[f"remembered_{y}"] = True
        recipe_specs[f"recipe_record_everyone_{t}_{y}.slim"] = d

# add the script name itself to the keys: this must come last!
for k in recipe_specs:
    recipe_specs[k][k] = True

def recipe_eq(*keys, exclude=None):
    """
    Return an iterator over those recipes whose spec contains the specified keys.
    If key is empty, return all of them.
    If exclude is given exclude recipes with the specified keys.
    """
    if isinstance(exclude, str) or exclude is None:
        exclude = [exclude]
    return (
        k for k, v in recipe_specs.items()
        if (all(u not in v for u in exclude) and all(kk in v for kk in keys))
    )

# These SLiM scripts read in an existing trees file; the "input" gives a key in the
# recipe_specs array that will produce a "ts" file suitable for input
restarted_recipe_specs = {
    "restart_nucleotides_WF.slim":   {
        "recipe_nucleotides_WF.slim": {"WF": True, "nucleotides": True, "no_op": True},
    },
    "restart_nucleotides_nonWF.slim":   {
        "recipe_nucleotides_nonWF.slim": {"nonWF": True, "nucleotides": True, "no_op": True},
    },
    "restart_and_run_WF.slim":    {
        "recipe_init_mutated_WF.slim": {"WF": True}, 
    },
    "restart_and_run_nonWF.slim": {
        "recipe_init_mutated_nonWF.slim": {"nonWF": True},
    },
    "restart_and_remove_subpop_nonWF.slim":    {
        "recipe_init_mutated_nonWF.slim": {"nonWF": True, "remove_subpop": True},
    },
    # recipes that read in and write out immediately ("no_op")
    "restart_WF.slim": {
        "recipe_WF.slim": {"WF": True, "no_op": True},
        "recipe_resettable_WF.slim": {"WF": True, "no_op": True, "resettable": True},
    },
    "restart_nonWF.slim": {
        "recipe_nonWF.slim": {"nonWF": True, "no_op": True},
        "recipe_resettable_nonWF.slim": {"nonWF": True, "no_op": True, "resettable": True},
    },
}

# add the script name itself to the keys: this must come last!
for k in restarted_recipe_specs:
    for j in restarted_recipe_specs[k]:
        d = restarted_recipe_specs[k][j]
        d[k] = True
        d[j] = True

def restarted_recipe_eq(*keys):
    """
    Return an iterator over a tuple of (restarted_recipe_name, input_recipe_name) giving
    those restarted recipes whose spec contains the specified keys.
    If key is empty, return all of them.
    """
    for k, v in restarted_recipe_specs.items():
        for ik, iv in v.items():
            if all(kk in iv for kk in keys):
                assert ik in recipe_specs  # we need a valid input recipe
                yield (k, ik)

