"""specify which tests to run on each of the recipes in ./test_recipes"""

import os

# possible attributes to simulation scripts are
#  WF, nonWF
#  nucleotides
#  everyone: records everyone ever
#  pedigree: writes out accompanying info file containing the pedigree
#  remembered_early: remembering and saving the ts happens during early
#  multipop: has more than one population
# All files are of the form `tests/test_recipes/{key}`
recipe_specs = {
    "recipe_nonWF.slim":                       {"nonWF": True, "pedigree": True},
    "recipe_long_nonWF.slim":                  {"nonWF": True},
    "recipe_WF.slim":                          {"WF": True, "pedigree": True},
    "recipe_long_WF.slim":                     {"WF": True},
    "recipe_WF_migration.slim":                {"WF": True, "pedigree": True, "multipop": True},
    "recipe_nonWF_early_early.slim":           {"nonWF": True, "pedigree": True, "remembered_early": True},
    "recipe_nonWF_late_early.slim":            {"nonWF": True, "pedigree": True, "remembered_early": True, "begun_late": True},
    "recipe_WF_early_early.slim":              {"WF": True, "pedigree": True, "remembered_early": True},
    "recipe_WF_late_early.slim":               {"WF": True, "pedigree": True, "remembered_early": True, "begun_late": True},
    "recipe_nonWF_early_late.slim":            {"nonWF": True, "pedigree": True},
    "recipe_nonWF_late_late.slim":             {"nonWF": True, "pedigree": True, "begun_late": True},
    "recipe_WF_early_late.slim":               {"WF": True, "pedigree": True},
    "recipe_WF_late_late.slim":                {"WF": True, "pedigree": True, "begun_late": True},
    "recipe_nucleotides_WF.slim":              {"WF": True, "pedigree": True, "nucleotides": True},
    "recipe_nucleotides_nonWF.slim":           {"nonWF": True, "pedigree": True, "nucleotides": True},
    "recipe_nucleotides_plus_others.slim":     {"WF": True, "pedigree": True, "nucleotides": True, "non-nucleotides": True, "adds_mutations": True},
    "recipe_long_nucleotides.slim":            {"WF": True, "nucleotides": True},
    "recipe_roots.slim":                       {"WF": True, "pedigree": True},
    "recipe_nonWF_selfing.slim":               {"nonWF": True, "pedigree": True},
    "recipe_nonWF_nonstacked.slim":            {"nonWF": True, "nonstacked": True},
    "recipe_init_mutated_WF.slim":             {"WF": True, "init_mutated": True},
    "recipe_init_mutated_nonWF.slim":          {"nonWF": True, "init_mutated": True},
    "recipe_with_metadata.slim":               {"user_metadata": True},
    "recipe_record_everyone_WF_late.slim":     {"WF": True, "pedigree": True, "everyone": True},
    "recipe_record_everyone_WF_early.slim":    {"WF": True, "pedigree": True, "everyone": True, "remembered_early": True},
    "recipe_retain_everyone_WF_late.slim":     {"WF": True, "pedigree": True, "retained": True},
    "recipe_retain_sometimes_WF_late.slim":    {"WF": True, "pedigree": True, "retained": True},
    "recipe_retain_unary_WF_late.slim":        {"WF": True, "pedigree": True, "retained": True, "retainCoalescentOnly": False},
    "recipe_record_everyone_nonWF_late.slim":  {"nonWF": True, "pedigree": True, "everyone": True},
    "recipe_record_everyone_nonWF_early.slim": {"nonWF": True, "pedigree": True, "everyone": True, "remembered_early": True},
    "recipe_retain_everyone_nonWF_late.slim":  {"nonWF": True, "pedigree": True, "retained": True},
    "recipe_retain_sometimes_nonWF_late.slim": {"nonWF": True, "pedigree": True, "retained": True},
    "recipe_retain_unary_nonWF_late.slim":     {"nonWF": True, "pedigree": True, "retained": True, "retainCoalescentOnly": False},
    "recipe_remember_and_retain.slim":         {"nonWF": True, "pedigree": True, "retained": True},
}

def recipe_eq(*keys, exclude=None):
    """
    Return an iterator over those recipes whose spec contains the specified keys.
    If key is empty, return all of them.
    If exclude is given exclude recipes with the specified keys.
    """
    if exclude is None:
        return (
            k for k, v in recipe_specs.items()
            if all(kk in v for kk in keys)
        )
    else:
        return (
            k for k, v in recipe_specs.items()
            if (exclude not in v and all(kk in v for kk in keys))
        )

# These SLiM scripts read in an existing trees file; the "input" gives a key in the
# recipe_specs array that will produce a "ts" file suitable for input
restarted_recipe_specs = {
    "restart_nucleotides_WF.slim":   {"WF": True, "nucleotides": True, "no_op": True, "input": "recipe_nucleotides_WF.slim"},
    "restart_nucleotides_nonWF.slim":   {"nonWF": True, "nucleotides": True, "no_op": True, "input": "recipe_nucleotides_nonWF.slim"},
    "restart_and_run_WF.slim":    {"WF": True, "input": "recipe_init_mutated_WF.slim"},
    "restart_and_run_nonWF.slim": {"nonWF": True, "input": "recipe_init_mutated_nonWF.slim"},
}
for t in ("WF", "nonWF"):
    # recipes that read in and write out immediately ("no_op")
    value = {t: True, "no_op": True, "input": f"recipe_{t}.slim"}
    restarted_recipe_specs[f'restart_{t}.slim'] = value

def restarted_recipe_eq(*keys):
    """
    Return an iterator over a tuple of (restarted_recipe_name, input_recipe_name) giving
    those restarted recipes whose spec contains the specified keys.
    If key is empty, return all of them.
    """
    for k, v in restarted_recipe_specs.items():
        if all(kk in v for kk in keys):
            assert v["input"] in recipe_specs  # we need a valid input recipe
            yield (k, v["input"])

