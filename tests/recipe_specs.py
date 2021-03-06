"""specify which tests to run on each of the recipes in ./examples"""

import os

# possible attributes to simulation scripts are
#  WF, nonWF
#  nucleotides
#  everyone: records everyone ever
#  pedigree: writes out accompanying info file containing the pedigree
#  remembered_early: remembering and saving the ts happens during early
#  multipop: has more than one population
# All files are of the form `tests/examples/{key}.slim`
basic_recipe_specs = {
    "recipe_nonWF.slim":              {"nonWF": True, "pedigree": True},
    "recipe_long_nonWF.slim":         {"nonWF": True},
    "recipe_WF.slim":                 {"WF": True, "pedigree": True},
    "recipe_long_WF.slim":            {"WF": True},
    "recipe_WF_migration.slim":       {"WF": True, "pedigree": True, "multipop": True},
    "recipe_nonWF_early.slim":        {"nonWF": True, "pedigree": True, "remembered_early": True},
    "recipe_WF_early.slim":           {"WF": True, "pedigree": True, "remembered_early": True},
    "recipe_nucleotides.slim":        {"WF": True, "pedigree": True, "nucleotides": True},
    "recipe_long_nucleotides.slim":   {"WF": True, "nucleotides": True},
    "recipe_roots.slim":              {"WF": True, "pedigree": True},
    "recipe_nonWF_selfing.slim":      {"nonWF": True, "pedigree": True},
    "recipe_init_mutated_WF.slim":    {"WF": True, "init_mutated": True},
    "recipe_init_mutated_nonWF.slim": {"nonWF": True, "init_mutated": True},
    "recipe_with_metadata.slim":      {"user_metadata": True},
}

# Some more complicated ones
for t in ("WF", "nonWF"):
    for s in ("early", "late"):
        value = {t: True, "everyone": True, "pedigree": True}
        if s == 'early':
            value['remembered_early'] = True
        basic_recipe_specs[f'recipe_record_everyone_{t}_{s}.slim'] = value

def basic_recipe_eq(*keys, exclude_if_has_key=None):
    """
    Return an iterator over those recipes whose spec contains the specified keys.
    If key is empty, return all of them.
    If exclude_key is given exclude recipes with the specified key
    """
    if exclude_if_has_key is None:
        return (k for k, v in basic_recipe_specs.items() if all(kk in v for kk in keys))
    else:
        return (
            k for k, v in basic_recipe_specs.items()
            if exclude_if_has_key not in v and all(kk in v for kk in keys)
        )

# These SLiM scripts read in an existing trees file; the "input" gives a key in the
# basic_recipe_specs array that will produce a "ts" file suitable for input
restarted_recipe_specs = {
    "restart_nucleotides.slim":   {"WF": True, "nucleotides": True, "no_op": True, "input": "recipe_nucleotides.slim"},
    #"restart_and_run_WF.slim":    {"WF": True, "input": "recipe_init_mutated.slim"},
    #"restart_and_run_nonWF.slim": {"nonWF": True, "input": "recipe_init_mutated.slim"},
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
            assert v["input"] in basic_recipe_specs  # we need a valid input recipe
            yield (k, v["input"])

