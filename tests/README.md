Install the dependencies using

```
pip3 install -r requirements/CI/requirements.txt
```

SLiM recipes used for testing are in `test_recipes/`: the specifications for which
tests to run on them are specified in `recipe_specs.py`

To run the entire test suite: 

```
python3 -m pytest
```

To restrict to individual tests, the syntax is:

```
python3 -m pytest tests/test_tree_sequence.py::test_slim_generation
```


To restrict to specific recipes in the `test_recipes` dir, the syntax is:

```
python3 -m pytest -k recipe_nonWF.slim
```

Note that this can be combined with restricting to specific tests, and the value
following the `-k` is treated as a match, so e.g. `-k recipe_WF` will match all recipes
containing the string `recipe_WF`