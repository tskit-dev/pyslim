---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell}
:tags: [remove-cell]
import pyslim, tskit, msprime
from IPython.display import SVG
import numpy as np

ts = tskit.load("example_sim.trees")
tables = ts.tables
```

```{eval-rst}
.. currentmodule:: pyslim
```


(sec_python_api)=

# Python API

This page provides detailed documentation for the methods and classes
available in pyslim.
Here is a quick reference to some of the methods:

```{eval-rst}
.. autosummary::

  recapitate
  annotate
  individuals_alive_at
  individual_ages
  individual_ages_at
  remove_vacant
  restore_vacant
  has_vacant_samples
  node_is_vacant
  slim_time
  convert_alleles
  generate_nucleotides
  population_size
  set_slim_state
  default_slim_metadata
  update
```


## Editing or adding to tree sequences

``pyslim`` provides tools for transforming tree sequences:


```{eval-rst}
.. autofunction:: recapitate
```

```{eval-rst}
.. autofunction::  convert_alleles
```

```{eval-rst}
.. autofunction::  generate_nucleotides
```

```{eval-rst}
.. autofunction::  update
```

```{eval-rst}
.. autofunction::  remove_vacant
```

```{eval-rst}
.. autofunction::  restore_vacant
```

```{eval-rst}
.. autofunction::  set_slim_state
```

## Summarizing tree sequences

Additionally, ``pyslim`` contains the following methods:

```{eval-rst}
.. autofunction::  individuals_alive_at
```

```{eval-rst}
.. autofunction::  individual_ages
```

```{eval-rst}
.. autofunction::  individual_ages_at
```

```{eval-rst}
.. autofunction::  individual_parents
```

```{eval-rst}
.. autofunction::  has_individual_parents
```

```{eval-rst}
.. autofunction::  population_size
```

## Utilities

```{eval-rst}
.. autofunction::  slim_time
```

```{eval-rst}
.. autofunction::  next_slim_mutation_id
```

```{eval-rst}
.. autofunction::  has_vacant_samples
```

```{eval-rst}
.. autofunction::  node_is_vacant
```


## Metadata

SLiM-specific metadata is made visible to the user by ``.metadata`` properties.
For instance:
```{code-cell}
ts.individual(4).metadata
```
shows that the fifth individual in the tree sequence was given pedigree ID ``495999`` by SLiM,
had parents with pedigree IDs ``493739`` and ``494784``,
was age 10 at the time that they died (or the simulation ended),
lived in subpopulation 1,
was female (because ``sex`` matches ``pyslim.INDIVIDUAL_TYPE_FEMALE``, below),
and has no additional metadata flags.


### Annotation

These two functions will add default SLiM metadata to a tree sequence (or the
underlying tables), which can then be modified and loaded into SLiM.

```{eval-rst}
.. autofunction:: pyslim.annotate
```

```{eval-rst}
.. autofunction:: pyslim.annotate_tables
```



(sec_constants_and_flags)=

## Constants and flags

```{eval-rst}
.. autodata:: NUCLEOTIDES
```

This is a flag used in `node.flags` (see {meth}`.remove_vacant`):

```{eval-rst}
.. autodata:: NODE_IS_VACANT_SAMPLE
```


These are the possible values for ``individual.metadata["sex"]``:

```{eval-rst}
.. autodata:: INDIVIDUAL_TYPE_HERMAPHRODITE

.. autodata:: INDIVIDUAL_TYPE_FEMALE

.. autodata:: INDIVIDUAL_TYPE_MALE
```

This is a flag used in ``individual.metadata["flags"]``
(*note*: this is *not* used in `individual.flags`, which is different!):

```{eval-rst}
.. autodata:: INDIVIDUAL_FLAG_MIGRATED
```

Finally, these are used in ``individual.flags``:

```{eval-rst}
.. autodata:: INDIVIDUAL_ALIVE

.. autodata:: INDIVIDUAL_REMEMBERED

.. autodata:: INDIVIDUAL_RETAINED
```

