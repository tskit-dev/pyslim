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
  individual_ages_at
  slim_time
  convert_alleles
  generate_nucleotides
  population_size
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

## Summarizing tree sequences

Additionally, ``pyslim`` contains the following methods:

```{eval-rst}
.. autofunction::  individuals_alive_at
```

```{eval-rst}
.. autofunction::  individual_ages_at
```

```{eval-rst}
.. autofunction::  slim_time
```

```{eval-rst}
.. autofunction::  population_size
```


## Metadata

SLiM-specific metadata is made visible to the user by ``.metadata`` properties.
For instance:
```{code-cell}
ts.node(4).metadata
```
shows that the fifth node in the tree sequence was given pedigree ID ``982740`` by SLiM,
is *not* a null genome, and has ``genome_type`` zero, which corresponds to an autosome 
(see below).


### Annotation

These two functions will add default SLiM metadata to a tree sequence (or the
underlying tables), which can then be modified and loaded into SLiM.

:::{eval-rst}
.. autofunction:: pyslim.annotate
:::

:::{eval-rst}
.. autofunction:: pyslim.annotate_tables
:::



(sec_constants_and_flags)=

## Constants and flags


:::{eval-rst}
.. autodata:: NUCLEOTIDES
:::

These flags are the possible values for ``node.metadata["genome_type"]``:

:::{eval-rst}
.. autodata:: GENOME_TYPE_AUTOSOME

.. autodata:: GENOME_TYPE_X

.. autodata:: GENOME_TYPE_Y
:::


These flags are the possible values for ``individual.metadata["sex"]``:

:::{eval-rst}
.. autodata:: INDIVIDUAL_TYPE_HERMAPHRODITE

.. autodata:: INDIVIDUAL_TYPE_FEMALE

.. autodata:: INDIVIDUAL_TYPE_MALE
:::

This is a flag used in ``individual.metadata["flags"]``:
:::{eval-rst}
.. data:: INDIVIDUAL_FLAG_MIGRATED == 0x01
:::

Finally, these are used in ``individual.flags``:

:::{eval-rst}
.. autodata:: INDIVIDUAL_ALIVE

.. autodata:: INDIVIDUAL_REMEMBERED

.. autodata:: INDIVIDUAL_RETAINED
:::

