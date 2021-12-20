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

ts = pyslim.load("example_sim.trees")
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
  convert_alleles
  generate_nucleotides
  update_tables
  population_size
  default_slim_metadata
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
.. autofunction::  update_tables
```

## Summarizing tree sequences

Additionally, ``pyslim`` contains the following methods:

```{eval-rst}
.. autofunction::  population_size
```


## Additions to the tree sequence

The {class}`tskit.TreeSequence` class represents a sequence of trees
that describes the genealogical history of a population.
Here, we describe pyslim's additions to this class.


### The SLiM Tree Sequence


:::{eval-rst}
.. autoclass:: pyslim.SlimTreeSequence()
    :members:
:::


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
.. autofunction:: pyslim.annotate_defaults
:::

:::{eval-rst}
.. autofunction:: pyslim.annotate_defaults_tables
:::


### Provenances

:::{eval-rst}
.. autoclass:: pyslim.ProvenanceMetadata
   :members:
:::

:::{eval-rst}
.. autofunction:: pyslim.get_provenance
:::


(sec_constants_and_flags)=

## Constants and flags


:::{eval-rst}
.. data:: NUCLEOTIDES == ['A', 'C', 'G', 'T']

   Nucleotide states in nucleotide models are encoded as integers (0, 1, 2, 3),
   so a nucleotide encoded as ``k`` refers to nucleotide
   ``pyslim.NUCLEOTIDES[k]``.
:::

These flags are the possible values for ``node.metadata["genome_type"]``:

:::{eval-rst}
.. data:: GENOME_TYPE_AUTOSOME == 0

.. data:: GENOME_TYPE_X == 1

.. data:: GENOME_TYPE_Y == 2
:::


These flags are the possible values for ``individual.metadata["sex"]``:

:::{eval-rst}
.. data:: INDIVIDUAL_TYPE_HERMAPHRODITE == -1

.. data:: INDIVIDUAL_TYPE_FEMALE == 0

.. data:: INDIVIDUAL_TYPE_MALE == 1
:::

This is a flag used in ``individual.metadata["flags"]``:
:::{eval-rst}
.. data:: INDIVIDUAL_FLAG_MIGRATED == 0x01
:::

Finally, these are used in ``individual.flags``:

:::{eval-rst}
.. data:: INDIVIDUAL_ALIVE == 2**16

.. data:: INDIVIDUAL_REMEMBERED == 2**17

.. data:: INDIVIDUAL_RETAINED == 2**18
:::

