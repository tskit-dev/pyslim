from __future__ import print_function
from __future__ import division

import numpy as np

#: Used in ``individual.flags`` to denote the individual is alive
#: when the tree sequence was written out.
INDIVIDUAL_ALIVE = np.uint32(2**16)

#: Used in ``individual.flags`` to denote the individual was
#: marked as "remembered".
INDIVIDUAL_REMEMBERED = np.uint32(2**17)

#: Used in ``individual.flags`` to denote the individual was
#: marked as "retained".
INDIVIDUAL_RETAINED = np.uint32(2**18)

# deprecated but keep it around for backwards compatibility
# (also, it means effectively the same thing as RETAINED)
INDIVIDUAL_FIRST_GEN = INDIVIDUAL_RETAINED

#: Mutation metadata records the nucleotide as an integer,
#: translated to ACGT by indexing this array,
#: so a nucleotide value of ``k`` actually means NUCLEOTIDES[k].
NUCLEOTIDES = ['A', 'C', 'G', 'T']

#: A value used in node metadata ("genome_type") to indicate the node is an autosome.
#: **DEPRECATED.**
GENOME_TYPE_AUTOSOME = 0

#: A value used in node metadata ("genome_type") to indicate the node is an X chromosome.
#: **DEPRECATED.**
GENOME_TYPE_X = 1

#: A value used in node metadata ("genome_type") to indicate the node is a Y chromosome.
#: **DEPRECATED.**
GENOME_TYPE_Y = 2

#: A value used in individual metadata ("sex") to indicate the individual is a hermaphrodite.
INDIVIDUAL_TYPE_HERMAPHRODITE = -1

#: A value used in individual metadata ("sex") to indicate the individual is a male.
INDIVIDUAL_TYPE_FEMALE = 0

#: A value used in individual metadata ("sex") to indicate the individual is a female.
INDIVIDUAL_TYPE_MALE = 1

#: An individual flag indicating the individual is a migrant.
INDIVIDUAL_FLAG_MIGRATED = np.uint32(1 << 1)

#: This flag exists because SLiM expects certain vacant nodes (=haplosomes)
#: to be marked as samples (those vacant nodes corresponding to alive individuals),
#: but including these as samples causes problems for certain operations in tskit.
#: So, if {meth}`.remove_vacant` is used to remove the 'sample' flags from those
#: vacant nodes, this flag is applied so that the sample flag can be easily put back
#: (by {meth}`.restore_vacant`). The flag does *not* mean simply that this is a
#: vacant sample node (indeed, if this flag is set then the `tskit.NODE_IS_SAMPLE`
#: flag is not expected to be set).
NODE_IS_VACANT_SAMPLE = np.uint32(1 << 16)

from pyslim.slim_metadata import *       # NOQA
from pyslim.slim_tree_sequence import *  # NOQA
from pyslim.provenance import *          # NOQA
from pyslim.methods import *             # NOQA
from pyslim.spatial import *             # NOQA
from pyslim._version import pyslim_version as __version__

