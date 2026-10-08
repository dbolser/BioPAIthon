# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Package for working with phylogenetic trees.

See Also: http://biopython.org/wiki/Phylo

"""

__all__ = [
    "convert",
    "draw",
    "draw_ascii",
    "parse",
    "read",
    "to_igraph",
    "to_networkx",
    "write",
]

import importlib
import typing

from Bio.Phylo._io import convert
from Bio.Phylo._io import parse
from Bio.Phylo._io import read
from Bio.Phylo._io import write
from Bio.Phylo._utils import draw
from Bio.Phylo._utils import draw_ascii
from Bio.Phylo._utils import to_igraph
from Bio.Phylo._utils import to_networkx

# The tree format modules are imported on demand, so that "import Bio.Phylo"
# does not pay for every format and its dependencies (PhyloXML pulls in
# Bio.Align and hence NumPy, NexusIO pulls in Bio.Nexus, and CDAOIO needs
# rdflib).  They remain attributes of this package, as when Bio.Phylo._io
# imported them all: type checkers see real imports, and at run time the
# module __getattr__ below (PEP 562) imports each on first access.
if typing.TYPE_CHECKING:
    from Bio.Phylo import _cdao_owl
    from Bio.Phylo import CDAO
    from Bio.Phylo import CDAOIO
    from Bio.Phylo import NeXML
    from Bio.Phylo import NeXMLIO
    from Bio.Phylo import Newick
    from Bio.Phylo import NewickIO
    from Bio.Phylo import NexusIO
    from Bio.Phylo import PhyloXML
    from Bio.Phylo import PhyloXMLIO
else:
    _submodule_names = frozenset(
        [
            "_cdao_owl",
            "CDAO",
            "CDAOIO",  # needs rdflib
            "NeXML",
            "NeXMLIO",
            "Newick",
            "NewickIO",
            "NexusIO",
            "PhyloXML",
            "PhyloXMLIO",
        ]
    )

    def __getattr__(name):
        """Import tree format submodules of Bio.Phylo on first attribute access.

        A submodule that cannot be imported (CDAOIO without rdflib) is
        reported as a missing attribute, as when it was imported eagerly, so
        hasattr() returns False for it.
        """
        message = f"module {__name__!r} has no attribute {name!r}"
        if name in _submodule_names:
            try:
                return importlib.import_module(f"{__name__}.{name}")
            except ImportError as err:
                raise AttributeError(message) from err
        raise AttributeError(message)
