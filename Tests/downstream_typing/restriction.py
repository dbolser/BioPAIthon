# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.Restriction.

Type-checked with the mypy.ini beside it, never run.

The enzymes are only built when Bio.Restriction is imported, so a type
checker cannot see them; a ``__getattr__`` that only type checkers see gives
them the type Any. The other names the package exports keep their real types.
"""

from typing import Any

from typing_extensions import assert_type

from Bio.Restriction import AllEnzymes
from Bio.Restriction import Analysis
from Bio.Restriction import EcoRI
from Bio.Restriction import RestrictionBatch
from Bio.Restriction.Restriction import BamHI
from Bio.Seq import Seq

EcoRI.search(Seq("GAATTC"))

assert_type(EcoRI, Any)
assert_type(BamHI, Any)
assert_type(Analysis, type[Analysis])
assert_type(AllEnzymes, RestrictionBatch)
