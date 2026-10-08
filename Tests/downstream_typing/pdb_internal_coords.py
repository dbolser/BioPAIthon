# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.PDB.internal_coords.

Type-checked with the mypy.ini beside it, never run.
"""

from Bio.PDB import internal_coords
from Bio.PDB.internal_coords import MissingAtomError
from Bio.PDB.internal_coords import NoSuchName  # type: ignore[attr-defined]

# The deprecated class still type-checks; only the runtime warns.
try:
    pass
except (MissingAtomError, internal_coords.MissingAtomError):
    pass

# The module's runtime __getattr__ must not make other names acceptable.
typo = internal_coords.NoSuchName  # type: ignore[attr-defined]
