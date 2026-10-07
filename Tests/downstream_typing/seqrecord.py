# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.Seq and Bio.SeqRecord.

Type-checked with the mypy.ini beside it, never run.
"""

from typing_extensions import assert_type

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

record = SeqRecord(Seq("ACGT"), id="test")

assert_type(Seq("A"), Seq)
assert_type(record, SeqRecord)
assert_type(record[0], str)
assert_type(record[1:], SeqRecord)

SeqRecord(Seq("A"), id=1)  # type: ignore[arg-type]
