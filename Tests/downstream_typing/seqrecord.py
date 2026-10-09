# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.Seq and Bio.SeqRecord.

Type-checked with the mypy.ini beside it, never run.
"""

from typing import Any

import numpy as np
from typing_extensions import assert_type

from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord

record = SeqRecord(Seq("ACGT"), id="test")

assert_type(Seq("A"), Seq)
assert_type(record, SeqRecord)
assert_type(record[0], str)
assert_type(record[1:], SeqRecord)
# Any integer type indexes a record, as it does a Seq.
assert_type(record[np.int64(0)], str)

SeqRecord(Seq("A"), id=1)  # type: ignore[arg-type]

# The seq and id are rarely None (see SeqRecord.id), so they are used without
# narrowing first; wrong uses are still errors.
assert_type(record.seq, Seq | MutableSeq | Any)
record.seq.translate()
record.id.upper()
record.seq.frobnicate()  # type: ignore[union-attr]
n: int = record.id  # type: ignore[assignment]
# Hence extracting a feature from the record's seq, the commonest call, gives
# a Seq.
assert_type(SeqFeature(SimpleLocation(0, 2)).extract(record.seq), Seq)

# Records derived from a record are SeqRecords.
assert_type(record.reverse_complement(), SeqRecord)
assert_type(record.upper(), SeqRecord)

# Annotation values may be lists, as GenBank and EMBL store them.
record.annotations["taxonomy"].append("x")
SeqRecord(Seq("A"), annotations={"taxonomy": ["Bacteria"]})

assert_type(record.count("A"), int)
assert_type(record.isupper(), bool)
assert_type(record.islower(), bool)
