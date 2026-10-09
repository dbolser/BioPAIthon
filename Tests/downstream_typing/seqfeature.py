# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.SeqFeature.

Type-checked with the mypy.ini beside it, never run.
"""

from typing import Any

from typing_extensions import assert_type

from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import OneOfPosition
from Bio.SeqFeature import Position
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import UncertainPosition
from Bio.SeqFeature import UnknownPosition
from Bio.SeqFeature import WithinPosition
from Bio.SeqRecord import SeqRecord

location = SimpleLocation(1, 5)
f = SeqFeature(location, type="CDS")
record = SeqRecord(Seq("ACGT"))

# A coordinate is an int, or rarely an UnknownPosition, which is not one; it
# needs no narrowing for arithmetic, but a wrong use is still an error.
assert_type(location.start, int | Any)
assert_type(location.end, int | Any)
location.start + 1
s: str = location.start  # type: ignore[assignment]
if isinstance(location.start, BeforePosition):
    assert_type(location.start, BeforePosition)
assert_type((location + SimpleLocation(7, 9)).start, int | Any)

# The location is None only on a feature still being built, so it needs no
# narrowing either; a missing attribute is still an error.
f.location.start + 1
f.location.strand
f.location.frobnicate()  # type: ignore[union-attr]
if isinstance(f.location, CompoundLocation):
    assert_type(f.location.operator, str)
f.qualifiers["note"].append("x")

# The id stays None on a swiss feature with no FTId, so it does need a check.
SeqFeature(location, id=None)
assert_type(f.id, str | None)
f.id.upper()  # type: ignore[union-attr]

# extract gives a Seq for a Seq or MutableSeq, else the type it is given, and
# follows the references when they are given.
assert_type(f.extract("ACGT"), str)
assert_type(f.extract(Seq("ACGT")), Seq)
assert_type(f.extract(MutableSeq("ACGT")), Seq)
assert_type(f.extract(record), SeqRecord)
assert_type(f.extract(Seq("A"), references={"x": Seq("ACGT")}), Seq)
assert_type(f.extract("A", references={"x": "ACGT"}), str)
assert_type(f.extract(record, references={"x": record}), SeqRecord)
assert_type(location.extract(Seq("ACGT")), Seq)
assert_type((location + SimpleLocation(7, 9)).extract(record), SeqRecord)
f.extract(5)  # type: ignore[call-overload]

# Mixing kinds works at run time, but what comes back depends on which parts
# of the location have a ref, so the type is the union of all three.
assert_type(f.extract(Seq("A"), references={"x": record}), str | Seq | SeqRecord)

# translate needs a Seq, MutableSeq or SeqRecord; a str raises TypeError.
assert_type(f.translate(Seq("ATGTAA")), Seq)
assert_type(f.translate(MutableSeq("ATGTAA")), Seq)
assert_type(f.translate(record), SeqRecord)
f.translate("ATGTAA")  # type: ignore[call-overload]

# Adding locations joins them; adding an int shifts them.
assert_type(location + SimpleLocation(7, 9), CompoundLocation)
assert_type(location + 3, SimpleLocation)
assert_type(3 + location, SimpleLocation)
assert_type(location - 1, SimpleLocation)
assert_type(SimpleLocation(7, 9) + (location + SimpleLocation(1, 2)), CompoundLocation)

# Positions are ints, and shifting one keeps its class; UnknownPosition is not
# an int.
int(ExactPosition(5))
assert_type(ExactPosition(5) + 1, ExactPosition)
assert_type(UncertainPosition(5) + 1, UncertainPosition)
assert_type(WithinPosition(10, left=10, right=13) + 1, WithinPosition)
assert_type(BeforePosition(5) + 1, BeforePosition)
assert_type(AfterPosition(5) + 1, AfterPosition)
assert_type(OneOfPosition(5, [ExactPosition(5), ExactPosition(8)]) + 1, OneOfPosition)
assert_type(Position.fromstring("<5"), Position)
int(UnknownPosition())  # type: ignore[call-overload]
SimpleLocation(UnknownPosition(), 5)
SimpleLocation("1", 5)  # type: ignore[arg-type]

# Only None is accepted for the obsolete sub_features.
SeqFeature(location, sub_features=[f])  # type: ignore[arg-type]
