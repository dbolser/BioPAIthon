# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.Seq.

Type-checked with the mypy.ini beside it, never run.
"""

from typing_extensions import assert_type

from Bio.Seq import back_transcribe
from Bio.Seq import complement
from Bio.Seq import complement_rna
from Bio.Seq import MutableSeq
from Bio.Seq import reverse_complement
from Bio.Seq import reverse_complement_rna
from Bio.Seq import Seq
from Bio.Seq import transcribe
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord

seq = Seq("ACGT")
mutable_seq = MutableSeq("ACGT")
record = SeqRecord(Seq("ACGT"))

# Methods return the class they are called on.
assert_type(seq.reverse_complement(), Seq)
assert_type(mutable_seq.reverse_complement(), MutableSeq)
assert_type(seq.translate(), Seq)
assert_type(mutable_seq.translate(), MutableSeq)
assert_type(mutable_seq.upper(inplace=True), MutableSeq)
assert_type(mutable_seq.replace("A", "T"), MutableSeq)
assert_type(seq.strip(), Seq)

# Indexing gives a letter; slicing gives the same class.
assert_type(seq[0], str)
assert_type(mutable_seq[0], str)
assert_type(seq[1:], Seq)
assert_type(mutable_seq[1:], MutableSeq)

# split and rsplit always give Seq objects.
assert_type(MutableSeq("A-C").split("-"), list[Seq])
assert_type(seq.rsplit("G"), list[Seq])

# Concatenation keeps the class of the left operand.
assert_type(seq + mutable_seq, Seq)
assert_type(mutable_seq + seq, MutableSeq)
assert_type(seq + "A", Seq)
assert_type("A" + seq, Seq)
assert_type("A" + mutable_seq, MutableSeq)
assert_type(seq * 2, Seq)
assert_type(2 * mutable_seq, MutableSeq)

# transcribe, back_transcribe and translate give a Seq for any sequence.
assert_type(transcribe("ATG"), str)
assert_type(transcribe(seq), Seq)
assert_type(transcribe(mutable_seq), Seq)
assert_type(back_transcribe("AUG"), str)
assert_type(back_transcribe(seq), Seq)
assert_type(back_transcribe(mutable_seq), Seq)
assert_type(translate("ATG"), str)
assert_type(translate(seq), Seq)
assert_type(translate(mutable_seq), Seq)

# The complement functions keep the class they are given.
assert_type(complement("ATG"), str)
assert_type(complement(seq), Seq)
assert_type(complement(mutable_seq), MutableSeq)
assert_type(complement_rna("AUG"), str)
assert_type(complement_rna(seq), Seq)
assert_type(complement_rna(mutable_seq), MutableSeq)
assert_type(reverse_complement("ATG"), str)
assert_type(reverse_complement(seq), Seq)
assert_type(reverse_complement(mutable_seq), MutableSeq)
assert_type(reverse_complement_rna("AUG"), str)
assert_type(reverse_complement_rna(seq), Seq)
assert_type(reverse_complement_rna(mutable_seq), MutableSeq)

# Of the module functions, only reverse_complement takes a SeqRecord.
assert_type(reverse_complement(record), SeqRecord)
transcribe(record)  # type: ignore[call-overload]
translate(record)  # type: ignore[call-overload]
complement(record)  # type: ignore[call-overload]

# Wrong types are still errors.
x: int = Seq("A").translate()  # type: ignore[assignment]
seq[0] = "A"  # type: ignore[index]
mutable_seq[0] = seq  # type: ignore[call-overload]
seq.frobnicate()  # type: ignore[attr-defined]
