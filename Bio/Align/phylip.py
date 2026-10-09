# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the alignment format for input files for PHYLIP tools.

You are expected to use this module via the Bio.Align functions.

The "phylip" format cuts names at 10 characters.  The "phylip-relaxed"
format, which RAxML and PhyML read, allows names of any length that contain
no whitespace, and separates each name from its sequence by whitespace.
"""

from itertools import chain
from types import SimpleNamespace

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_PHYLIP_ID_WIDTH = 10


def _parse_header(line):
    """Return the number of sequences and their length from a header line (PRIVATE)."""
    words = line.split()
    if len(words) == 2:
        try:
            return int(words[0]), int(words[1])
        except ValueError:
            pass
    raise ValueError("Expected two integers in the first line, received '%s'" % line)


def _sanitize_name(name):
    """Remove the characters PHYLIP does not allow in a name (PRIVATE)."""
    name = name.strip()
    for char in "[](),":
        name = name.replace(char, "")
    for char in ":;":
        name = name.replace(char, "|")
    return name


class AlignmentWriter(interfaces.AlignmentWriter):
    """Clustalw alignment writer."""

    fmt = "PHYLIP"

    def _format_names(self, alignment):
        """Return the name of each sequence, padded to start its line (PRIVATE)."""
        names = []
        for record in alignment.sequences:
            try:
                name = record.id
            except AttributeError:
                name = ""
            else:
                name = _sanitize_name(name)[:_PHYLIP_ID_WIDTH]
            names.append(name.ljust(_PHYLIP_ID_WIDTH))
        return names

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Phylip format."""
        names = self._format_names(alignment)

        lines = []
        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")
        line = "%d %d\n" % (nseqs, length)
        lines.append(line)

        # From experimentation, the use of tabs is not understood by the
        # EMBOSS suite.  The nature of the expected white space is not
        # defined in the PHYLIP documentation, simply "These are in free
        # format, separated by blanks".  We'll use spaces to keep EMBOSS
        # happy.
        for name, sequence in zip(names, alignment):
            # Write the entire sequence to one line
            line = name + sequence + "\n"
            lines.append(line)
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Reads a Phylip alignment file and returns an Alignment iterator.

    Record names are limited to at most 10 characters.

    The parser determines from the file contents if the file format is
    sequential or interleaved, and parses the file accordingly.

    A file may hold several alignments one after another, each starting with
    its own header line, as PHYLIP's seqboot writes them. The number of
    columns in each header says where that alignment ends.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    """

    fmt = "PHYLIP"

    def _read_header(self, stream):
        line = stream.readline()
        if not line:
            raise ValueError("Empty file.") from None
        self._number_of_seqs, self._length_of_seqs = _parse_header(line)

    def _split_id(self, line):
        """Return the name at the start of a line and the residues after it (PRIVATE).

        The first 10 characters are the name, and the rest is sequence.
        """
        name = line[:_PHYLIP_ID_WIDTH].strip()
        seq = line[_PHYLIP_ID_WIDTH:].strip().replace(" ", "")
        return name, seq

    def _parse_interleaved_first_block(self, lines, seqs, names):
        for line in lines:
            name, seq = self._split_id(line)
            names.append(name)
            seqs.append([seq])

    def _parse_interleaved_other_blocks(self, stream, seqs):
        i = 0
        # Stop when the last row is complete, as another alignment may follow.
        length = len(seqs[-1][0])
        while length < self._length_of_seqs:
            line = stream.readline()
            if not line:
                break
            line = line.rstrip()
            if not line:
                if i != self._number_of_seqs:
                    raise ValueError(
                        f"Expected {self._number_of_seqs} sequence lines in each "
                        f"block of the interleaved alignment, found {i}"
                    )
                i = 0
            else:
                seq = line.replace(" ", "")
                seqs[i].append(seq)
                i += 1
                if i == self._number_of_seqs:
                    length += len(seq)
        if i != 0 and i != self._number_of_seqs:
            raise ValueError("Unexpected file format")

    def _parse_sequential(self, lines, seqs, names):
        length = 0
        for line in lines:
            if length == 0:
                name, seq = self._split_id(line)
                names.append(name)
                seqs.append([])
            else:
                seq = line.strip().replace(" ", "")
            seqs[-1].append(seq)
            length += len(seq)
            if length == self._length_of_seqs:
                if len(names) == self._number_of_seqs:
                    # The last sequence is complete; another alignment may follow.
                    break
                length = 0

    def _read_file(self, stream):
        names = []
        seqs = []
        lines = [stream.readline() for i in range(self._number_of_seqs)]
        if all(len(self._split_id(line)[1]) == self._length_of_seqs for line in lines):
            # One line per sequence; anything after it is another alignment.
            self._parse_interleaved_first_block(lines, seqs, names)
            return names, seqs
        line = stream.readline()
        if line.rstrip():
            # sequential file format
            lines.append(line)
            lines = chain(lines, iter(stream.readline, ""))
            self._parse_sequential(lines, seqs, names)
        else:
            # interleaved file format
            self._parse_interleaved_first_block(lines, seqs, names)
            self._parse_interleaved_other_blocks(stream, seqs)
        return names, seqs

    def _read_next_alignment(self, stream):
        try:
            self._number_of_seqs
        except AttributeError:
            # Not the first alignment; skip blank lines to the next header.
            line = stream.readline()
            while line.isspace():
                line = stream.readline()
            if not line:
                return
            self._number_of_seqs, self._length_of_seqs = _parse_header(line)
        names, seqs = self._read_file(stream)

        seqs = ["".join(seq) for seq in seqs]
        if len(seqs) != self._number_of_seqs:
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (len(seqs), self._number_of_seqs)
            )
        for seq in seqs:
            if len(seq) != self._length_of_seqs:
                raise ValueError(
                    "Expected all sequences to have length %d; found %d"
                    % (self._length_of_seqs, len(seq))
                )
            if "." in seq:
                raise ValueError("PHYLIP format no longer allows dots in sequence")

        seqs = [seq.encode() for seq in seqs]
        seqs, coordinates = Alignment.parse_printed_alignment(seqs)
        records = [
            SeqRecord(Seq(seq), id=name, description="")
            for (name, seq) in zip(names, seqs)
        ]
        alignment = Alignment(records, coordinates)
        del self._number_of_seqs
        del self._length_of_seqs
        return alignment


class RelaxedAlignmentWriter(AlignmentWriter):
    """Relaxed PHYLIP alignment writer.

    Names are written in full, padded with spaces to one more than the length
    of the longest name.  A name cannot contain whitespace, every sequence
    needs a name, and no two sequences can have the same name.  As for the
    "phylip" format, the characters ``[](),`` are removed from names, and
    ``:`` and ``;`` become ``|``.
    """

    def _format_names(self, alignment):
        """Return the name of each sequence, padded to start its line (PRIVATE)."""
        names = []
        for record in alignment.sequences:
            original = getattr(record, "id", None) or ""
            name = original.strip()
            if any(char.isspace() for char in name):
                raise ValueError(f"Whitespace not allowed in identifier: {name}")
            name = _sanitize_name(name)
            if not name:
                raise ValueError("Relaxed PHYLIP needs a name for every sequence")
            if name in names:
                raise ValueError(f"Repeated name {name!r} (originally {original!r})")
            names.append(name)
        width = max(map(len, names), default=0) + 1
        return [name.ljust(width) for name in names]


class RelaxedAlignmentIterator(AlignmentIterator):
    """Reads a relaxed PHYLIP alignment file and returns an Alignment iterator.

    Each name is the first word of its line, ending at the first whitespace,
    so it may have any length but cannot contain whitespace.  Everything after
    the name is sequence.  This is the relaxed PHYLIP that RAxML and PhyML read,
    and what Bio.AlignIO reads as "phylip-relaxed".

    The layout (sequential or interleaved), and where each alignment ends in
    a file holding several, are found as for the "phylip" format.
    """

    def _split_id(self, line):
        """Return the name at the start of a line and the residues after it (PRIVATE).

        The name ends at the first whitespace, and the rest is sequence.
        """
        name, *rest = line.split(None, 1) or [""]
        seq = rest[0].strip().replace(" ", "") if rest else ""
        return name, seq


# Bio.Align's registry maps the "phylip-relaxed" format to this object.
_relaxed = SimpleNamespace(
    AlignmentIterator=RelaxedAlignmentIterator,
    AlignmentWriter=RelaxedAlignmentWriter,
)
