# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the alignment format for input files for PHYLIP tools.

You are expected to use this module via the Bio.Align functions.
"""

from itertools import chain

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


class AlignmentWriter(interfaces.AlignmentWriter):
    """Clustalw alignment writer."""

    fmt = "PHYLIP"

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Phylip format."""
        names = []
        for record in alignment.sequences:
            try:
                name = record.id
            except AttributeError:
                name = ""
            else:
                name = name.strip()
                for char in "[](),":
                    name = name.replace(char, "")
                for char in ":;":
                    name = name.replace(char, "|")
                name = name[:_PHYLIP_ID_WIDTH]
            names.append(name)

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
            line = name[:_PHYLIP_ID_WIDTH].ljust(_PHYLIP_ID_WIDTH) + sequence + "\n"
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

    def _parse_interleaved_first_block(self, lines, seqs, names):
        for line in lines:
            line = line.rstrip()
            name = line[:_PHYLIP_ID_WIDTH].strip()
            seq = line[_PHYLIP_ID_WIDTH:].strip().replace(" ", "")
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
                line = line.rstrip()
                name = line[:_PHYLIP_ID_WIDTH].strip()
                seq = line[_PHYLIP_ID_WIDTH:].strip()
                names.append(name)
                seqs.append([])
            else:
                seq = line.strip()
            seq = seq.replace(" ", "")
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
        if all(
            len(line[_PHYLIP_ID_WIDTH:].replace(" ", "").strip())
            == self._length_of_seqs
            for line in lines
        ):
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
