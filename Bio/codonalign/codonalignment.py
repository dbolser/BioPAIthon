# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com). All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with Codon Alignment.

CodonAlignment class is inherited from MultipleSeqAlignment class. This is
the core class to deal with codon alignment in biopython.
"""

import warnings

from Bio import BiopythonWarning
from Bio.Align import MultipleSeqAlignment
from Bio.Align.analysis import _mktest
from Bio.codonalign.codonseq import _get_codon_list
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import CodonSeq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord


class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.

    >>> from Bio.SeqRecord import SeqRecord
    >>> a = SeqRecord(CodonSeq("AAAACGTCG"), id="Alpha")
    >>> b = SeqRecord(CodonSeq("AAA---TCG"), id="Beta")
    >>> c = SeqRecord(CodonSeq("AAAAGGTGG"), id="Gamma")
    >>> print(CodonAlignment([a, b, c]))
    CodonAlignment with 3 rows and 9 columns (3 codons)
    AAAACGTCG Alpha
    AAA---TCG Beta
    AAAAGGTGG Gamma

    """

    def __init__(self, records="", name=None):
        """Initialize the class."""
        MultipleSeqAlignment.__init__(self, records)

        # check the type of the alignment to be nucleotide
        for rec in self:
            if not isinstance(rec.seq, CodonSeq):
                raise TypeError(
                    "CodonSeq objects are expected in each SeqRecord in CodonAlignment"
                )

        if self.get_alignment_length() % 3 != 0:
            raise ValueError(
                "Alignment length is not a multiple of "
                "three (i.e. a whole number of codons)"
            )

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        This output is indicated to be readable, but large alignment
        is shown truncated. A maximum of 20 rows (sequences) and
        60 columns (20 codons) are shown, with the record identifiers.
        This should fit nicely on a single screen. e.g.

        """
        rows = len(self._records)
        lines = [
            "CodonAlignment with %i rows and %i columns (%i codons)"
            % (rows, self.get_alignment_length(), self.get_aln_length())
        ]

        if rows <= 60:
            lines.extend([self._str_line(rec, length=60) for rec in self._records])
        else:
            lines.extend([self._str_line(rec, length=60) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1], length=60))
        return "\n".join(lines)

    def __getitem__(self, index):
        """Return a CodonAlignment object for single indexing."""
        if isinstance(index, int):
            return self._records[index]
        elif isinstance(index, slice):
            return CodonAlignment(self._records[index])
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, int):
            return self._records[row_index][col_index]
        elif isinstance(col_index, int):
            return "".join(str(rec[col_index]) for rec in self._records[row_index])
        else:
            return MultipleSeqAlignment(
                rec[col_index] for rec in self._records[row_index]
            )

    def __add__(self, other):
        """Combine two codonalignments with the same number of rows by adding them.

        The method also allows to combine a CodonAlignment object with a
        MultipleSeqAlignment object. The following rules apply:

            * CodonAlignment + CodonAlignment -> CodonAlignment
            * CodonAlignment + MultipleSeqAlignment -> MultipleSeqAlignment
        """
        if isinstance(other, CodonAlignment):
            if len(self) != len(other):
                raise ValueError(
                    "When adding two alignments they must have the same length"
                    " (i.e. same number or rows)"
                )
            warnings.warn(
                "Please make sure the two CodonAlignment objects are sharing the same codon table. This is not checked by Biopython.",
                BiopythonWarning,
            )
            merged = (
                SeqRecord(seq=CodonSeq(left.seq + right.seq))
                for left, right in zip(self, other)
            )
            return CodonAlignment(merged)
        elif isinstance(other, MultipleSeqAlignment):
            if len(self) != len(other):
                raise ValueError(
                    "When adding two alignments they must have the same length"
                    " (i.e. same number or rows)"
                )
            return self.toMultipleSeqAlignment() + other
        else:
            raise TypeError(
                "Only CodonAlignment or MultipleSeqAlignment object can be"
                f" added with a CodonAlignment object. {object(other)} detected."
            )

    def get_aln_length(self):
        """Get alignment length."""
        return self.get_alignment_length() // 3

    def toMultipleSeqAlignment(self):
        """Convert the CodonAlignment to a MultipleSeqAlignment.

        Return a MultipleSeqAlignment containing all the
        SeqRecord in the CodonAlignment using Seq to store
        sequences
        """
        alignments = [SeqRecord(rec.seq.toSeq(), id=rec.id) for rec in self._records]
        return MultipleSeqAlignment(alignments)

    def get_dn_ds_matrix(self, method="NG86", codon_table=None):
        """Available methods include NG86, LWL85, YN00 and ML.

        Argument:
         - method       - Available methods include NG86, LWL85, YN00 and ML.
         - codon_table  - Codon table to use for forward translation.

        """
        from Bio.Phylo.TreeConstruction import DistanceMatrix as DM

        if codon_table is None:
            codon_table = CodonTable.generic_by_id[1]
        names = [i.id for i in self._records]
        size = len(self._records)
        dn_matrix = []
        ds_matrix = []
        for i in range(size):
            dn_matrix.append([])
            ds_matrix.append([])
            for j in range(i + 1):
                if i != j:
                    dn, ds = cal_dn_ds(
                        self._records[i],
                        self._records[j],
                        method=method,
                        codon_table=codon_table,
                    )
                    dn_matrix[i].append(dn)
                    ds_matrix[i].append(ds)
                else:
                    dn_matrix[i].append(0.0)
                    ds_matrix[i].append(0.0)
        dn_dm = DM(names, matrix=dn_matrix)
        ds_dm = DM(names, matrix=ds_matrix)
        return dn_dm, ds_dm

    def get_dn_ds_tree(
        self, dn_ds_method="NG86", tree_method="UPGMA", codon_table=None
    ):
        """Construct dn tree and ds tree.

        Argument:
         - dn_ds_method - Available methods include NG86, LWL85, YN00 and ML.
         - tree_method  - Available methods include UPGMA and NJ.

        """
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

        if codon_table is None:
            codon_table = CodonTable.generic_by_id[1]
        dn_dm, ds_dm = self.get_dn_ds_matrix(
            method=dn_ds_method, codon_table=codon_table
        )
        dn_constructor = DistanceTreeConstructor()
        ds_constructor = DistanceTreeConstructor()
        if tree_method == "UPGMA":
            dn_tree = dn_constructor.upgma(dn_dm)
            ds_tree = ds_constructor.upgma(ds_dm)
        elif tree_method == "NJ":
            dn_tree = dn_constructor.nj(dn_dm)
            ds_tree = ds_constructor.nj(ds_dm)
        else:
            raise RuntimeError(
                f"Unknown tree method ({tree_method}). Only NJ and UPGMA are accepted."
            )
        return dn_tree, ds_tree

    @classmethod
    def from_msa(cls, align):
        """Convert a MultipleSeqAlignment to CodonAlignment.

        Function to convert a MultipleSeqAlignment to CodonAlignment.
        It is the user's responsibility to ensure all the requirement
        needed by CodonAlignment is met.
        """
        rec = [SeqRecord(CodonSeq(str(i.seq)), id=i.id) for i in align._records]
        return cls(rec)


def mktest(codon_alns, codon_table=None, alpha=0.05):
    """McDonald-Kreitman test for neutrality.

    Implement the McDonald-Kreitman test for neutrality (PMID: 1904993)
    This method counts changes rather than sites
    (http://mkt.uab.es/mkt/help_mkt.asp).

    Arguments:
     - codon_alns  - list of CodonAlignment to compare (each
       CodonAlignment object corresponds to gene sampled from a species)

    Return the p-value of test result. Codon positions with a gap in any
    sequence are skipped, as in ``Bio.Align.analysis.mktest``.
    """
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    if not all(isinstance(i, CodonAlignment) for i in codon_alns):
        raise TypeError("mktest accepts CodonAlignment list.")
    codon_aln_len = [i.get_alignment_length() for i in codon_alns]
    if len(set(codon_aln_len)) != 1:
        raise RuntimeError(
            "CodonAlignment object for mktest should be of equal length."
        )
    codon_num = codon_aln_len[0] // 3
    # prepare codon_lst
    codon_lst = []
    for codon_aln in codon_alns:
        codon_lst.append([])
        for i in codon_aln:
            codon_lst[-1].append(_get_codon_list(i.seq))
    codon_set = []
    for i in range(codon_num):
        uniq_codons = []
        for j in codon_lst:
            uniq_codon = {k[i] for k in j}
            uniq_codons.append(uniq_codon)
        if not any("-" in codon for codons in uniq_codons for codon in codons):
            codon_set.append(uniq_codons)
    return _mktest(codon_set, codon_table)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
