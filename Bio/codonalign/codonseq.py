# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com). All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with coding sequence.

CodonSeq class is inherited from Seq class. This is the core class to
deal with sequences in CodonAlignment in biopython.

"""

from Bio.Align.analysis import _lwl85
from Bio.Align.analysis import _ml
from Bio.Align.analysis import _ng86
from Bio.Align.analysis import _yn00
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class CodonSeq(Seq):
    """CodonSeq is designed to be within the SeqRecords of a CodonAlignment class.

    CodonSeq is useful as it allows the user to specify
    reading frame when translate CodonSeq

    CodonSeq also accepts codon style slice by calling
    get_codon() method.

    **Important:** Ungapped CodonSeq can be any length if you
    specify the rf_table. Gapped CodonSeq should be a
    multiple of three.

    >>> codonseq = CodonSeq("AAATTTGGGCCAAATTT", rf_table=(0,3,6,8,11,14))
    >>> print(codonseq.translate())
    KFGAKF

    test get_full_rf_table method

    >>> p = CodonSeq('AAATTTCCCGG-TGGGTTTAA', rf_table=(0, 3, 6, 9, 11, 14, 17))
    >>> full_rf_table = p.get_full_rf_table()
    >>> print(full_rf_table)
    [0, 3, 6, 9, 12, 15, 18]
    >>> print(p.translate(rf_table=full_rf_table, ungap_seq=False))
    KFPPWV*
    >>> p = CodonSeq('AAATTTCCCGGGAA-TTTTAA', rf_table=(0, 3, 6, 9, 14, 17))
    >>> print(p.get_full_rf_table())
    [0, 3, 6, 9, 12.0, 15, 18]
    >>> p = CodonSeq('AAA------------TAA', rf_table=(0, 3))
    >>> print(p.get_full_rf_table())
    [0, 3.0, 6.0, 9.0, 12.0, 15]

    """

    def __init__(self, data="", gap_char="-", rf_table=None):
        """Initialize the class."""
        # rf_table should be a tuple or list indicating the every
        # codon position along the sequence. For example:
        # sequence = 'AAATTTGGGCCAAATTT'
        # rf_table = (0, 3, 6, 8, 11, 14)
        # the translated protein sequences will be
        # AAA TTT GGG GCC AAA TTT
        #  K   F   G   A   K   F
        # Notice: rf_table applies to ungapped sequence. If there
        #   are gaps in the sequence, they will be discarded. This
        #   feature ensures the rf_table is independent of where the
        #   codon sequence appears in the alignment

        Seq.__init__(self, data.upper())
        self.gap_char = gap_char

        # check the length of the alignment to be a triple
        if rf_table is None:
            length = len(self)
            if length % 3 != 0:
                raise ValueError(
                    "Sequence length is not a multiple of "
                    "three (i.e. a whole number of codons)"
                )
            self.rf_table = list(range(0, length - self.count(gap_char), 3))
        else:
            # if gap_char in self:
            #    assert  len(self) % 3 == 0, \
            #            "Gapped sequence length is not a triple number"
            if not isinstance(rf_table, (tuple, list)):
                raise TypeError("rf_table should be a tuple or list object")
            if not all(isinstance(i, int) for i in rf_table):
                raise TypeError(
                    "Elements in rf_table should be int "
                    "that specify the codon positions of "
                    "the sequence"
                )
            self.rf_table = rf_table

    def get_codon(self, index):
        """Get the index codon from the sequence."""
        if len({i % 3 for i in self.rf_table}) != 1:
            raise RuntimeError(
                "frameshift detected. CodonSeq object is not able to deal with "
                "codon sequence with frameshift. Please use normal slice option."
            )
        if isinstance(index, int):
            if index != -1:
                return str(self[index * 3 : (index + 1) * 3])
            else:
                return str(self[index * 3 :])
        else:
            # This slice ensures that codon will always be the unit
            # in slicing (it won't change to other codon if you are
            # using reverse slicing such as [::-1]).
            # The idea of the code below is to first map the slice
            # to amino acid sequence and then transform it into
            # codon sequence.
            aa_index = range(len(self) // 3)

            def cslice(p):
                aa_slice = aa_index[p]
                codon_slice = ""
                for i in aa_slice:
                    codon_slice += self[i * 3 : i * 3 + 3]
                return str(codon_slice)

            codon_slice = cslice(index)
            return CodonSeq(codon_slice)

    def get_codon_num(self):
        """Return the number of codons in the CodonSeq."""
        return len(self.rf_table)

    def translate(
        self, codon_table=None, stop_symbol="*", rf_table=None, ungap_seq=True
    ):
        """Translate the CodonSeq based on the reading frame in rf_table.

        It is possible for the user to specify
        a rf_table at this point. If you want to include
        gaps in the translated sequence, this is the only
        way. ungap_seq should be set to true for this
        purpose.
        """
        if codon_table is None:
            codon_table = CodonTable.generic_by_id[1]
        amino_acids = []
        if ungap_seq:
            tr_seq = str(self).replace(self.gap_char, "")
        else:
            tr_seq = str(self)
        if rf_table is None:
            rf_table = self.rf_table
        p = -1  # initiation
        for i in rf_table:
            if isinstance(i, float):
                amino_acids.append("-")
                continue
            # elif '---' == tr_seq[i:i+3]:
            #    amino_acids.append('-')
            #    continue
            elif "-" in tr_seq[i : i + 3]:
                # considering two types of frameshift
                if p == -1 or p - i == 3:
                    p = i
                    codon = tr_seq[i : i + 6].replace("-", "")[:3]
                elif p - i > 3:
                    codon = tr_seq[i : i + 3]
                    p = i
            else:
                # normal condition without gaps
                codon = tr_seq[i : i + 3]
                p = i
            if codon in codon_table.stop_codons:
                amino_acids.append(stop_symbol)
                continue
            try:
                amino_acids.append(codon_table.forward_table[codon])
            except KeyError:
                raise RuntimeError(
                    f"Unknown codon detected ({codon}). Did you"
                    " forget to specify the ungap_seq argument?"
                )
        return "".join(amino_acids)

    def toSeq(self):
        """Convert DNA to seq object."""
        return Seq(str(self))

    def get_full_rf_table(self):
        """Return full rf_table of the CodonSeq records.

        A full rf_table is different from a normal rf_table in that
        it translate gaps in CodonSeq. It is helpful to construct
        alignment containing frameshift.
        """
        ungap_seq = str(self).replace("-", "")
        relative_pos = [self.rf_table[0]]
        for i in range(1, len(self.rf_table[1:]) + 1):
            relative_pos.append(self.rf_table[i] - self.rf_table[i - 1])
        full_rf_table = []
        codon_num = 0
        for i in range(0, len(self), 3):
            if self[i : i + 3] == self.gap_char * 3:
                full_rf_table.append(i + 0.0)
            elif relative_pos[codon_num] == 0:
                full_rf_table.append(i)
                codon_num += 1
            elif relative_pos[codon_num] in (-1, -2):
                # check the gap status of previous codon
                gap_stat = 3 - self.count("-", i - 3, i)
                if gap_stat == 3:
                    full_rf_table.append(i + relative_pos[codon_num])
                elif gap_stat == 2:
                    full_rf_table.append(i + 1 + relative_pos[codon_num])
                elif gap_stat == 1:
                    full_rf_table.append(i + 2 + relative_pos[codon_num])
                codon_num += 1
            elif relative_pos[codon_num] > 0:
                full_rf_table.append(i + 0.0)
            try:
                this_len = 3 - self.count("-", i, i + 3)
                relative_pos[codon_num] -= this_len
            except Exception:  # TODO: IndexError?
                # we probably reached the last codon
                pass
        return full_rf_table

    def full_translate(self, codon_table=None, stop_symbol="*"):
        """Apply full translation with gaps considered."""
        if codon_table is None:
            codon_table = CodonTable.generic_by_id[1]
        full_rf_table = self.get_full_rf_table()
        return self.translate(
            codon_table=codon_table,
            stop_symbol=stop_symbol,
            rf_table=full_rf_table,
            ungap_seq=False,
        )

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s)."""
        if len(gap) != 1 or not isinstance(gap, str):
            raise ValueError(f"Unexpected gap character, {gap!r}")
        return CodonSeq(str(self).replace(gap, ""), rf_table=self.rf_table)

    @classmethod
    def from_seq(cls, seq, rf_table=None):
        """Get codon sequence from sequence data."""
        if rf_table is None:
            return cls(str(seq))
        else:
            return cls(str(seq), rf_table=rf_table)


def _get_codon_list(codonseq):
    """List of codons according to full_rf_table for counting (PRIVATE)."""
    # if not isinstance(codonseq, CodonSeq):
    #    raise TypeError("_get_codon_list accept a CodonSeq object "
    #                    "({0} detected)".format(type(codonseq)))
    full_rf_table = codonseq.get_full_rf_table()
    codon_lst = []
    for i, k in enumerate(full_rf_table):
        if isinstance(k, int):
            start = k
            try:
                end = int(full_rf_table[i + 1])
            except IndexError:
                end = start + 3
            this_codon = str(codonseq[start:end])
            if len(this_codon) == 3:
                codon_lst.append(this_codon)
            else:
                codon_lst.append(str(this_codon.ungap()))
        elif str(codonseq[int(k) : int(k) + 3]) == "---":
            codon_lst.append("---")
        else:
            # this may be problematic, as normally no codon should
            # fall into this condition
            codon_lst.append(codonseq[int(k) : int(k) + 3])
    return codon_lst


def cal_dn_ds(codon_seq1, codon_seq2, method="NG86", codon_table=None, k=1, cfreq=None):
    """Calculate dN and dS of the given two sequences.

    Available methods:
        - NG86  - `Nei and Gojobori (1986)`_ (PMID 3444411).
        - LWL85 - `Li et al. (1985)`_ (PMID 3916709).
        - ML    - `Goldman and Yang (1994)`_ (PMID 7968486).
        - YN00  - `Yang and Nielsen (2000)`_ (PMID 10666704).

    .. _`Nei and Gojobori (1986)`: http://www.ncbi.nlm.nih.gov/pubmed/3444411
    .. _`Li et al. (1985)`: http://www.ncbi.nlm.nih.gov/pubmed/3916709
    .. _`Goldman and Yang (1994)`: http://mbe.oxfordjournals.org/content/11/5/725
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    Arguments:
     - codon_seq1 - CodonSeq or or SeqRecord that contains a CodonSeq
     - codon_seq2 - CodonSeq or or SeqRecord that contains a CodonSeq
     - w  - transition/transversion ratio
     - cfreq - Current codon frequency vector can only be specified
       when you are using ML method. Possible ways of
       getting cfreq are: F1x4, F3x4 and F61.

    """
    if isinstance(codon_seq1, CodonSeq) and isinstance(codon_seq2, CodonSeq):
        pass
    elif isinstance(codon_seq1, SeqRecord) and isinstance(codon_seq2, SeqRecord):
        codon_seq1 = codon_seq1.seq
        codon_seq2 = codon_seq2.seq
    else:
        raise TypeError(
            "cal_dn_ds accepts two CodonSeq objects or SeqRecord "
            "that contains CodonSeq as its seq!"
        )
    if len(codon_seq1.get_full_rf_table()) != len(codon_seq2.get_full_rf_table()):
        raise RuntimeError(
            f"full_rf_table length of seq1 ({len(codon_seq1.get_full_rf_table())})"
            f" and seq2 ({len(codon_seq2.get_full_rf_table())}) are not the same"
        )
    if cfreq is None:
        cfreq = "F3x4"
    elif cfreq is not None and method != "ML":
        raise RuntimeError("cfreq can only be specified when you are using ML method")
    if cfreq not in ("F1x4", "F3x4", "F61"):
        import warnings

        warnings.warn(
            f"Unknown cfreq ({cfreq}). "
            "Only F1x4, F3x4 and F61 are acceptable. Used F3x4 in the following."
        )
        cfreq = "F3x4"
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    seq1_codon_lst = _get_codon_list(codon_seq1)
    seq2_codon_lst = _get_codon_list(codon_seq2)
    # remove gaps in seq_codon_lst
    seq1 = []
    seq2 = []
    for i, j in zip(seq1_codon_lst, seq2_codon_lst):
        if ("-" not in i) and ("-" not in j):
            seq1.append(i)
            seq2.append(j)
    # The methods themselves are shared with Bio.Align.analysis.calculate_dn_ds
    if method == "ML":
        return _ml(seq1, seq2, cfreq, codon_table)
    elif method == "NG86":
        return _ng86(seq1, seq2, k, codon_table)
    elif method == "LWL85":
        return _lwl85(seq1, seq2, codon_table)
    elif method == "YN00":
        return _yn00(seq1, seq2, codon_table)
    # As when this looked the method up in a dict
    raise KeyError(method)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
