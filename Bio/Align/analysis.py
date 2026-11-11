# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com). All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for performing calculations on codon alignments."""

import sys
import warnings
from collections import Counter
from collections import defaultdict
from heapq import heapify
from heapq import heappop
from heapq import heappush
from itertools import permutations
from math import erfc
from math import floor
from math import isnan
from math import log
from math import nan
from math import sqrt

import numpy as np

from Bio import BiopythonWarning
from Bio.Align import Alignment
from Bio.Data import CodonTable


def calculate_dn_ds(alignment, method="NG86", codon_table=None, k=1, cfreq=None):
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
     - k  - transition/transversion rate ratio
     - cfreq - Current codon frequency vector can only be specified
       when you are using ML method. Possible ways of
       getting cfreq are: F1x4, F3x4 and F61.

    Returns a tuple ``(dN, dS)``. Where the sequences are too divergent for
    NG86 or LWL85 to correct for multiple substitutions, a distance is
    undefined, and is nan, with a ``BiopythonWarning``: for NG86, when the
    proportion of synonymous or of nonsynonymous differences is 3/4 or more;
    for LWL85, when a class of sites has too many differences for its
    number of sites, which can happen with few codons even for similar
    sequences. YN00 instead falls back from F84 to K80 to JC69, as the yn00
    program of PAML does, with a distance of at most 99, and ML maximises a
    likelihood, which needs no such correction.

    """
    if cfreq is None:
        cfreq = "F3x4"
    elif cfreq is not None and method != "ML":
        raise ValueError("cfreq can only be specified when you are using ML method")
    elif cfreq not in ("F1x4", "F3x4", "F61"):
        raise ValueError("cfreq must be 'F1x4', 'F3x4', or 'F61'")
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    codons1 = []
    codons2 = []
    sequence1, sequence2 = alignment.sequences
    try:
        sequence1 = sequence1.seq  # stupid SeqRecord
    except AttributeError:
        pass
    sequence1 = str(sequence1)
    try:
        sequence2 = sequence2.seq  # stupid SeqRecord
    except AttributeError:
        pass
    sequence2 = str(sequence2)
    aligned1, aligned2 = alignment.aligned
    for block1, block2 in zip(aligned1, aligned2):
        start1, end1 = block1
        start2, end2 = block2
        n1 = end1 - start1
        n2 = end2 - start2
        assert n1 == n2
        if n1 < 3:  # frameshift event, not a codon
            continue
        codons1.extend(sequence1[i : i + 3] for i in range(start1, end1, 3))
        codons2.extend(sequence2[i : i + 3] for i in range(start2, end2, 3))
    bases = {"A", "T", "C", "G"}
    for codon1 in codons1:
        if not all(nucleotide in bases for nucleotide in codon1):
            raise ValueError(
                f"Unrecognized character in {codon1} in the target sequence"
                " (Codons consist of A, T, C or G)"
            )
    for codon2 in codons2:
        if not all(nucleotide in bases for nucleotide in codon2):
            raise ValueError(
                f"Unrecognized character in {codon2} in the query sequence"
                " (Codons consist of A, T, C or G)"
            )
    if method == "ML":
        return _ml(codons1, codons2, cfreq, codon_table)
    elif method == "NG86":
        return _ng86(codons1, codons2, k, codon_table)
    elif method == "LWL85":
        return _lwl85(codons1, codons2, codon_table)
    elif method == "YN00":
        return _yn00(codons1, codons2, codon_table)
    else:
        raise ValueError(f"Unknown method '{method}'")


#################################################################
#              private functions for NG86 method
#################################################################


def _ng86(codons1, codons2, k, codon_table):
    """NG86 method main function (PRIVATE)."""
    S_sites1, N_sites1 = _count_site_NG86(codons1, codon_table=codon_table, k=k)
    S_sites2, N_sites2 = _count_site_NG86(codons2, codon_table=codon_table, k=k)
    S_sites = (S_sites1 + S_sites2) / 2.0
    N_sites = (N_sites1 + N_sites2) / 2.0
    SN = [0, 0]
    for codon1, codon2 in zip(codons1, codons2):
        SN = [
            m + n
            for m, n in zip(
                SN, _count_diff_NG86(codon1, codon2, codon_table=codon_table)
            )
        ]

    ps = SN[0] / S_sites
    pn = SN[1] / N_sites
    # The Jukes-Cantor correction needs a proportion of differences below
    # 3/4; at or above it the distance is undefined, and is nan. This was
    # -1, as in DistanceMatNG86() of PAML, which a ratio dN/dS took as 1.
    saturated = []
    if ps < 3 / 4:
        dS = abs(-3.0 / 4 * log(1 - 4.0 / 3 * ps))
    else:
        dS = nan
        saturated.append(f"synonymous differences, pS = {ps:g}, gives nan for dS")
    if pn < 3 / 4:
        dN = abs(-3.0 / 4 * log(1 - 4.0 / 3 * pn))
    else:
        dN = nan
        saturated.append(f"nonsynonymous differences, pN = {pn:g}, gives nan for dN")
    if saturated:
        warnings.warn(
            "NG86 cannot correct a proportion of differences of 3/4 or more "
            "for multiple substitutions: the proportion of "
            + ", and of ".join(saturated),
            BiopythonWarning,
        )
    return dN, dS


def _count_site_NG86(codons, codon_table, k=1):
    """Count synonymous and non-synonymous sites of a list of codons (PRIVATE).

    Arguments:
     - codons - A list of three letter codons.
     - k - transition/transversion rate ratio.

    """
    S_site = 0  # synonymous sites
    N_site = 0  # non-synonymous sites
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    bases = ("A", "T", "C", "G")
    for codon in codons:
        neighbor_codon = {"transition": [], "transversion": []}
        # classify neighbor codons
        codon = codon.replace("U", "T")
        for i, nucleotide in enumerate(codon):
            for base in bases:
                if nucleotide == base:
                    pass
                elif nucleotide in purine and base in purine:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transition"].append(this_codon)
                elif nucleotide in pyrimidine and base in pyrimidine:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transition"].append(this_codon)
                else:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transversion"].append(this_codon)
        # count synonymous and non-synonymous sites
        aa = codon_table.forward_table[codon]
        this_codon_N_site = this_codon_S_site = 0
        for neighbor in neighbor_codon["transition"]:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += k
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += k
            else:
                this_codon_N_site += k
        for neighbor in neighbor_codon["transversion"]:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += 1
            else:
                this_codon_N_site += 1
        norm_const = (this_codon_N_site + this_codon_S_site) / 3
        S_site += this_codon_S_site / norm_const
        N_site += this_codon_N_site / norm_const
    return (S_site, N_site)


def _count_diff_NG86(codon1, codon2, codon_table):
    """Count differences between two codons, three-letter string (PRIVATE).

    The function will take multiple pathways from codon1 to codon2
    into account.
    """
    # calculate_dn_ds checks its codons before getting here, but cal_dn_ds in
    # Bio.codonalign relies on these checks; without them a U would count as
    # a substitution from T.
    if not isinstance(codon1, str) or not isinstance(codon2, str):
        raise TypeError(
            "_count_diff_NG86 accepts string object to represent codon"
            f" ({type(codon1)}, {type(codon2)} detected)"
        )
    if len(codon1) != 3 or len(codon2) != 3:
        raise RuntimeError(
            "codon should be three letter string"
            f" ({len(codon1)}, {len(codon2)} detected)"
        )
    bases = ("A", "C", "G", "T")
    if not all(nucleotide in bases for nucleotide in codon1):
        raise RuntimeError(
            f"Unrecognized character detected in codon1 {codon1}"
            " (Codons consist of A, T, C or G)"
        )
    if not all(nucleotide in bases for nucleotide in codon2):
        raise RuntimeError(
            f"Unrecognized character detected in codon2 {codon2}"
            " (Codons consist of A, T, C or G)"
        )
    SN = [0, 0]  # synonymous and nonsynonymous counts
    if codon1 == codon2:
        return SN
    else:
        diff_pos = [
            i
            for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
            if nucleotide1 != nucleotide2
        ]

        def compare_codon(codon1, codon2, codon_table, weight=1):
            """Compare two codon accounting for different pathways."""
            sd = nd = 0
            if len(set(map(codon_table.forward_table.get, [codon1, codon2]))) == 1:
                sd += weight
            else:
                nd += weight
            return (sd, nd)

        if len(diff_pos) == 1:
            SN = [
                i + j
                for i, j in zip(
                    SN, compare_codon(codon1, codon2, codon_table=codon_table)
                )
            ]
        elif len(diff_pos) == 2:
            for i in diff_pos:
                temp_codon = codon1[:i] + codon2[i] + codon1[i + 1 :]
                SN = [
                    i + j
                    for i, j in zip(
                        SN,
                        compare_codon(
                            codon1, temp_codon, codon_table=codon_table, weight=0.5
                        ),
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN,
                        compare_codon(
                            temp_codon, codon2, codon_table=codon_table, weight=0.5
                        ),
                    )
                ]
        elif len(diff_pos) == 3:
            paths = list(permutations([0, 1, 2], 3))
            tmp_codon = []
            for index1, index2, index3 in paths:
                tmp1 = codon1[:index1] + codon2[index1] + codon1[index1 + 1 :]
                tmp2 = tmp1[:index2] + codon2[index2] + tmp1[index2 + 1 :]
                tmp_codon.append((tmp1, tmp2))
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(codon1, tmp1, codon_table, weight=0.5 / 3)
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(tmp1, tmp2, codon_table, weight=0.5 / 3)
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(tmp2, codon2, codon_table, weight=0.5 / 3)
                    )
                ]
    return SN


#################################################################
#               private functions for LWL85 method
#################################################################


def _lwl85(codons1, codons2, codon_table):
    """LWL85 method main function (PRIVATE).

    Nomenclature is according to Li et al. (1985), PMID 3916709.
    """
    codon_fold_dict = _get_codon_fold(codon_table)
    # count number of sites in different degenerate classes
    fold0 = [0, 0]
    fold2 = [0, 0]
    fold4 = [0, 0]
    for codon in codons1 + codons2:
        fold_num = codon_fold_dict[codon]
        for f in fold_num:
            if f == "0":
                fold0[0] += 1
            elif f == "2":
                fold2[0] += 1
            elif f == "4":
                fold4[0] += 1
    L = [sum(fold0) / 2.0, sum(fold2) / 2.0, sum(fold4) / 2.0]
    # count number of differences in different degenerate classes
    PQ = [0] * 6  # with P0, P2, P4, Q0, Q2, Q4 in each position
    for codon1, codon2 in zip(codons1, codons2):
        if codon1 == codon2:
            continue
        PQ = [
            i + j
            for i, j in zip(PQ, _diff_codon(codon1, codon2, fold_dict=codon_fold_dict))
        ]
    PQ = [i / j for i, j in zip(PQ, L * 2)]
    P = PQ[:3]
    Q = PQ[3:]
    # Kimura's correction for each class of sites needs 1 - 2P - Q and
    # 1 - 2Q to be positive. Where one is not, the class has too many
    # differences for its number of sites, A or both A and B are undefined,
    # and so is dS or dN, whichever uses them. They are nan, as in yn00.c of
    # PAML, which takes the logarithm without a check.
    A = []
    B = []
    undefined = []
    for name, i, j in zip(
        ("nondegenerate", "twofold degenerate", "fourfold degenerate"), P, Q
    ):
        if 1 - 2 * j > 0:
            B.append((1.0 / 2) * log(1.0 / (1 - 2 * j)))
        else:
            B.append(nan)
        if 1 - 2 * i - j > 0 and 1 - 2 * j > 0:
            A.append(
                (1.0 / 2) * log(1.0 / (1 - 2 * i - j))
                - (1.0 / 4) * log(1.0 / (1 - 2 * j))
            )
        else:
            A.append(nan)
            undefined.append(f"{name} sites (P = {i:g}, Q = {j:g})")
    dS = 3 * (L[1] * A[1] + L[2] * (A[2] + B[2])) / (L[1] + 3 * L[2])
    dN = 3 * (L[1] * B[1] + L[0] * (A[0] + B[0])) / (2 * L[1] + 3 * L[0])
    if undefined:
        nans = " and ".join(name for name, d in (("dN", dN), ("dS", dS)) if isnan(d))
        warnings.warn(
            "LWL85 cannot correct for multiple substitutions at the "
            + " or the ".join(undefined)
            + ", as 1 - 2P - Q or 1 - 2Q is not positive: there are too many "
            f"differences for the number of sites; returning nan for {nans}",
            BiopythonWarning,
        )
    return dN, dS


def _get_codon_fold(codon_table):
    """Classify different position in a codon into different folds (PRIVATE)."""
    fold_table = {}
    forward_table = codon_table.forward_table
    bases = {"A", "T", "C", "G"}
    for codon in forward_table:
        if "U" in codon:
            continue
        fold = ""
        codon_base_lst = list(codon)
        for i, base in enumerate(codon_base_lst):
            other_bases = bases - set(base)
            aa = []
            for other_base in other_bases:
                codon_base_lst[i] = other_base
                try:
                    aa.append(forward_table["".join(codon_base_lst)])
                except KeyError:
                    aa.append("stop")
            if aa.count(forward_table[codon]) == 0:
                fold += "0"
            elif aa.count(forward_table[codon]) in (1, 2):
                fold += "2"
            elif aa.count(forward_table[codon]) == 3:
                fold += "4"
            else:
                raise RuntimeError(
                    "Unknown Error, cannot assign the position to a fold"
                )
            codon_base_lst[i] = base
        fold_table[codon] = fold
    return fold_table


def _diff_codon(codon1, codon2, fold_dict):
    """Count number of different substitution types between two codons (PRIVATE).

    returns tuple (P0, P2, P4, Q0, Q2, Q4)

    Nomenclature is according to Li et al. (1958), PMID 3916709.
    """
    P0 = P2 = P4 = Q0 = Q2 = Q4 = 0
    fold_num = fold_dict[codon1]
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    for n, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2)):
        if nucleotide1 == nucleotide2:
            pass
        elif nucleotide1 in purine and nucleotide2 in purine:
            if fold_num[n] == "0":
                P0 += 1
            elif fold_num[n] == "2":
                P2 += 1
            elif fold_num[n] == "4":
                P4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
        elif nucleotide1 in pyrimidine and nucleotide2 in pyrimidine:
            if fold_num[n] == "0":
                P0 += 1
            elif fold_num[n] == "2":
                P2 += 1
            elif fold_num[n] == "4":
                P4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
        else:
            # nucleotide1 in purine and nucleotide2 in pyrimidine, or
            # nucleotide1 in pyrimidine and nucleotide2 in purine
            if fold_num[n] == "0":
                Q0 += 1
            elif fold_num[n] == "2":
                Q2 += 1
            elif fold_num[n] == "4":
                Q4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
    return (P0, P2, P4, Q0, Q2, Q4)


#################################################################
#               private functions for YN00 method
#################################################################


def _yn00(codons1, codons2, codon_table):
    """YN00 method main function (PRIVATE).

    Nomenclature is according to Yang and Nielsen (2000), PMID 10666704.
    Returns a tuple ``(dN, dS)``. If the estimates have not converged after 20
    iterations, the last ones are returned with a ``BiopythonWarning``. If the
    sequences have no synonymous or no nonsynonymous sites, both are nan, also
    with a ``BiopythonWarning``. A codon that the table gives as both a stop
    codon and a sense codon (TGA in table 27, for example) is taken to be the
    sense codon throughout, as ``Seq.translate()`` takes it.
    """
    from scipy.linalg import expm

    stop_codons = [
        codon
        for codon in codon_table.stop_codons
        if codon not in codon_table.forward_table
    ]
    if len(stop_codons) < len(codon_table.stop_codons):
        # Otherwise such a codon, which the sequences use as a sense codon,
        # would be a stop codon in the substitution matrix, with no path out
        codon_table = CodonTable.CodonTable(
            forward_table=codon_table.forward_table, stop_codons=stop_codons
        )
    fcodon = [
        {"A": 0, "G": 0, "C": 0, "T": 0},
        {"A": 0, "G": 0, "C": 0, "T": 0},
        {"A": 0, "G": 0, "C": 0, "T": 0},
    ]
    kappa = _get_kappa_YN00(codons1, codons2, codon_table)
    for codon in codons1 + codons2:
        # count sites at different codon position
        fcodon[0][codon[0]] += 1
        fcodon[1][codon[1]] += 1
        fcodon[2][codon[2]] += 1
    for i in range(3):
        tot = sum(fcodon[i].values())
        fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
    # Codon frequencies from the base frequencies at each codon position
    # (F3x4), as GetFreqs() in yn00.c of PAML computes them, scaled to sum to
    # one over the sense codons; stop codons have frequency zero.
    pi = {}
    for codon in codon_table.forward_table:
        if "U" not in codon:
            pi[codon] = fcodon[0][codon[0]] * fcodon[1][codon[1]] * fcodon[2][codon[2]]
    total = sum(pi.values())
    pi = {codon: freq / total for codon, freq in pi.items()}
    for codon in codon_table.stop_codons:
        if "U" not in codon:
            pi[codon] = 0
    # count synonymous sites and non-synonymous sites
    S_sites1, N_sites1, bfreqSN1 = _count_site_YN00(
        codons1, codons2, pi, k=kappa, codon_table=codon_table
    )
    S_sites2, N_sites2, bfreqSN2 = _count_site_YN00(
        codons2, codons1, pi, k=kappa, codon_table=codon_table
    )
    N_sites = (N_sites1 + N_sites2) / 2
    S_sites = (S_sites1 + S_sites2) / 2
    if not (S_sites1 + N_sites1 and S_sites2 + N_sites2 and S_sites and N_sites):
        # Without both kinds of sites omega = dN/dS is undefined, and the
        # weighting of the paths needs it; yn00.c gives nan here too
        warnings.warn(
            f"YN00 needs synonymous and nonsynonymous sites, and found "
            f"{S_sites:g} and {N_sites:g}; returning nan for dN and dS",
            BiopythonWarning,
        )
        return nan, nan
    bfreqSN = [{"A": 0, "T": 0, "C": 0, "G": 0}, {"A": 0, "T": 0, "C": 0, "G": 0}]
    for i in range(2):
        # Average over the sequences that have sites of this kind; yn00.c
        # averages over both, and so gives nan when one has none
        freqs = [f[i] for f in (bfreqSN1, bfreqSN2) if any(f[i].values())]
        for base in ("A", "T", "C", "G"):
            bfreqSN[i][base] = sum(f[base] for f in freqs) / len(freqs)
    # Start from the values yn00.c starts its first pair of sequences from
    t = 0.4
    w = 1
    tolerance = 1e-5
    dSdN_pre = [0, 0]
    for temp in range(20):
        # count synonymous and nonsynonymous differences under kappa, w, t
        codons = [
            codon
            for codon in list(codon_table.forward_table.keys())
            + codon_table.stop_codons
            if "U" not in codon
        ]
        Q = _get_Q(pi, kappa, w, codons, codon_table)
        P = expm(Q * t)
        TV = [0, 0, 0, 0]  # synonymous/nonsynonymous transition/transversion
        codon_npath = Counter(zip(codons1, codons2))
        for (nucleotide1, nucleotide2), count in codon_npath.items():
            tv = _count_diff_YN00(nucleotide1, nucleotide2, P, codons, codon_table)
            TV = [m + n * count for m, n in zip(TV, tv)]
        # according to the DistanceF84() function of yn00.c in paml,
        # the t (e.q. 10) appears in PMID: 10666704 is dS and dN
        _, dS = _distance_f84(S_sites, TV[0] / S_sites, TV[1] / S_sites, bfreqSN[0])
        _, dN = _distance_f84(N_sites, TV[2] / N_sites, TV[3] / N_sites, bfreqSN[1])
        dSdN = [dS, dN]
        t = dS * 3 * S_sites / (S_sites + N_sites) + dN * 3 * N_sites / (
            S_sites + N_sites
        )
        # As DistanceYN00() in yn00.c: w is at least 1e-5, and 99 if dS is 0
        w = 99 if dS < 1e-9 else max(1e-5, dN / dS)
        if all(abs(i - j) < tolerance for i, j in zip(dSdN, dSdN_pre)):
            return float(dN), float(dS)
        dSdN_pre = dSdN
    warnings.warn(
        f"YN00 did not converge in 20 iterations; returning the last "
        f"estimates (dN = {dN:g}, dS = {dS:g})",
        BiopythonWarning,
    )
    return float(dN), float(dS)


def _get_kappa_YN00(codons1, codons2, codon_table):
    """Estimate the transition to transversion rate ratio kappa (PRIVATE).

    A port of ``GetKappa()`` in ``yn00.c`` of PAML for two sequences. The
    positions that are non-degenerate in both codons give one estimate, and
    the third positions that are fourfold degenerate in both codons, which
    code for the same amino acid, give another. Each comes from the F84
    distance between the bases at those positions, with the frequencies of
    those bases. The two are averaged, weighted by their numbers of sites,
    leaving out an estimate that is not positive. If neither is, kappa is
    4.6, the value yn00.c starts from.
    """
    fold = _get_codon_fold(codon_table)
    amino_acids = codon_table.forward_table
    kappa = weight = 0
    index = {base: i for i, base in enumerate("TCAG")}
    for fold_class, positions in (("0", (0, 1, 2)), ("4", (2,))):
        # F[i, j] is the proportion of sites with bases i and j in the two
        # sequences, counted both ways round
        F = np.zeros((4, 4))
        for codon1, codon2 in zip(codons1, codons2):
            for i in positions:
                if fold[codon1][i] == fold[codon2][i] == fold_class and (
                    fold_class == "0" or amino_acids[codon1] == amino_acids[codon2]
                ):
                    F[index[codon1[i]], index[codon2[i]]] += 0.5
                    F[index[codon2[i]], index[codon1[i]]] += 0.5
        n = F.sum()
        if n == 0:
            continue
        F /= n
        # Computed as yn00.c does: with few sites, 1 - 2P - Q is often exactly
        # 0, where K80 fails, and computing P and Q another way can round it
        # to just above 0, which gives a kappa in the tens
        P = (F[0, 1] + F[2, 3]) * 2
        Q = 1 - np.trace(F) - P
        pi = dict(zip("TCAG", F.sum(axis=1)))
        kappa_class, _ = _distance_f84(n, P, Q, pi)
        if kappa_class > 0:
            kappa += kappa_class * n
            weight += n
    if weight == 0:
        return 4.6
    return kappa / weight


def _distance_f84(n, P, Q, pi):
    """Return kappa and the distance between two sequences under F84 (PRIVATE).

    A port of ``DistanceF84()`` in ``yn00.c`` of PAML, without the standard
    error. Where F84 fails it falls back to K80, and where K80 fails too, to
    JC69, as PAML does. Formulae and names are according to Yang and Nielsen
    (2000), PMID 10666704.

    Arguments:
     - n  - number of sites
     - P  - proportion of sites with a transitional difference
     - Q  - proportion of sites with a transversional difference
     - pi - base frequencies, a dict with keys "T", "C", "A" and "G"

    Returns a tuple ``(kappa, t)``: ``kappa`` is the HKY85 transition to
    transversion rate ratio, at most 999, or -1 if neither F84 nor K80 gives
    one; ``t`` is the distance, at most 99.
    """
    max_kappa = 999
    max_t = 99.0
    if P + Q > 1:
        return 1, max_t
    Y = pi["T"] + pi["C"]
    R = pi["A"] + pi["G"]
    tc = pi["T"] * pi["C"]
    ag = pi["A"] * pi["G"]
    kappa = -1
    fail_f84 = fail_k80 = False
    if Q < min(1e-10, 0.1 / n):
        fail_f84 = fail_k80 = True
    elif Y <= 0 or R <= 0 or (tc <= 0 and ag <= 0):
        fail_f84 = True
    else:
        A = tc / Y + ag / R
        B = tc + ag
        C = Y * R
        a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A)
        b = 1 - Q / (2 * C)
        if a <= 0 or b <= 0:
            fail_f84 = True
        else:
            a = -0.5 * log(a)  # this seems to be an error in YANG's original paper
            b = -0.5 * log(b)
            if b <= 0:
                fail_f84 = True
            else:
                kappa_f84 = a / b - 1
                t = 4 * b * (tc * (1 + kappa_f84 / Y) + ag * (1 + kappa_f84 / R) + C)
                kappa = (B + A * kappa_f84) / B
    if fail_f84 and not fail_k80:
        a = 1 - 2 * P - Q
        b = 1 - 2 * Q
        if a <= 0 or b <= 0:
            fail_k80 = True
        else:
            a = -log(a)
            b = -log(b)
            if b <= 0:
                fail_k80 = True
            else:
                kappa = (0.5 * a - 0.25 * b) / (0.25 * b)
                t = 0.5 * a + 0.25 * b
    if fail_k80:
        P += Q
        if P < 0.75:
            t = abs(-0.75 * log(1 - P * 4 / 3))  # abs() turns -0.0 into 0.0
        elif n > 1:
            # JC69 fails too. As yn00.c, take P just below 3/4, at
            # 3/4 (n - 1) / n, so that t = 3/4 log(n).
            t = 0.75 * log(n)
        else:
            # That would not be positive, so use the maximum, as the comment
            # on DistanceF84() in yn00.c says it does when JC69 fails.
            t = max_t
        t = min(t, max_t)
    return min(kappa, max_kappa), t


def _count_site_YN00(codons1, codons2, pi, k, codon_table):
    """Site counting method from Ina / Yang and Nielsen (PRIVATE).

    Method from `Ina (1995)`_ as modified by `Yang and Nielsen (2000)`_.
    This will return the total number of synonymous and nonsynonymous sites
    and base frequencies in each category. The function is equivalent to
    the ``CountSites()`` function in ``yn00.c`` of PAML.

    .. _`Ina (1995)`: https://doi.org/10.1007/BF00167113
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    """
    length = len(codons1)
    if length != len(codons2):
        raise ValueError(
            f"Expected the same number of aligned codons in both sequences, "
            f"found {length} and {len(codons2)}"
        )
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    bases = ("A", "T", "C", "G")
    codon_dict = codon_table.forward_table
    stop = codon_table.stop_codons
    codon_npath = Counter(zip(codons1, codons2))
    S_sites = N_sites = 0
    freqSN = [
        {"A": 0, "T": 0, "C": 0, "G": 0},  # synonymous
        {"A": 0, "T": 0, "C": 0, "G": 0},
    ]  # nonsynonymous
    for codon_pair, npath in codon_npath.items():
        codon = codon_pair[0]
        S = N = 0
        for pos in range(3):
            for base in bases:
                if codon[pos] == base:
                    continue
                neighbor_codon = codon[:pos] + base + codon[pos + 1 :]
                if neighbor_codon in stop:
                    continue
                weight = pi[neighbor_codon]
                if codon[pos] in pyrimidine and base in pyrimidine:
                    weight *= k
                elif codon[pos] in purine and base in purine:
                    weight *= k
                # As CountSites() in yn00.c, count the base frequencies with
                # this codon's own base, not the neighbouring codon's
                if codon_dict[codon] == codon_dict[neighbor_codon]:
                    S += weight
                    freqSN[0][codon[pos]] += weight * npath
                else:
                    N += weight
                    freqSN[1][codon[pos]] += weight * npath
        S_sites += S * npath
        N_sites += N * npath
    if S_sites + N_sites == 0:
        # Every neighbouring codon is a stop codon or has frequency zero
        return S_sites, N_sites, freqSN
    norm_const = 3 * length / (S_sites + N_sites)
    S_sites *= norm_const
    N_sites *= norm_const
    for i in freqSN:
        norm_const = sum(i.values())
        if norm_const:  # zero if the sequence has no sites of this kind
            for b in i:
                i[b] /= norm_const
    return S_sites, N_sites, freqSN


def _count_diff_YN00(codon1, codon2, P, codons, codon_table):
    """Count the differences between two codons (three-letter strings; PRIVATE).

    Returns a list of the numbers of synonymous transitions, synonymous
    transversions, nonsynonymous transitions and nonsynonymous transversions.
    Where the codons differ at more than one position, each order in which
    the changes could happen is a path; a path through a stop codon is
    impossible, and the others are weighted by their probabilities under
    ``P``, the matrix of codon substitution probabilities indexed as
    ``codons``. This is ``CountDiffs()`` in ``yn00.c`` of PAML.
    """
    TV = [0, 0, 0, 0]
    if codon1 == codon2:
        return TV
    diff_pos = [
        i
        for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
        if nucleotide1 != nucleotide2
    ]
    stop_codons = codon_table.stop_codons
    amino_acids = codon_table.forward_table
    transitions = ({"A", "G"}, {"C", "T"})
    paths = []  # probability and differences of each possible path
    for order in permutations(diff_pos):
        prob = 1
        diffs = [0, 0, 0, 0]
        codon = codon1
        for pos in order:
            next_codon = codon[:pos] + codon2[pos] + codon[pos + 1 :]
            if next_codon in stop_codons:
                break
            prob *= P[codons.index(codon), codons.index(next_codon)]
            nonsynonymous = amino_acids[codon] != amino_acids[next_codon]
            transversion = {codon[pos], codon2[pos]} not in transitions
            diffs[2 * nonsynonymous + transversion] += 1
            codon = next_codon
        else:  # no stop codon on the way
            paths.append((prob, diffs))
    if not paths:
        # Every path passes through a stop codon, which happens in some codon
        # tables. As yn00.c, count all differences as nonsynonymous, half a
        # transition and the rest transversions.
        return [0, 0, 0.5, len(diff_pos) - 0.5]
    total = sum(prob for prob, diffs in paths)
    for prob, diffs in paths:
        weight = prob / total
        TV = [m + n * weight for m, n in zip(TV, diffs)]
    return TV


#################################################################
#        private functions for Maximum Likelihood method
#################################################################


def _ml(codons1, codons2, cmethod, codon_table):
    """ML method main function (PRIVATE)."""
    from scipy.optimize import minimize

    pi = _get_pi(codons1, codons2, cmethod, codon_table=codon_table)
    codon_cnt = Counter(zip(codons1, codons2))
    codons = [
        codon
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons
        if "U" not in codon
    ]

    # apply optimization
    def func(
        params, pi=pi, codon_cnt=codon_cnt, codons=codons, codon_table=codon_table
    ):
        """Temporary function, params = [t, k, w]."""
        return -_likelihood_func(
            params[0],
            params[1],
            params[2],
            pi,
            codon_cnt,
            codons=codons,
            codon_table=codon_table,
        )

    # count sites
    opt_res = minimize(
        func,
        [1, 0.1, 2],
        method="L-BFGS-B",
        bounds=((1e-10, 20), (1e-10, 20), (1e-10, 10)),
        tol=1e-5,
    )
    t, k, w = opt_res.x
    Q = _get_Q(pi, k, w, codons, codon_table)
    Sd = Nd = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if i != j:
                try:
                    if (
                        codon_table.forward_table[codon1]
                        == codon_table.forward_table[codon2]
                    ):
                        # synonymous count
                        Sd += pi[codon1] * Q[i, j]
                    else:
                        # nonsynonymous count
                        Nd += pi[codon1] * Q[i, j]
                except KeyError:
                    # This is probably due to stop codons
                    pass
    Sd *= t
    Nd *= t

    # count differences (with w fixed to 1)
    def func_w1(
        params, pi=pi, codon_cnt=codon_cnt, codons=codons, codon_table=codon_table
    ):
        """Temporary function, params = [t, k]. w is fixed to 1."""
        return -_likelihood_func(
            params[0],
            params[1],
            1.0,
            pi,
            codon_cnt,
            codons=codons,
            codon_table=codon_table,
        )

    opt_res = minimize(
        func_w1,
        [1, 0.1],
        method="L-BFGS-B",
        bounds=((1e-10, 20), (1e-10, 20)),
        tol=1e-5,
    )
    t, k = opt_res.x
    w = 1.0
    Q = _get_Q(pi, k, w, codons, codon_table)
    rhoS = rhoN = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if i != j:
                try:
                    if (
                        codon_table.forward_table[codon1]
                        == codon_table.forward_table[codon2]
                    ):
                        # synonymous count
                        rhoS += pi[codon1] * Q[i, j]
                    else:
                        # nonsynonymous count
                        rhoN += pi[codon1] * Q[i, j]
                except KeyError:
                    # This is probably due to stop codons
                    pass
    rhoS *= 3
    rhoN *= 3
    dN = Nd / rhoN
    dS = Sd / rhoS
    return dN, dS


def _get_pi(codons1, codons2, cmethod, codon_table):
    """Obtain codon frequency dict (pi) from two codon list (PRIVATE).

    This function is designed for ML method. Available counting methods
    (cfreq) are F1x4, F3x4 and F64.
    """
    # TODO:
    # Stop codon should not be allowed according to Yang.
    # Try to modify this!
    pi = {}
    if cmethod == "F1x4":
        # A base absent from both sequences has frequency zero, as in F3x4.
        fcodon = {"A": 0, "G": 0, "C": 0, "T": 0}
        for codon in codons1 + codons2:
            for nucleotide in codon:
                fcodon[nucleotide] += 1
        tot = sum(fcodon.values())
        fcodon = {j: k / tot for j, k in fcodon.items()}
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = fcodon[codon[0]] * fcodon[codon[1]] * fcodon[codon[2]]
    elif cmethod == "F3x4":
        # three codon position
        fcodon = [
            {"A": 0, "G": 0, "C": 0, "T": 0},
            {"A": 0, "G": 0, "C": 0, "T": 0},
            {"A": 0, "G": 0, "C": 0, "T": 0},
        ]
        for codon in codons1 + codons2:
            fcodon[0][codon[0]] += 1
            fcodon[1][codon[1]] += 1
            fcodon[2][codon[2]] += 1
        for i in range(3):
            tot = sum(fcodon[i].values())
            fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = (
                    fcodon[0][codon[0]] * fcodon[1][codon[1]] * fcodon[2][codon[2]]
                )
    elif cmethod == "F61":
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = 0.1
        for codon in codons1 + codons2:
            pi[codon] += 1
        tot = sum(pi.values())
        pi = {j: k / tot for j, k in pi.items()}
    return pi


def _q(codon1, codon2, pi, k, w, codon_table):
    """Q matrix for codon substitution (PRIVATE).

    Arguments:
     - codon1, codon2  : three letter codon string
     - pi              : expected codon frequency
     - k               : transition/transversion ratio
     - w               : nonsynonymous/synonymous rate ratio
     - codon_table     : Bio.Data.CodonTable object

    """
    if codon1 == codon2:
        # diagonal elements is the sum of all other elements
        return 0
    if codon1 in codon_table.stop_codons or codon2 in codon_table.stop_codons:
        return 0
    if (codon1 not in pi) or (codon2 not in pi):
        return 0
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    diff = [
        (i, nucleotide1, nucleotide2)
        for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
        if nucleotide1 != nucleotide2
    ]
    if len(diff) >= 2:
        return 0
    if codon_table.forward_table[codon1] == codon_table.forward_table[codon2]:
        # synonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return k * pi[codon2]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return k * pi[codon2]
        else:
            # transversion
            return pi[codon2]
    else:
        # nonsynonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return w * k * pi[codon2]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return w * k * pi[codon2]
        else:
            # transversion
            return w * pi[codon2]


def _get_Q(pi, k, w, codons, codon_table):
    """Q matrix for codon substitution (PRIVATE)."""
    codon_num = len(codons)
    Q = np.zeros((codon_num, codon_num))
    for i1, codon1 in enumerate(codons):
        for i2, codon2 in enumerate(codons):
            if i1 != i2:
                Q[i1, i2] = _q(codon1, codon2, pi, k, w, codon_table=codon_table)
    nucl_substitutions = 0
    for i, codon in enumerate(codons):
        Q[i, i] = -sum(Q[i, :])
        try:
            nucl_substitutions += pi[codon] * (-Q[i, i])
        except KeyError:
            pass
    Q /= nucl_substitutions
    return Q


def _likelihood_func(t, k, w, pi, codon_cnt, codons, codon_table):
    """Likelihood function for ML method (PRIVATE)."""
    from scipy.linalg import expm

    Q = _get_Q(pi, k, w, codons, codon_table)
    P = expm(Q * t)
    likelihood = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if (codon1, codon2) in codon_cnt:
                if P[i, j] * pi[codon1] <= 0:
                    likelihood += codon_cnt[(codon1, codon2)] * 0
                else:
                    likelihood += codon_cnt[(codon1, codon2)] * log(
                        pi[codon1] * P[i, j]
                    )
    return likelihood


def _dn_ds_for_pair(name1, name2, function, *args, **kwargs):
    """Call function for a pair of a matrix, naming the pair in its warnings (PRIVATE).

    A warning from calculate_dn_ds or cal_dn_ds does not say which pair it
    is about, and the default filter of Python shows a repeated message only
    once, so a matrix of many pairs could show one warning for several.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = function(*args, **kwargs)
    for warning in caught:
        warnings.warn(
            f"{name1} and {name2}: {warning.message}", warning.category, stacklevel=3
        )
    return result


def calculate_dn_ds_matrix(alignment, method="NG86", codon_table=None):
    """Calculate dN and dS pairwise for the multiple alignment, and return as matrices.

    Argument:
     - method       - Available methods include NG86, LWL85, YN00 and ML.
     - codon_table  - Codon table to use for forward translation.

    A distance that ``calculate_dn_ds`` finds undefined, because the pair is
    too divergent for the method, is nan in the matrix, with a
    ``BiopythonWarning`` for each such pair that starts with the names of
    the two sequences; the other distances are still calculated.

    """
    from Bio.Phylo.TreeConstruction import DistanceMatrix

    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    sequences = alignment.sequences
    coordinates = alignment.coordinates
    names = [record.id for record in sequences]
    size = len(names)
    dn_matrix = []
    ds_matrix = []
    for i in range(size):
        dn_matrix.append([])
        ds_matrix.append([])
        for j in range(i):
            pairwise_sequences = [sequences[i], sequences[j]]
            pairwise_coordinates = coordinates[(i, j), :]
            pairwise_alignment = Alignment(pairwise_sequences, pairwise_coordinates)
            dn, ds = _dn_ds_for_pair(
                names[i],
                names[j],
                calculate_dn_ds,
                pairwise_alignment,
                method=method,
                codon_table=codon_table,
            )
            dn_matrix[i].append(dn)
            ds_matrix[i].append(ds)
        dn_matrix[i].append(0.0)
        ds_matrix[i].append(0.0)
    dn_dm = DistanceMatrix(names, matrix=dn_matrix)
    ds_dm = DistanceMatrix(names, matrix=ds_matrix)
    return dn_dm, ds_dm


def mktest(alignment, species=None, codon_table=None):
    """McDonald-Kreitman test for neutrality.

    Implement the McDonald-Kreitman test for neutrality (PMID: 1904993)
    This method counts changes rather than sites
    (http://mkt.uab.es/mkt/help_mkt.asp).

    Arguments:
     - alignment    - Alignment of gene nucleotide sequences to compare.
     - species      - List of the species ID for each sequence in the alignment.
       Typically, the species ID is the species name as a string, or an integer.
     - codon_table  - Codon table to use for forward translation.

    Return the p-value of test result.
    """
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    unique_species = set(species)
    sequences = []
    for sequence in alignment.sequences:
        try:
            sequence = sequence.seq
        except AttributeError:
            pass
        sequence = str(sequence)
        sequences.append(sequence)
    columns = []
    starts = sys.maxsize
    for ends in alignment.coordinates.transpose():
        step = min(ends - starts)
        for j in range(0, step, 3):
            codons = {key: set() for key in unique_species}
            for key, sequence, start in zip(species, sequences, starts):
                codons[key].add(sequence[start + j : start + j + 3])
            columns.append(list(codons.values()))
        starts = ends
    return _mktest(columns, codon_table)


def _mktest(columns, codon_table):
    """Count substitutions in codon columns and apply the G test (PRIVATE).

    This is the McDonald-Kreitman test proper, shared by ``mktest`` above and
    by ``Bio.codonalign.mktest``. Each column describes one codon position of
    the alignment, as a list holding, for each species, the set of codons seen
    in that species. A column with only one codon is invariant and skipped.
    The changes in a column count as fixed if every species has one codon,
    and as polymorphic otherwise.
    """
    G, nonsyn_G = _get_codon2codon_matrix(codon_table=codon_table)
    syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = 0, 0, 0, 0
    for column in columns:
        all_codons = set().union(*column)
        if len(all_codons) == 1:
            continue
        nonsyn = _count_replacement(all_codons, nonsyn_G)
        syn = _count_replacement(all_codons, G) - nonsyn
        if all(len(codons) == 1 for codons in column):
            # fixed
            nonsyn_fix += nonsyn
            syn_fix += syn
        else:
            # not fixed
            nonsyn_poly += nonsyn
            syn_poly += syn
    return _G_test([syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])


def _get_codon2codon_matrix(codon_table):
    """Get codon codon substitution matrix (PRIVATE).

    Elements in the matrix are number of synonymous and nonsynonymous
    substitutions required for the substitution.
    """
    bases = ("A", "T", "C", "G")
    codons = [
        codon
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons
        if "U" not in codon
    ]
    # set up codon_dict considering stop codons
    codon_dict = codon_table.forward_table.copy()
    for stop in codon_table.stop_codons:
        codon_dict[stop] = "stop"
    # count site
    num = len(codons)
    G = {}  # graph for substitution
    nonsyn_G = {}  # graph for nonsynonymous substitution
    graph = {}
    graph_nonsyn = {}
    for i, codon in enumerate(codons):
        graph[codon] = {}
        graph_nonsyn[codon] = {}
        for p in range(3):
            for base in bases:
                tmp_codon = codon[0:p] + base + codon[p + 1 :]
                if codon_dict[codon] != codon_dict[tmp_codon]:
                    graph_nonsyn[codon][tmp_codon] = 1
                    graph[codon][tmp_codon] = 1
                else:
                    if codon != tmp_codon:
                        graph_nonsyn[codon][tmp_codon] = 0.1
                        graph[codon][tmp_codon] = 1
    for codon1 in codons:
        nonsyn_G[codon1] = {}
        G[codon1] = {}
        for codon2 in codons:
            if codon1 == codon2:
                nonsyn_G[codon1][codon2] = 0
                G[codon1][codon2] = 0
            else:
                nonsyn_G[codon1][codon2] = _dijkstra(graph_nonsyn, codon1, codon2)
                G[codon1][codon2] = _dijkstra(graph, codon1, codon2)
    return G, nonsyn_G


def _dijkstra(graph, start, end):
    """Dijkstra's algorithm Python implementation (PRIVATE).

    Algorithm adapted from
    http://thomas.pelletier.im/2010/02/dijkstras-algorithm-python-implementation/.
    However, an obvious bug in::

        if D[child_node] >(<) D[node] + child_value:

    is fixed.
    This function will return the distance between start and end.

    Arguments:
     - graph: Dictionary of dictionary (keys are vertices).
     - start: Start vertex.
     - end: End vertex.

    Output:
       List of vertices from the beginning to the end.

    """
    D = {}  # Final distances dict
    P = {}  # Predecessor dict
    # Fill the dicts with default values
    for node in graph.keys():
        D[node] = 100  # Vertices are unreachable
        P[node] = ""  # Vertices have no predecessors
    D[start] = 0  # The start vertex needs no move
    unseen_nodes = list(graph.keys())  # All nodes are unseen
    while len(unseen_nodes) > 0:
        # Select the node with the lowest value in D (final distance)
        shortest = None
        node = ""
        for temp_node in unseen_nodes:
            if shortest is None:
                shortest = D[temp_node]
                node = temp_node
            elif D[temp_node] < shortest:
                shortest = D[temp_node]
                node = temp_node
        # Remove the selected node from unseen_nodes
        unseen_nodes.remove(node)
        # For each child (ie: connected vertex) of the current node
        for child_node, child_value in graph[node].items():
            if D[child_node] > D[node] + child_value:
                D[child_node] = D[node] + child_value
                # To go to child_node, you have to go through node
                P[child_node] = node
        if node == end:
            break
    # Set a clean path
    path = []
    # We begin from the end
    node = end
    distance = 0
    # While we are not arrived at the beginning
    while not (node == start):
        if path.count(node) == 0:
            path.insert(0, node)  # Insert the predecessor of the current node
            node = P[node]  # The current node becomes its predecessor
        else:
            break
    path.insert(0, start)  # Finally, insert the start vertex
    for i in range(len(path) - 1):
        distance += graph[path[i]][path[i + 1]]
    return distance


def _count_replacement(codons, G):
    """Count replacement needed for a given codon_set (PRIVATE)."""
    if len(codons) == 1:
        return 0, 0
    elif len(codons) == 2:
        codons = list(codons)
        return floor(G[codons[0]][codons[1]])
    else:
        subgraph = {
            codon1: {codon2: G[codon1][codon2] for codon2 in codons if codon1 != codon2}
            for codon1 in codons
        }
        return _prim(subgraph)


def _prim(G):
    """Prim's algorithm to find minimum spanning tree (PRIVATE).

    Code is adapted from
    http://programmingpraxis.com/2010/04/09/minimum-spanning-tree-prims-algorithm/
    """
    nodes = []
    edges = []
    for i in G.keys():
        nodes.append(i)
        for j in G[i]:
            if (i, j, G[i][j]) not in edges and (j, i, G[i][j]) not in edges:
                edges.append((i, j, G[i][j]))
    conn = defaultdict(list)
    for n1, n2, c in edges:
        conn[n1].append((c, n1, n2))
        conn[n2].append((c, n2, n1))
    mst = []  # minimum spanning tree
    used = {nodes[0]}
    usable_edges = conn[nodes[0]][:]
    heapify(usable_edges)
    while usable_edges:
        cost, n1, n2 = heappop(usable_edges)
        if n2 not in used:
            used.add(n2)
            mst.append((n1, n2, cost))
            for e in conn[n2]:
                if e[2] not in used:
                    heappush(usable_edges, e)
    length = 0
    for p in mst:
        length += floor(p[2])
    return length


def _G_test(site_counts):
    """G test for 2x2 contingency table (PRIVATE).

    Arguments:
     - site_counts - [syn_fix, nonsyn_fix, syn_poly, nonsyn_poly]

    Return the upper tail probability of ``G = 2 sum(O ln(O/E))`` under a
    chi-square distribution with one degree of freedom.

    Empty cells are handled as follows. An observed count of zero contributes
    nothing to the sum, following the usual convention that ``lim x->0 of
    x ln(x) = 0``. That is enough to keep the sum defined, because an expected
    count can only be zero when its row or column total is zero, and then the
    observed count in that cell must be zero as well. A table with an empty
    row or column is degenerate: every remaining observed count then equals
    its expected count exactly, so ``G`` is zero and the returned p-value is
    one. The same value is returned for an all-zero table, for which the
    expected counts are not defined at all; this is the table produced by an
    alignment with no substitutions to count.

    >>> print("%0.6f" % _G_test([17, 7, 42, 2]))
    0.004924
    >>> print("%0.6f" % _G_test([0, 0, 0, 0]))
    1.000000
    >>> print("%0.6f" % _G_test([3, 5, 0, 0]))
    1.000000
    """
    # TODO:
    #   Apply continuity correction for Chi-square test.
    G = 0
    tot = sum(site_counts)
    if tot == 0:
        # No counts at all; the expected counts are undefined, and there is
        # nothing for the test to detect.
        return 1.0
    tot_syn = site_counts[0] + site_counts[2]
    tot_non = site_counts[1] + site_counts[3]
    tot_fix = sum(site_counts[:2])
    tot_poly = sum(site_counts[2:])
    exp = [
        tot_fix * tot_syn / tot,
        tot_fix * tot_non / tot,
        tot_poly * tot_syn / tot,
        tot_poly * tot_non / tot,
    ]
    for obs, ex in zip(site_counts, exp):
        if obs:
            # ex is strictly positive whenever obs is, as obs contributes to
            # both the row total and the column total that define ex.
            G += obs * log(obs / ex)
    # with only 1 degree of freedom for a 2x2 table,
    # the cumulative chi-square distribution reduces to a simple form:
    return erfc(sqrt(G))


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
