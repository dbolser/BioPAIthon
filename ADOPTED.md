# Work adopted from upstream

Biopython has 142 open pull requests, some of them years old. A number are
finished, reviewed and approved, and are held up only because a maintainer
asked whether an AI tool was involved. This fork does not ask that question, so
it can take them.

This file records every one, who wrote it, what was taken, and what was left
behind. Nothing here is this fork's work.

## How this is done, and why the history is kept

Adopted commits are **cherry-picked, not retyped**. Git records the original
author in the commit's `Author` field and this fork only as `Committer`, so
`git log --format='%an'` names the person who actually wrote the code and
`git blame` keeps pointing at them. Authors are added to `CONTRIB.rst`.

The licence basis is upstream's own pull request template, which every
contributor ticks:

> I hereby agree to dual licence this and any previous contributions under
> both the *Biopython License Agreement* **AND** the *BSD 3-Clause License*.

That makes the work BSD-3 and adoptable here. **The box is checked on each
pull request before anything is taken from it** — a few upstream PRs leave it
unticked, and those are not eligible however good they are.

Adopting is not endorsing. Where a pull request grew a change under review that
this fork does not want, only the wanted commits are taken, and the omission is
recorded below with the reason.

**An adoption pull request adds its own entry here.** This file went four
adoptions out of date before anyone noticed, which is easy to do and bad in a
specific way: `CONTRIB.rst` and `git log` kept the authorship, but the record of
*what was taken and what was left behind* is only in this file, and a partial
adoption that nobody wrote down looks like an endorsement of the whole pull
request.

## Adopted

Newest first.

### [#5244](https://github.com/biopython/biopython/pull/5244) — Gian Carlo D Dumpit

`Bio.Align.sam.AlignmentWriter` hard-coded `tLen = 0`, so no SAM file the
library wrote ever carried a template length, even when the parser had read
one in (upstream issue #5232). The writer now honours an explicit
`alignment.tlen` attribute — which the parser has always set for a nonzero
TLEN, so round-trips keep the value — and, when the alignments are supplied as a list or tuple,
computes TLEN for properly paired reads: same query name, paired flag (0x1)
set, both mates mapped on the same reference, each record's `pnext` matching
its mate's start; unmapped (0x4/0x8), secondary (0x100) and supplementary
(0x800) records are excluded from pairing. The value is the span from the
leftmost mapped base to the rightmost, signed `+` for the leftmost mate and
`-` for the rightmost as the SAM specification requires, with a tie on start
position broken by end coordinate.

**The single commit was taken whole; nothing was omitted.** Two limitations
discussed in review remain exactly as upstream left them: alignments supplied
as an iterator are still written streaming with TLEN 0 unless `.tlen` is set,
because pairing needs the mate's coordinates; and the pairing pass keeps one
small tuple per paired alignment in memory for the duration of the write —
the maintainer sketched a name-sorted streaming alternative for very large
files that was never implemented. The author was already in `CONTRIB.rst`,
having contributed upstream #5220, which arrived here with Biopython 1.88.

### [#3812](https://github.com/biopython/biopython/pull/3812) — Dominique Sydow

`HSExposureCB._get_cb` returned `self._get_gly_cb_vector(r2), 0.0` for a
glycine. `_get_gly_cb_vector` returns `None` when the residue is missing N, C
or CA, so the tuple became `(None, 0.0)` — which is not `None`, so the caller's
`if result is None: continue` does not fire, and `pcb` goes on to be used as a
vector. Returning `None` outright lets the existing guard do its job.

The cost of not doing so is the whole calculation, not one residue. On
`Tests/PDB/a_structure.pdb` with the N removed from a single glycine,
`HSExposureCB` populates **2 of 86 residues** before dying with
`AttributeError: 'NoneType' object has no attribute 'norm'`. With the fix, 84.

**One of the six commits taken.** The other five add a `skip_residues=False`
keyword to `_AbstractHSExposure`, `HSExposureCA` and `HSExposureCB`, and make
the default behaviour `raise KeyError` where the code currently skips the
residue and carries on. That is a design decision rather than a bug fix: every
existing caller that today gets a slightly short result set would start getting
an exception instead, and the argument for that is about API taste, not
correctness. The narrow fix already delivers what
`IMPROVEMENTS.md`'s harvest table asked for. If the fork later wants the
opt-in strictness, those commits are still there to take.

The `except Exception` in `_get_gly_cb_vector` is also left as it is; the
pull request narrows it to `KeyError` and adds a warning, which is a separate
improvement worth its own change.

### [#5127](https://github.com/biopython/biopython/pull/5127) — Michiel de Hoon

`multi_coord_space` computes the polar angle as `arccos(np.divide(p[:, 2], r,
where=r != 0))`. Without an `out=` argument, `np.divide` allocates its result
with `np.empty`, so every entry the `where` mask skips — every hedron whose
second atom sits on its first — keeps whatever was in that memory. The polar
angle, and the transform built from it, are then read from uninitialised heap.
Passing `out=np.ones_like(r)` makes the quotient 1 there, so the angle is zero,
which is what the scalar sibling `get_spherical_coordinates` already returns for
`r == 0`.

**One of the two commits taken.** The second, `2e4221b44` ("update"), changes
four call sites in `coord_space` from `a1[0]` to `a1[0][0]` and `sc[2]` to
`sc[2][0]`, following a stale comment that records the indexing from before the
callers switched to flat arrays. Applied here it fails immediately:

```
IndexError: invalid index to scalar variable.
```

`Tests/test_PDB_vectors.py` passes `np.array([2.0, 0.0, 2.0, 1.0])`, so `a1[0]`
is already a scalar. This is not a difference between the fork and upstream —
upstream `master` at `e136be720` has the same flat-array test and the same
`a1[0]` call sites, so that commit would break upstream's own suite as well.
They already know: peterjc posted mypy's output on the pull request on
2025-12-27, naming the same lines, and the thread has not moved since. See
`UPSTREAM.md` for the two details it does not have.

### [#4866](https://github.com/biopython/biopython/pull/4866) — Manuel Lera-Ramirez

Fixes upstream issue
[#4865](https://github.com/biopython/biopython/issues/4865). Every
`Bio.SeqUtils.MeltingTemp.Tm_*` method began with `seq = str(seq)`. For a
`SeqRecord`, `str()` is the multi-line human-readable summary rather than the
sequence, and `_check` then stripped that summary down to the characters in the
method's base set and melted the result — no error, no warning, a plausible
wrong number. Five commits, all taken; three corrections for this tree sit on
top as separate commits. Taken in fork PR
[#45](https://github.com/dbolser/BioPAIthon/pull/45).

### [#4450](https://github.com/biopython/biopython/pull/4450) — Gail Bartlett

Closes upstream issue
[#4449](https://github.com/biopython/biopython/issues/4449).
`Bio/PDB/parse_pdb_header.py::_format_date` looked the month up with a bare
`all_months.index(pdb_date[3:6])`, so any PDB file whose `HEADER` or `REVDAT`
line carries a non-English month abbreviation killed the entire parse, with an
exception naming no file, no line and no field. Six commits, all taken, plus
three corrections on top.

peterjc approved it on 2023-09-25 and asked for a rebase on 2025-01-03; the
author never returned. Taken in fork PR
[#44](https://github.com/dbolser/BioPAIthon/pull/44).

### [#3897](https://github.com/biopython/biopython/pull/3897) — Bart Grosman

Closes upstream issue
[#3896](https://github.com/biopython/biopython/issues/3896), open since April
2022. `DisorderedEntityWrapper.copy()` takes a shallow copy, empties
`child_dict` and re-adds copies of the children through `disordered_add` —
which only calls `disordered_select` when a child's occupancy is *strictly
greater* than `self.last_occupancy`. The shallow copy inherited
`last_occupancy` already at the maximum, so no child was ever selected and the
copy's `selected_child` still pointed into the original: reading or
transforming the copy read and transformed the original's coordinates. Five
commits taken, plus corrections on top. Taken in fork PR
[#43](https://github.com/dbolser/BioPAIthon/pull/43).

### [#5181](https://github.com/biopython/biopython/pull/5181) — Ernest Provo

`Motif.calculate_consensus` raised `UnboundLocalError` on a column where every
count is zero: `consensus_letter` is only assigned inside a `count > maximum`
comparison that such a column never satisfies, and the degenerate-column test
then divided by a zero total.

**One of the four commits taken.** Under review the pull request grew alphabet
validation in `Motif.__init__` that rejects any letter outside the declared
alphabet, and the author noted the consequence himself on 2026-06-05:
`Motif("ACGT", Alignment([Seq("ACGN"), Seq("ACGT")]))` starts raising.
Ambiguity codes in a DNA alignment are routine input, so that half would break
working code to fix an unrelated problem. Taken in fork PR
[#42](https://github.com/dbolser/BioPAIthon/pull/42).

### [#5175](https://github.com/biopython/biopython/pull/5175) — Abdel ATIA

Replaces nine `assert` statements that validate input in production code with
explicit `if`/`raise`, in `Bio/ExPASy/Enzyme.py`, `Bio/ExPASy/ScanProsite.py`,
`Bio/SeqIO/TabIO.py` and `Bio/motifs/alignace.py`. Under `python -O` every one
of those checks vanishes, so a malformed file produces silently wrong output
rather than an error — the same defect `IMPROVEMENTS.md` §0.8 describes across
343 sites. Four commits, all taken, tests included.

Approved by mdehoon on 2026-03-18 and peterjc on 2026-04-27 ("The changes look
good, thank you"). The requested `assertRaises` tests were added on 2026-05-01.
The thread then stopped at "is something funny about your recent comments - are
you using a coding agent LLM AI tool?" and has not moved since.
