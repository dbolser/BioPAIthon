# Upstream reports

Most of what this fork fixes is also present in Biopython, because most of it
predates the fork. Where that is true the defect is worth reporting upstream,
so that people who use Biopython benefit whether or not they ever hear of this
project.

This file records what has been reported, what is queued, and what was checked
and deliberately not reported. It is a working record rather than a plan of
campaign: entries leave the queue as often because they turn out to be already
known as because they get filed.

## How reports are made

- **As information, not as pull requests.** A report carries a reproduction, a
  root cause and a suggested fix, and offers a PR rather than opening one.
  Upstream sets its own contribution terms and this fork does not test them.
- **One at a time.** A queue of verified defects is not a reason to file a
  queue of issues. Reports are spaced, and the reception of the last one
  informs whether to send the next.
- **Verified immediately before filing.** Upstream is active. Every report is
  re-checked against current `master` on the day it is sent, because a defect
  found months ago may have been fixed since.
- **Calibrated.** Each report states what was run and what was only read, and
  names what was not checked. Where a fix is a behaviour choice rather than a
  correction, the report says so and leaves the choice upstream.
- **Attributed.** Reports produced with AI assistance say so.

## Reported

| upstream | date | subject | fixed here |
|---|---|---|---|
| [#3771](https://github.com/biopython/biopython/issues/3771) | 2026-07-30 | Use-after-free in the `cpairwise2` sequence conversion; root cause for a crash open since 2021 | [#33](https://github.com/dbolser/BioPAIthon/pull/33) |
| [#5269](https://github.com/biopython/biopython/issues/5269) | 2026-07-30 | `run_cealign` segfaults on a zero fragment size, short coordinate entries, or tuples | no |

## Queued

Verified unreported in the GitHub tracker and still present on upstream
`master`. Ordered by how self-evident the defect is, not by severity.

| # | subject | note |
|---|---|---|
| 1 | `Bio.File.as_handle` re-catches a `TypeError` from the caller's own block | Fixed upstream in 2018 for [#1544](https://github.com/biopython/biopython/issues/1544), reverted by `f09b169a4` in 2020, no regression test. Comment on #1544 rather than a new issue. |
| 2 | `_get_pi` builds `forward_table.keys() + stop_codons` | `TypeError`; the `F1x4` and `F61` codon frequency models are unusable. |
| 3 | `Alignment.format()` reports any `AttributeError` as "not yet implemented" | Upstream's own `Bio.Align.write()` already does this correctly. |
| 4 | `CodonAdaptationIndex` rejects `MutableSeq` | The docstring promises it is accepted. |
| 5 | GAF by-protein iterators never yield the last protein | Silent data loss. Fixing it changes output for every caller. |
| 6 | `Seq.search()` never matches at the final position | Off by one since the method shipped. |
| 7 | `six_frame_translations` misaligns the reverse frames | Correct only when `len(seq) % 3 == 1`. Display only. |
| 8 | `_G_test` divides by zero on any empty cell | `mktest` crashes on identical sequences. |
| 9 | `nt_search(seq, "")` never terminates | The fix is an API choice; report, do not prescribe. |
| 10 | `Tm_GC(valueset=0)` raises `UnboundLocalError` | The guard admits a value no branch handles. |
| 11 | `str()` on an `index_db` dictionary raises `AttributeError` | Upstream carries a `TODO` about this exact question; lead with it. |
| 12 | `SeqIO.index()` raises `IndexError` on a bare `>` FASTA header | `parse()` reads the same record as `id=""`. |
| 13 | Two `except ValueError` clauses in `SeqIO/_index.py` cannot fire | File the `tab` and `genbank` halves separately; cite [#1344](https://github.com/biopython/biopython/issues/1344). |
| 14 | `reverse_complement()` cannot sort features with an `UnknownPosition` | Cite [#1772](https://github.com/biopython/biopython/issues/1772), same cause, different path. |
| 15 | `ProteinAnalysis.flexibility()` reads the wrong window centre | Comment on [#4170](https://github.com/biopython/biopython/pull/4170), which fixes the adjacent window-count bug only. |

## Needs work before it can be filed

- **`translate()` drops an explicit `gap=` for `Seq` input.** Open PR
  [#4994](https://github.com/biopython/biopython/pull/4994) changes the
  module-level default and does not forward the argument, so this survives it.
  Any report must say so or it will be closed as a duplicate.
- **`SeqIO.index()` masks a proxy's `TypeError`.** The only reproduction is a
  patched SFF header. Needs a natural malformed input first.
- **Location truthiness.** The dead Python 2 `__nonzero__` methods and the
  stale docstring are real, but the resolution is a behaviour choice, and
  either direction touches `Bio/SeqIO/SnapGeneIO.py:180`.
- **`cpairwise2.rint` parses `int` with `"l"`.** A genuine LP64 stack
  overwrite, but reachable only from the deprecated `Bio.pairwise2`.

## Held pending a disclosure decision

Two findings are driven by untrusted input rather than programmer error, and
upstream's `.github/SECURITY.md` asks for those to go to the module maintainer
rather than the public tracker. Nothing is filed for either.

- **`Bio/PDB/bcifhelpermodule.c` writes past the output buffer.** The unpack
  helpers bound their loop by the input only; the destination is sized from
  the file's own `srcSize`. A `.bcif` file whose two fields disagree overflows
  the array.
- **`Bio/Align/_aligncore.c` returns heap contents.** `feed()` adds an
  unvalidated `offset` to the buffer, and the length-mismatch guard is a
  tautology, so `Alignment.parse_printed_alignment()` on a ragged block
  returns process memory in the coordinates array instead of raising.

## Checked and not reported

- **`ProteinAnalysis.flexibility()` window count.** Already open as
  [#4170](https://github.com/biopython/biopython/pull/4170) since 2022.
- **The `_prim` spanning-tree seed.** Already open as
  [#2301](https://github.com/biopython/biopython/pull/2301) since 2019, with
  the same one-line fix and maintainer confirmation.
- **`Bio.PDB.ccealign` hanging on specific structures**
  ([#4888](https://github.com/biopython/biopython/issues/4888)). Reproduced on
  1.83 and confirmed fixed from 1.84, but the change was a wholesale rewrite of
  the module rather than a single commit, so there is nothing useful to add.
- **The zero-width-space crash** ([#3387](https://github.com/biopython/biopython/issues/3387),
  closed). Plausibly the same defect as #3771 but not reproducible on any
  version tried, so no claim was made.

## Limits of the survey

The duplicate check covered the GitHub issue and pull request tracker only. It
did not cover the `biopython-dev` mailing list or the pre-GitHub Redmine and
Bugzilla archives, so "unreported" here means "not in the GitHub tracker".
