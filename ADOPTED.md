# Work adopted from upstream

Biopython had 142 open pull requests on 2026-08-03, when this file was
started (162 on 2026-10-22), some of them years old. A number are
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

### [#5391](https://github.com/biopython/biopython/pull/5391) — Rubén Chaves

Fixes upstream issue
[#1383](https://github.com/biopython/biopython/issues/1383), open since 2017.
`PDBIO.save` wrote every `TER` record 81 columns wide, one more than a PDB
record, and in its default numbering gave the `TER` the serial number that
the first atom of the next chain then reused: `ATOM 371`, `TER 372`,
`ATOM 372` for `Tests/PDB/2BEG.pdb`. The fix pads `TER` to 80 columns and
advances the counter after it, the second patch proposed in the 2017 thread.
Measured here, writing 2BEG back out now reproduces the deposited serials of
all 1860 `ATOM` and `TER` records, where 1488 differed before. The pull
request also documents `preserve_atom_numbering` in the `save` docstring.

**The single commit taken, plus two corrections on top.** Its `NEWS.rst`
paragraph was written for upstream's in-progress 1.89 section, which in this
tree is preserved upstream history, so the cherry-pick moves it into the
fork's own section and leaves upstream's contributor list alone; the author
is in `CONTRIB.rst`. The first correction makes the new test read its fixture
through `support.DATA`: it opened `PDB/1SSU_mod.pdb` relative to the current
directory, which fails when pytest starts outside `Tests/`, as CI does.

The second correction rewords the end of the new `preserve_atom_numbering`
docstring. Upstream's text said that the stored serials of files read with
`PDBParser` leave room for each `TER`. Files from the PDB leave that room
after each chain's polymer, but `PDBIO` writes the chain's waters and ligands
before its `TER`, so the `TER` can take a serial that the file gave to a
later atom. Of the test fixtures, 1LCD, 2XHE and 7DDO do this: 1LCD has
`TER 1030` after the last water of chain B and `HETATM 1030` for the first
water of chain C. Biopython 1.88 and this fork's main write the same serials;
only the documentation claimed more.

What changes for users, checked by writing 1SSU_mod, 2BEG and 1LCD with and
without the patch: every `TER` line is one trailing space shorter. With
`preserve_atom_numbering=True` nothing else differs. In the default
numbering, every serial after a model's first `TER` also moves up by one per
`TER` before it, and nothing else differs. One consequence: in the default
numbering, a model whose atoms fit in the five-column serial field only
because `TER` shared their numbers now raises `PDBIOException`, as a model
with more than 99,999 atoms already did. A synthetic structure of 99,999
atoms in two chains raises where it used to be written. A fork commit on top
also makes a `TER` whose own serial would pass 99,999 raise
`PDBIOException`, instead of writing it six digits wide and so 81 columns
(the Codex review on #275 found that the extra `TER` serials could push the
last one over the limit).

JoaoRodrigues approved the pull request on 2026-11-06 and asked whether it
was written with AI; there is no reply on the thread as of 2026-11-09. His
one inline comment, on the `preserve_atom_numbering` docstring, reads: "This
will cause the next atom after the TER to have the same atom number as the
TER. I think this is OK, weird, but OK." That happens with
`preserve_atom_numbering=True` when the stored serials leave no room for the
`TER`, as in the fixtures above; it is what the second correction documents.
The default numbering no longer does it.

**Related, not taken here:**
[#5381](https://github.com/biopython/biopython/pull/5381) by Shubhodeep
Chatterjee adds a `write_ter` keyword to `PDBIO.save` (upstream issue
[#4661](https://github.com/biopython/biopython/issues/4661)). It edits the
same `TER` block and is a separate decision. Whoever takes it should know
that it no longer applies as its author tested it. Its test is also named
`test_pdbio_write_ter`, so one of the two needs renaming. Its test also
asserts that omitting `TER` leaves every other line unchanged, which was
true only while `TER` did not advance the serial. With this adoption,
dropping the `TER` either closes the gaps in the numbering, changing the
other lines, or leaves gaps where no `TER` is, and that choice has to be
made.

### [#5366](https://github.com/biopython/biopython/pull/5366) — popcoder565

Fixes upstream issue
[#4644](https://github.com/biopython/biopython/issues/4644). When a PDB
record's element column is missing or not a known element,
`Atom._assign_element` guesses the element from the atom name. A name
starting with a letter took the element of that letter only if it had more
than two characters and all of those after the second were digits (`HE21`);
any other was looked up whole as an element symbol. So four character
hydrogen names such as `HNZ1` and `HNE2`, which some force fields write for
lysine, came out as `X`. The fix tests only the last character, since no
element symbol ends in a digit.

**The single commit taken** (`bb836bd8da`, authored as "berth"), with its
`NEWS.rst` paragraph moved out of upstream's in-progress Biopython 1.89
section, which this tree keeps as upstream history, into the fork's own.
There was no review discussion. It supersedes
[#4645](https://github.com/biopython/biopython/pull/4645) by eunos-1128,
closed unmerged in 2024 without a test, whose `or` still sent `HNZ1` down the
whole-name branch. The author is not in `CONTRIB.rst` and ticked the box for
being listed already or not wishing to be, so is not added.

One correction sits on top, in `NEWS.rst` and the tests. The pull request
describes only the hydrogen case, and says no guessed element changes for any
file under `Tests/PDB`. That holds for the files as parsed. Only atoms whose
element column is missing or holds something that is not an element, such as
`X`, are guessed at all; in the fixtures that is four atoms, in
`a_structure.pdb` and `occupancy.pdb`, whose `CA` and `N` names give the same
element either way. With the element ignored, 25 of the 305 distinct fixture
names change, all unpadded mmCIF names such as `C1`, which go from `X` to the
right element. Some exotic names go from `X` to a first-letter guess: `CLA1`
carbon, `FEB2` fluorine, `BRC1` boron. The old rule already guessed that way
for `CL12`, but it is not strictly an improvement. No name that gave an element
before gives a different one or `X` now, checked over all 1,823,580 names of
one to four characters from A–Z and 0–9, padded and unpadded. The fork's NEWS
paragraph says so, and a new test pins the unpadded case.

### [#4390](https://github.com/biopython/biopython/pull/4390) — Thomas Holder

With `auth_residues=False`, `MMCIFParser` and `FastMMCIFParser` number
residues from `_atom_site.label_seq_id` but still attach the insertion code
from `_atom_site.pdbx_PDB_ins_code`, which is an author-scheme field — label
numbering has no insertion codes. The resulting ids mix the two schemes: in
`Tests/PDB/4ZHL.cif` chain U, the residue that is `(' ', 37, 'A')` in author
numbering came out as `(' ', 23, 'A')`, an id that exists in neither scheme
(the label scheme has a plain 23; the author scheme has 37 with icode A). The
fix blanks the insertion-code list when label numbering is in use, so that
residue is `(' ', 23, ' ')`.

**The single commit taken in full.** The pull request carried no review
discussion and grew no extras. Its test exercises only `MMCIFParser`, although
the commit patches both parsers; this fork adds a matching `FastMMCIFParser`
regression test in a follow-up commit.

### [#4634](https://github.com/biopython/biopython/pull/4634) — Brandon Seah

`Bio.Phylo` had no way to resolve polytomies, which many downstream tools
require — they insist on strictly bifurcating trees, while others (FastTree,
for one) happily produce multifurcations. The pull request adds
`TreeMixin.resolve_polytomies(shuffle=False, recursive=True)`, which converts
each polytomy into zero-length bifurcations, optionally shuffling the child
order first.

**All three commits taken**, including the author adding himself to
`CONTRIB.rst` — the only review comment (peterjc asking the test to compare
the resolved tree's string representation) was addressed by the author in the
third commit, and nothing was left behind. Corrections sit on top in separate
commits, and one changes behaviour for wide polytomies: the adopted code
resolved a polytomy into a caterpillar, one child per level, built by a
recursion that hit Python's recursion limit at ~1000 children — a star tree
being precisely the input the feature targets — and even built another way, a
chain that deep defeats the library's own recursive machinery
(`is_bifurcating`, the Newick writer, preorder `find_clades`). Polytomies now
resolve into a *balanced* hierarchy, ~log2(n) deep, preserving child order
left to right when `shuffle=False`. A three-child polytomy resolves exactly as
upstream's code resolved it; wider ones differ in shape (both are valid
zero-length resolutions), and the expected string in the adopted unit test is
updated to match. The public method also gained the doctest the house style
asks of new API.

### [#4918](https://github.com/biopython/biopython/pull/4918) — Erik Whiting

Fixes upstream issue
[#4879](https://github.com/biopython/biopython/issues/4879). The LOCUS line
checks in `Bio/GenBank/Scanner.py` accepted only molecule types containing
`DNA` or `RNA`, so `NA` — a nucleic acid of unspecified type, listed in
section 3.4.4.2 of the GenBank release notes and written by e.g. SnapGene
exports — was rejected with `ValueError: LOCUS line does not contain valid
sequence type`. The one adopted commit adds `NA` as an exact-match
alternative, deliberately not a substring test so that `99NA` and friends
stay rejected.

**The single commit taken, plus one correction on top.** As veghp pointed out
on the pull request, the scanner change alone does not make such files parse:
`_FeatureConsumer.record_end` in `Bio/GenBank/__init__.py` still raised
`Could not determine molecule_type for seq_type NA`. The correction commit
teaches `record_end` the `NA` branch, makes the same one-line change to the
writer's LOCUS checks in `Bio/SeqIO/InsdcIO.py` (by their own comment "copied
from Bio.GenBank.Scanner") so the record round-trips, and adds the regression
tests the upstream pull request was asked for but never grew.

### [#3911](https://github.com/biopython/biopython/pull/3911) — Erik Whiting

Fixes upstream issue
[#3885](https://github.com/biopython/biopython/issues/3885). SnapGene files
can carry literal newlines in feature qualifier values, and the SnapGene
parser passes them through to the `SeqRecord`. `GenBankWriter.write_record`
then wrote the qualifier across two lines with the continuation unindented —
a file the GenBank parser itself rejects with `ValueError`. The fix replaces
each newline with a space and warns with the qualifier key, record ID and
feature index, the form peterjc's review asked for. The related
[#4073](https://github.com/biopython/biopython/pull/4073) bundles the same
writer fix with a `Scanner.py` rewrite that changes line-joining for every
GenBank file read; this narrower shape is the one wanted.

**One commit taken** (`570c02c79`) — the branch's only other content is a
June 2022 merge of upstream master, long since part of this tree's history.
Three corrections sit on top as separate commits. First, the adopted check
tested `"\n" in feature.qualifiers[k]` behind an `isinstance(..., Iterable)`
guard, and qualifier values are usually *lists* of strings — that is what
every `SeqIO` parser produces, the SnapGene parser included — for which that
expression is an element-membership test, so the shape the motivating issue
actually hits was never matched; bare strings and lists or tuples of strings
are now handled explicitly. Second, the adopted code assigned the cleaned
values back into `feature.qualifiers`, so writing a record permanently
changed the caller's `SeqRecord`; the cleaned values now go to the file only,
via a shallow-copied feature. Third, replacing only `\n` left `\r` behind —
and a bare carriage return read back in text mode is a line break under
universal newlines, regenerating the defect — so CRLF collapses to one space
and a lone CR is treated like LF. The regression test covers all of this and
re-parses the writer's output. The EMBL writer, which shares the defect, was
outside the pull request's scope and left unchanged; the cleaning has since
moved to the shared `_InsdcWriter` base, so both writers do it.

### [#4938](https://github.com/biopython/biopython/pull/4938) — Will Tyler

`PDBList.retrieve_pdb_file(..., file_format="mmtf")` built its URL on
`http://mmtf.rcsb.org`, a host the RCSB decommissioned when it retired the
MMTF format in 2024 — while BinaryCIF, the format the RCSB points MMTF users
at and which `Bio.PDB.binary_cif` already parses, could not be fetched at all.
The pull request removes the mmtf option (requesting it now raises
`ValueError` naming the surviving formats) and adds `file_format="bcif"`,
fetched from `https://models.rcsb.org/`.

**Two of the three commits taken.** The third adds a paragraph to the
release notes upstream was accumulating for Biopython 1.86; in this tree that
part of `NEWS.rst` is preserved upstream history, so the fork's own in-progress
section carries the paragraph instead. There was no review discussion, so
nothing review-grown to leave behind.

Three corrections sit on top as separate commits: the new bcif branch guarded
obsolete requests with a bare `assert`, which `python -O` strips — it now
raises `ValueError`, per the same §0.8 stance the fork adopted
[#5175](https://github.com/biopython/biopython/pull/5175) for; the pull
request's two format-error tests had their `assertEqual` indented inside the
`assertRaises` block after the raising line, so they never ran, and both
expected strings were stale (running them revealed the supported-format list
now names `bcif`); and a new offline test module,
`Tests/test_PDB_PDBList_offline.py`, pins down the constructed URLs with the
download mocked, which no existing PDBList test could do without the network.

### [#5121](https://github.com/biopython/biopython/pull/5121) — Oz Nova

`MMCIF2Dict._splitline` walked every line character by character through the
full CIF 1.1 tokenizer, though the lines that dominate an mmCIF file — the
atom records — contain no quote and no comment character. A line without
`'`, `"` or `#` is now handed to `str.split`; any line that could need the
real tokenizer still gets it. Measured here, `MMCIF2Dict` reads
`Tests/PDB/6WG6.cif` in 181 ms instead of 556 ms and
`MMCIFParser.get_structure` drops from 755 ms to 364 ms (medians of 15
runs); parse output is identical on all 19 mmCIF fixtures in `Tests/PDB`,
compared at both the dictionary and the built-structure level.

**All three commits taken**: the fast path, the author's benchmark script
(`Scripts/Performance/benchmark_mmcif_parsing.py`), and his `CONTRIB.rst`
line. Nothing left behind. peterjc wrote on 2025-12-25 that it "seems good
to me - a worthwhile speedup for minimal change"; the thread has not moved
since. One recorded caveat: `str.split()` splits on any whitespace, where
the tokenizer recognises only space and tab — a difference visible only on
a line whose unquoted token contains `\v`, `\f` or non-ASCII whitespace,
none of which CIF 1.1's ASCII character set allows.

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
[#45](https://github.com/BioAGI-org/biop-ai/pull/45).

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
[#44](https://github.com/BioAGI-org/biop-ai/pull/44).

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
[#43](https://github.com/BioAGI-org/biop-ai/pull/43).

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
[#42](https://github.com/BioAGI-org/biop-ai/pull/42).

### [#5175](https://github.com/biopython/biopython/pull/5175) — Abdel ATIA

Replaces nine `assert` statements that validate input in production code with
explicit `if`/`raise`, in `Bio/ExPASy/Enzyme.py`, `Bio/ExPASy/ScanProsite.py`,
`Bio/SeqIO/TabIO.py` and `Bio/motifs/alignace.py`. Under `python -O` every one
of those checks vanishes, so a malformed file produces silently wrong output
rather than an error — the same defect `IMPROVEMENTS.md` §0.8 describes across
343 sites. Four commits, all taken, tests included.

Approved by mdehoon on 2026-03-18. peterjc's formal review the same day
requested changes (run black, follow the pull request template); on 2026-04-27
he commented "The changes look good, thank you" and asked for `assertRaises`
tests, which were added on 2026-05-01. The same day he asked "are you using a
coding agent LLM AI tool?", and on 2026-09-01 asked again, noting that any AI
use should at least be declared. That question has not been answered, and the
pull request is still open.
