# Upstream reports

Most of what this fork fixes is also present in Biopython, because most of it
predates the fork. Where that is true the defect is worth reporting upstream,
so that people who use Biopython benefit whether or not they ever hear of this
project.

This file records what has been reported, what is queued, what upstream already
knows and this fork has fixed, and what was checked and deliberately not
reported. It is a working record rather than a plan of
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
- **Human-written.** The report is written by a person in their own words, not
  generated text relayed under a name; the findings behind it may come from AI
  work. See [How this fork engages upstream](#how-this-fork-engages-upstream).

## Reported

| upstream | date | subject | fixed here |
|---|---|---|---|
| [#3771](https://github.com/biopython/biopython/issues/3771) | 2026-07-30 | Use-after-free in the `cpairwise2` sequence conversion; root cause for a crash open since 2021 | [#33](https://github.com/BioAGI-org/biop-ai/pull/33) |
| [#5269](https://github.com/biopython/biopython/issues/5269) | 2026-07-30 | `run_cealign` segfaults on a zero fragment size, short coordinate entries, or tuples | [#39](https://github.com/BioAGI-org/biop-ai/pull/39) |
| [#5272](https://github.com/biopython/biopython/issues/5272) | 2026-08-03 | `PrintedAlignmentParser.feed()` does not bound `offset`, and its length check cannot fire | [#2](https://github.com/BioAGI-org/biop-ai/pull/2) |
| email | 2026-08-03 | Out-of-bounds write in `bcifhelpermodule.c` — to the `Bio/PDB` owners under their `SECURITY.md` | [#3](https://github.com/BioAGI-org/biop-ai/pull/3) |
| [#5252 review](https://github.com/biopython/biopython/pull/5252#issuecomment-5171465860) | 2026-08-03 | Review of their in-flight endian fix: the output array is still allocated little-endian | [#36](https://github.com/BioAGI-org/biop-ai/pull/36) |

## How this fork engages upstream

This fork's own tree welcomes AI-written contributions — that is the reason it
exists, and `ADOPTED.md` is full of them. Reaching *into* Biopython is a
separate act, and the fork chooses to make it on Biopython's terms rather than
its own, because the recipient's comfort is the whole value of the contact.
These are the fork's standing choices for that outbound contact, not rules
imposed on its own work.

- **Outbound writing is a person's own words.** An issue, a comment or an email
  sent to Biopython is written by Dan himself and not signed as AI-authored.
  Biopython has said it prefers discussion text to be human-written (peterjc on
  [#5252](https://github.com/biopython/biopython/pull/5252#issuecomment-5171465860):
  "all written text on the discussions ... should be human"), and meeting that
  preference is a courtesy the contact depends on. The *findings* behind the
  writing may still come from AI work — upstream said in the same breath that it
  expects bug reports found with AI tools — but rewriting rather than
  transcribing is the line: generated prose under a human name would keep the
  letter and lose the point.
- **Code goes to upstream as a description, not a patch.** A report says what a
  fix must guarantee and leaves upstream to write it; where a patch already
  exists here it is offered for a clean-room reimplementation rather than
  submitted. This is the deliberate inverse of how the fork treats the same
  work internally, where the patch is adopted directly with its author's name
  on it (`ADOPTED.md`). Upstream has not decided whether it accepts
  AI-written code: the no-AI draft
  [#5241](https://github.com/biopython/biopython/pull/5241) and the permissive
  draft [#5388](https://github.com/biopython/biopython/pull/5388) are both
  open. What it has decided, in
  [#5334](https://github.com/biopython/biopython/pull/5334) (merged
  2026-11-06), is that a pull request must declare any AI use and that all
  communications must be written by a human, machine translation into English
  excepted. Until a policy is agreed, only AI-free pull requests are being
  merged.

## Queued

Verified unreported in the GitHub tracker and still present on upstream
`master`. Ordered by how self-evident the defect is, not by severity.

| # | subject | note |
|---|---|---|
| 1 | `Record.__str__` corrupts a string-valued `wgs` into `A-B-C-D-...` | New in 1.88, and the incomplete half of [#5274](https://github.com/biopython/biopython/issues/5274). [#5275](https://github.com/biopython/biopython/pull/5275) fixed the parsed-list path but left `__init__` at `self.wgs = ""`, so a string assignment now serialises one character at a time. `wgs_scafld` was corrected to `[]` in the same commit. Comment on #5275. **Fixed here in [#62](https://github.com/BioAGI-org/biop-ai/pull/62)**, so the report describes a defect this tree no longer has. |
| 2 | `RestrictionBatch.format` docstring still promises expression evaluation | Documentation only. The 1.88 `eval` removal changed the behaviour and left the docstring saying "If y can be evaluated to a RestrictionType". **Fixed here in [#62](https://github.com/BioAGI-org/biop-ai/pull/62).** |
| 3 | `Bio.File.as_handle` re-catches a `TypeError` from the caller's own block | Fixed upstream in 2018 for [#1544](https://github.com/biopython/biopython/issues/1544), reverted by `f09b169a4` in 2020, no regression test. Comment on #1544 rather than a new issue. |
| 4 | `_get_pi` builds `forward_table.keys() + stop_codons` | `TypeError`; the `F1x4` and `F61` codon frequency models are unusable. |
| 5 | `Alignment.format()` reports any `AttributeError` as "not yet implemented" | Upstream's own `Bio.Align.write()` already does this correctly. |
| 6 | `CodonAdaptationIndex` rejects `MutableSeq` | The docstring promises it is accepted. |
| 7 | GAF by-protein iterators never yield the last protein | Silent data loss. Fixing it changes output for every caller. |
| 8 | `Seq.search()` never matches at the final position | Off by one since the method shipped. |
| 9 | `six_frame_translations` misaligns the reverse frames | Correct only when `len(seq) % 3 == 1`. Display only. |
| 10 | `_G_test` divides by zero on any empty cell | `mktest` crashes on identical sequences. |
| 11 | *(Filed upstream by a third party; see [Checked and not reported](#checked-and-not-reported).)* | |
| 12 | `Tm_GC(valueset=0)` raises `UnboundLocalError` | The guard admits a value no branch handles. |
| 13 | `str()` on an `index_db` dictionary raises `AttributeError` | Upstream carries a `TODO` about this exact question; lead with it. |
| 14 | `SeqIO.index()` raises `IndexError` on a bare `>` FASTA header | `parse()` reads the same record as `id=""`. |
| 15 | Two `except ValueError` clauses in `SeqIO/_index.py` cannot fire | File the `tab` and `genbank` halves separately; cite [#1344](https://github.com/biopython/biopython/issues/1344). |
| 16 | `reverse_complement()` cannot sort features with an `UnknownPosition` | Cite [#1772](https://github.com/biopython/biopython/issues/1772), same cause, different path. |
| 17 | `ProteinAnalysis.flexibility()` reads the wrong window centre | Comment on [#4170](https://github.com/biopython/biopython/pull/4170), which fixes the adjacent window-count bug only. |
| 18 | `neighbor_simple_search` re-sorts the KD-tree's shared point list | `DataPoint_sort` at `kdtrees.c:1226` reorders the list the built tree's node ranges index into, so a later `search()` on the same `KDTree` silently misses points — measured 623 → 397 hits on 2,000 points. Silent wrong answer. **Fixed here in [#74](https://github.com/BioAGI-org/biop-ai/pull/74)**, which sorts a private copy so the tree is immutable after construction; the fork's fix is bundled with a larger rewrite, so describe the defect, not the patch. |
| 19 | `SeqIO.index()` on pretty-printed UniProt XML is silently empty | The `b"<entry "` marker at `_index.py:198` is matched at line start only; an indented `<entry>` never matches, `parse()` reads the same file fine. `Tests/SwissProt/multiligand.xml` in both trees demonstrates. **Fixed here in [#88](https://github.com/BioAGI-org/biop-ai/pull/88)** as a side effect of moving id rules onto the parsers — the narrow fix upstream is a `lstrip` or regex allowance. |
| 20 | FOGSAA sizes its matrix in `int` arithmetic | `PyMem_Calloc((nA+1) * (nB+1), ...)` at `_pairwisealigner.c:6075` and the matching `MATRIX()` index macro wrap for two sequences ≳46,341 nt, under-allocate, and corrupt the heap — segfault from pure Python, reproduced in this tree pre-fix. **Fixed here in [#67](https://github.com/BioAGI-org/biop-ai/pull/67)** with checked `size_t` sizing. Caller-driven, not file-driven, so the public tracker is the right venue (same reasoning as #5272). |
| 21 | [#3911](https://github.com/biopython/biopython/pull/3911)'s newline check never matches list-valued qualifiers | `"\n" in feature.qualifiers[k]` is element membership on a list, so the list-of-strings shape the motivating SnapGene case produces is untouched by the fix. Comment on #3911 rather than a new issue. Fork adopted the PR and corrected this in [#80](https://github.com/BioAGI-org/biop-ai/pull/80). |
| 22 | `System.stochiometry()` has never worked | Three defects in one method, all verbatim from the 2001 import: `[] * n` is `[]`, rows created as the integer `0`, and coefficients written with no row index. Any non-empty system raises `IndexError`; nothing in either tree calls it, which is how it survived 24 years. **Fixed here in [#97](https://github.com/BioAGI-org/biop-ai/pull/97)**; found by mypy's `check_untyped_defs` during [#96](https://github.com/BioAGI-org/biop-ai/pull/96). |
| 23 | `complement`, `complement_rna` and `reverse_complement_rna` promise `SeqRecord` support and raise `AttributeError` | Since [#3648](https://github.com/biopython/biopython/pull/3648) (1.80) all four module functions in `Bio/Seq.py` have a `SeqRecord` branch and a docstring saying "returns a new SeqRecord object", but `SeqRecord` only has `reverse_complement`, so the other three call methods that never existed. Still on `master` `372c71069`. Whether to implement `SeqRecord` output or reject it is upstream's call, so report both; a clear `TypeError` now does not block making it work later. **Fixed here in [#147](https://github.com/BioAGI-org/biop-ai/pull/147)** with the `TypeError`. |
| 24 | `Bio.Align.parse` cannot read stdin or a pipe | Iterating calls `rewind()`, which seeks to 0 (`AlignmentsAbstractBaseClass.__iter__` in `Bio/Align/__init__.py`, and `Bio/Align/interfaces.py:197-198`). `list()` calls `__len__`, which reads the stream and then seeks back (`interfaces.py:89-106`). The fasta, a2m and msf readers also call `tell()` (`fasta.py:75`, `a2m.py:91`, `msf.py:33`). So `list(Align.parse(f, fmt))` on a pipe raises `io.UnsupportedOperation`, where `AlignIO.parse` reads it. On Windows a pipe reports `seekable()` and `seek()` does nothing ([python/cpython#86768](https://github.com/python/cpython/issues/86768), fixed only on CPython `main` in August 2026), so there `list()` can silently return nothing. That last part was seen with a wrapper that mimics the Windows pipe, not on Windows itself. **Fixed here in [#182](https://github.com/BioAGI-org/biop-ai/pull/182)** (`IMPROVEMENTS.md` §1.1). |
| 25 | `Bio.Align`'s Stockholm reader silently misreads interleaved files | It makes one row per sequence line (`Bio/Align/stockholm.py:357-384`). A file with two blocks of equal width reads as 4×5 instead of 2×10, with no error. Other block layouts raise (`Tests/Stockholm/simple.sth`), and duplicate names are accepted as extra rows. `AlignIO` reads these files. The docstring says wrap-around is not supported, so the defect to report is the silent misread, not the missing feature. Background on strictness: [#3747](https://github.com/biopython/biopython/issues/3747). **Fixed here in [#197](https://github.com/BioAGI-org/biop-ai/pull/197).** |
| 26 | `Bio.Align`'s PHYLIP reader reads one alignment per file | It reads the rest of the stream as one alignment (`Bio/Align/phylip.py:134-148`). So a file holding several alignments, such as seqboot output or what `Align.write` makes from several alignments, raises `ValueError` ("Found 5 records in this alignment, told to expect 2"). `AlignIO` reads both. **Fixed here in [#184](https://github.com/BioAGI-org/biop-ai/pull/184).** |
| 27 | `AlignIO`'s MAF reader turns `.` into the reference letter, which contradicts the `s` line's size | `Bio/AlignIO/MafIO.py:181-194` replaces each `.` with the first sequence's letter. This is deliberate (`9470216c4`, 2011) and pinned by `test_reading_alignments_maf1`. In `Tests/MAF/humor.maf`, mm3's `s` line says size 3,424, but `AlignIO`'s row has 5,166 residues. `Bio.Align` reads the `.` as gaps and gets 3,424. A behaviour choice: describe both readings and leave the choice upstream. Not changed here. |
| 28 | `Bio.Align`'s Stockholm reader shifts per-residue `#=GR` annotations | It deletes every `.` from a `#=GR` line (`Bio/Align/stockholm.py:295`) instead of the row's gap columns. So an annotation with a `.` at a residue moves by one place up to the next gap: in `Tests/Stockholm/example_nonstandardannotations.sth`, CRN_DROME's `nonstandardgr` gets `X` at the residue the file marks `.`. When the letters no longer fit the residues it raises `TypeError`. That includes `Bio.Align`'s own output, because its writer (`_format_record`, `stockholm.py:598`) puts the annotation's letters at the residues and `.` at the gaps. **Fixed here in [#183](https://github.com/BioAGI-org/biop-ai/pull/183).** |
| 29 | `AlignIO`'s PHYLIP readers ignore the header's length | `PhylipIterator` parses it (`Bio/AlignIO/PhylipIO.py:206`) and never uses it, so `phylip` and `phylip-relaxed` can misread silently. `AlignIO.read(StringIO("2 10\ns1 AAAAAAAAAA\ns2 CCCCCCCCCC\n"), "phylip")` gives 2 rows of 3 residues, with ids `s1 AAAAAAA` and `s2 CCCCCCC`. Names repeated in a later block end up inside the sequence (`AAAAs1AAAA`). `Bio.Align` raises on both. **Not fixed here**: it is an `AlignIO` change, outside §1.1's new scope. |
| 30 | `PairwiseAligner.score()` and `align()` overwrite the caller's sequence array | With a substitution matrix that has an alphabet, such as `scoring="blastp"`, `_map_indices` writes the matrix indices back into the caller's buffer (`indices[i] = index`, `_pairwisealigner.c:7280` on `master` `590285be7`), and the Python layer passes a NumPy array straight through. Given an `int32` array of the letter codes of `HEAGAWGHE`, the first `score()` succeeds and leaves the array as `[7 4 0 6 0 22 6 7 4]`; a second call raises "sequence contains letters not in the alphabet", and so does passing one array as both sequences, on the first call, because it is mapped twice. Reproduced in this tree before the fix; on `master` the code is unchanged, read but not run. The corruption is not memory safety, so the public tracker. The same change closes a contrived out-of-bounds read: a gap function that writes an out-of-range letter into the array after the check makes the Waterman-Smith-Beyer kernel index outside the matrix (an ASan SEGV in this tree before the fix). That is caller-driven too, so it can go in the same report (as with 20). **Fixed here in [#181](https://github.com/BioAGI-org/biop-ai/pull/181)**, which checks and maps private copies of the sequences; describe the defect, not the patch. |
| 31 | Two `test_PDB_PDBList` tests call `unittest.patch`, which does not exist | `test_retrieve_pdb_file_not_existing` and `test_retrieve_pdb_file_bad_url` (`Tests/test_PDB_PDBList.py:166,181` on `master` `fc0c15678`) raise `AttributeError` whenever they run; the module is online-only, so CI never runs them. Two more faults hide behind it: the message they look for is `print`ed to stdout, not stderr, and the `" http://something.wrong "` server makes `http.client` raise `InvalidURL`, which `retrieve_pdb_file` does not catch. Test-only. **Fixed here in [#206](https://github.com/BioAGI-org/biop-ai/pull/206).** |
| 32 | `GfaIO`'s tag-type check is a substring test | `parts[1] not in "AifZJHB"` ([`Bio/SeqIO/GfaIO.py:65`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/SeqIO/GfaIO.py#L65) on `master` `fc0c15678`), so an empty type (`AB::x`) or several letters that appear together in `AifZJHB` (`AB:if:x`) pass with no warning and skip every per-type value check. Both `gfa1` and `gfa2`. Tracker checked 2026-10-14: nothing filed; the check came in [#4598](https://github.com/biopython/biopython/pull/4598). **Fixed here in [#208](https://github.com/BioAGI-org/biop-ai/pull/208).** |
| 33 | `CodonAligner.align()` frees an uninitialised pointer for fewer than 3 nucleotides | `PathGenerator* paths;` ([`Bio/Align/_codonaligner.c:900`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/Align/_codonaligner.c#L900) on `master` `fc0c15678`) is never set before the new `nB < 3` check from [#5369](https://github.com/biopython/biopython/pull/5369) jumps to `exit`, which runs `Py_XDECREF(paths)` (line 1046). So `CodonAligner().align("M", "A")` decrefs stack garbage; gcc `-Wall -O2` warns `'paths' may be used uninitialized`, and a build with `-ftrivial-auto-var-init=pattern` segfaults. Caller-driven, so the public tracker. Tracker checked 2026-10-14: nothing filed. **Fixed here in [#210](https://github.com/BioAGI-org/biop-ai/pull/210)** (`paths = NULL`). **Fixed upstream independently** by the same one-line change in [#5387](https://github.com/biopython/biopython/pull/5387) (2026-10-22), merged here in the 2026-10-27 sync; nothing left to report. |
| 34 | `test_join_MutableSeq_mixed` calls `join` before `assertRaises` sees it, so the check tests nothing | [`Tests/test_Seq_objs.py:964`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Tests/test_Seq_objs.py#L964) on `master` `fc0c15678`. It passes only because the `MutableSeq` result is not callable, and the `TypeError` it claims is wrong: mixed `Seq`/`MutableSeq` items are documented as accepted. Test-only. Tracker checked 2026-10-14: nothing filed. **Fixed here in [#218](https://github.com/BioAGI-org/biop-ai/pull/218).** |
| 35 | `PDBList.retrieve_pdb_file` leaves an empty structure file after a non-gzip download, and lets `InvalidURL` escape | It decompresses outside its `try` ([`Bio/PDB/PDBList.py:346`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/PDB/PDBList.py#L346) on `master` `fc0c15678`), so a captive-portal page raises `BadGzipFile` after `open(final_file, "wb")` has created the file, and the next call returns it as "Structure exists". It catches only `OSError` (line 339), so `http.client.InvalidURL` escapes instead of returning `None` (the cause behind row 31). Cite [#4920](https://github.com/biopython/biopython/issues/4920), whose mmtf traceback is this decompress path. Tracker checked 2026-10-14: nothing filed on either half. **Fixed here in [#214](https://github.com/BioAGI-org/biop-ai/pull/214).** |
| 36 | `Bio.PDB.binary_cif` indexes the StringArray string table with -1, the null marker | [`Bio/PDB/binary_cif.py:155`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/PDB/binary_cif.py#L155) on `master` `fc0c15678`. A BinaryCIF string column stores a null (mmCIF `.` or `?`) as lookup -1. On a column with no strings at all, such as the all-null `_atom_site.pdbx_PDB_ins_code` that 191 of 200 random current entries have (`label_alt_id`: 99 of 200), `unique_strings[lookups]` raises `IndexError`; on a column with some strings it silently substitutes the last string for each null. 196 of 200 randomly chosen current RCSB `.bcif` files fail to parse. Of the 4 that parse, 2 (9WQR, 9ZHI) have those columns encoded as integer zeros rather than strings, and 2 (4AX9, 5HI3) have real altlocs and insertion codes, so no string table is empty and every null takes the last letter: in 4AX9, 475 residues get insertion code `K`, and atoms merge. The three fixtures in `Tests/PDB` contain no null strings, so the suite never saw it. Reproduce: `BinaryCIFParser().get_structure("1ADB", "1ADB.bcif")` on `https://models.rcsb.org/1ADB.bcif`. Tracker checked 2026-10-16 for "BinaryCIF" and "binary_cif": nothing filed ([#5235](https://github.com/biopython/biopython/issues/5235) is the separate byte-order issue). **Fixed here in [#221](https://github.com/BioAGI-org/biop-ai/pull/221).** |
| 37 | `NewickIO` reads a branch length after two spaces or a tab as the clade's name | The edge-length token allows one optional space after the colon, `\:\ ?` ([`Bio/Phylo/NewickIO.py:29`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/Phylo/NewickIO.py#L29) on `master` `fc0c15678`, unchanged at `c911b1210`). After two spaces, a tab, or a line break and indentation, the colon is skipped and the length is read as an unquoted label, with no error: `((a:  0.1,b:\t0.2),c);` gives leaves named `0.1` and `0.2` and no branch lengths, and on an inner clade, `(a,b)0.95:  0.3` gives confidence 0.3 and loses the 0.95. The Newick format allows blanks, tabs and newlines anywhere outside a label or branch length. The one space is the `\ ?` that [#1821](https://github.com/biopython/biopython/pull/1821) added for PAML's `: ` ([#1782](https://github.com/biopython/biopython/issues/1782)). NEXUS input is not affected: `Bio.Nexus` collapses runs of whitespace before the tree is parsed. Tracker checked 2026-10-18: nothing else filed; [#4134](https://github.com/biopython/biopython/issues/4134) (whitespace inside an unquoted label) is a different defect. **Fixed here in [#226](https://github.com/BioAGI-org/biop-ai/pull/226).** |
| 38 | A Gotoh local alignment with score 0 raises `IndexError` on `[0]`, while `len()` and iteration find its alignments | `GOTOH_LOCAL_ALIGN` raises `maximum` only by more than epsilon, so it can stay 0 while cells scoring exactly epsilon are end points; then [`if (maximum == 0) M[0][0].path = DONE;`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L5480) (`master` `f68b89603`) starts the generator exhausted. `len()` counts the end points, and `PathGenerator_reset`, which iteration calls, clears `DONE`. `PairwiseAligner(mode="local", mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1, epsilon=1).align("GAT", "TA")`: `len()` 2, `list()` 2 alignments, `[0]` `IndexError`. Smith-Waterman sets `NONE` there and has no end points at score 0, so it is consistent. Waterman-Smith-Beyer local has the same line ([L5982](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L5982)) and the same defect: with gap score functions `-2 - 0.5 * n` and `epsilon=1`, `align("A", "A")` has `len()` 1 and `align("GAT", "TA")` 2, and `[0]` raises `IndexError` for both. Both reproduced on PyPI 1.88. Tracker checked 2026-10-20 ("PairwiseAligner", "local alignment", "PairwiseAligner IndexError", "PairwiseAligner epsilon", and for Waterman-Smith-Beyer "Waterman-Smith-Beyer", "Waterman", "gap score function"): nothing filed. **Fixed here in [#238](https://github.com/BioAGI-org/biop-ai/pull/238)** for Gotoh local, and in [#242](https://github.com/BioAGI-org/biop-ai/pull/242) for Waterman-Smith-Beyer local. The same `maximum` ([`SELECT_TRACE_GOTOH_LOCAL_ALIGN`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L4283)) also makes the reported score lag `score()` by up to epsilon, and as the first rise unmarks nothing, end points found before it stay: with the scores above, `"ACGT"` against itself gives score 4 and one path, `A` with `A`; `"AACAA"` against `"AAGAA"` gives score 2 (`score()` 3) and six paths. Separately, `GOTOH_LOCAL_SCORE` takes its maximum from Ix and Iy too, and lets a gap scoring 0 go on ([`SELECT_SCORE_LOCAL3`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L4161)), so with a positive gap score `score()` reports a score no local alignment has: `match_score=1, mismatch_score=-1, open_gap_score=-0.5, extend_gap_score=0.5`, `"AAA"` against `"CCC"`: `score()` 0.5, `align()` 0. Both reproduced on PyPI 1.88. Waterman-Smith-Beyer local has the same lag ([L4393](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L4393)). **Fixed here in [#244](https://github.com/BioAGI-org/biop-ai/pull/244)** for Gotoh local, as #239 did for Smith-Waterman. In Waterman-Smith-Beyer local, with `mismatch_score=-1`, gap score functions `-2 - (n - 1)` and `epsilon=1`, `"GAT"` against `"TA"` scores 0 (`score()` 1), `"ACGT"` against itself gives one path, `A` with `A`, and `"AACAA"` against `"AAGAA"` scores 2 (`score()` 3) with six paths. Its `score()` agrees with brute force, as local mode now rejects positive gap scores. Reproduced on PyPI 1.88. **Fixed here in [#245](https://github.com/BioAGI-org/biop-ai/pull/245)** for Waterman-Smith-Beyer local. |
| 39 | `PairwiseAligner` counts past an overflow, and Smith-Waterman's score, count and paths disagree | Two defects in [`Bio/Align/_pairwisealigner.c`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c) on `master` `f68b89603`, both run on PyPI 1.86 and 1.88. `SAFE_ADD` (line 32) checks only the sum for the overflow marker, -1, so an overflowed count added to one that had not gives a wrong count instead of `OverflowError`, in every algorithm: global, all scores 0, `"A" * 37` against `"C" * 37` gives `len()` 8001394856851978725, where the true count is about 2e27. Smith-Waterman (`SELECT_TRACE_SMITH_WATERMAN_HVD`, line 4213) takes its score from every cell, including those a gap ends in, but marks end points only where a diagonal move ends, and unmarks earlier end points only on a rise of more than epsilon. So with a positive gap score `align()` and `score()` report a score no alignment has, the stale end points stay, and `PathGenerator_next_smithwaterman` raises "Unexpected trace" at one left without its diagonal trace (`"AAA"` against `"CAACC"`, gap score 0.5: score 2.5, the first path scores 1, `[1]` raises). With `epsilon=1`, `"CC"` against itself yields three paths scoring 1 for a score of 2. Tracker checked 2026-10-20: nothing filed. Report as two issues; whether local alignments may start with a gap when gaps score positive is upstream's call (see the fork PR). This fork has since chosen to reject positive gap scores in local mode (`ValueError`, see `DEPRECATED.rst`); offer that as one option, not as the fix. **Fixed here in [#239](https://github.com/BioAGI-org/biop-ai/pull/239).** |
| 40 | The EMBL writer has the same stray-line-break defect as [#3885](https://github.com/biopython/biopython/issues/3885), and [#3911](https://github.com/biopython/biopython/pull/3911) does not touch it | `EmblWriter.write_record` passes each feature straight to `_write_feature` ([`Bio/SeqIO/InsdcIO.py:1522`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/SeqIO/InsdcIO.py#L1522) on `master` `fc0c15678`), so a qualifier value containing a newline, such as the SnapGene parser passes through, is written across two lines with the continuation unindented, and `SeqIO.read` of the writer's own output raises `ValueError: Problem with 'CDS' feature`. A record with `qualifiers={"note": ["line one\nline two"]}` written as `"embl"` or `"imgt"` reproduces it. #3911 fixes `GenBankWriter.write_record` only; the fix belongs on the shared `_InsdcWriter` base. Tracker checked 2026-10-22 ("EMBL line break", "EMBL writer newline qualifier", "EmblWriter", "qualifier newline"): nothing filed beyond #3885 itself, which names GenBank only. Comment on #3885 rather than a new issue. **Fixed here in [#246](https://github.com/BioAGI-org/biop-ai/pull/246)** by moving the cleaning to `_InsdcWriter`, so both writers share it. |
| 41 | The Tutorial's codon-alignment section shows the HIV rows under the *E. coli* heading, checks `len(alignments)` after building `protein_alignments`, cites the wrong mnemonics table, has two typos, and labels LWL85's dN/dS numbers as `YN00` | Documentation only, all from [#5376](https://github.com/biopython/biopython/pull/5376) (2026-10-16). The copied *E. coli* block has an extra blank line after `.. code:: pycon`, so `Tests/test_Tutorial.py` never extracts it and upstream's CI never ran it; the `len(alignments)` line passes only because the HIV iterator is still in scope. Still on `master` `f10e88c58` on 2026-10-28; tracker searched for `alignment_ecoli`, `codonalign-mnemonics` and a `None` from `calculate_dn_ds`, nothing. Detail in [#269](https://github.com/BioAGI-org/biop-ai/issues/269), which mentions #5376 so the cross-reference shows on that PR's timeline. **Fixed here in [#270](https://github.com/BioAGI-org/biop-ai/pull/270)**, together with five older slips in the same section (`calculate_dN_dS`, "and input", `calculate_dn_ds_matrix` "only takes `codon_table`", "these two sequences", "Mcdonald"). |
| 42 | `calculate_dn_ds(method="YN00")` returns `None` when its iteration does not converge, and two defects make it fail on the Tutorial's own alignment | On the 16-codon *E. coli*/*Shigella* alignment: (1) the fold-class normalisation `for i, j in zip(fold0_cnt, fold4_cnt)` stops at the shorter dict, so with no `A` among the 4-fold sites `fold0_cnt["A"]` stays a raw count and the initial kappa is -35.0 (0-fold estimate -42.4); (2) with that fixed (kappa 0.98) the three-difference pair TAT/ATC has every substitution path through a codon with zero observed count, `sum(path_prob)` is 0 at `analysis.py:764`, numpy warns `invalid value encountered in scalar divide`, and every pass is NaN; (3) the convergence test never passes and `_yn00` falls off the end of its 20-pass loop, returning `None` with no error or warning of its own. The same `path_prob` warning is quoted in [#1611](https://github.com/biopython/biopython/issues/1611); see also the two `YN00` lines in the October list below. One report once checked against PAML's `yn00.c`. Detail in [#269](https://github.com/BioAGI-org/biop-ai/issues/269). **Fixed here in [#290](https://github.com/BioAGI-org/biop-ai/pull/290)**, a port of PAML's `yn00.c`: the Tutorial pair now gives `yn00`'s dN 0.0785 and dS 0.3384, in either sequence order. Report as a comment on #1611 (see [Known upstream, fixed here](#known-upstream-fixed-here)). |
| 43 | `WellRecord.__str__` prints the time as the signal and the signal as the time | `__iter__` yields `(time, signal)`, and the minimum and maximum lines format that pair as `"signal %.2f at time %.2f"` ([`Bio/phenotype/phen_micro.py:835-836`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/phenotype/phen_micro.py#L835) on `master` `1ef21971a`). A well with signal 29 at time 0.25 prints "Minimum signal 0.25 at time 29.00"; the class docstring's own example shows the swap. Run here on 2026-11-07; nothing in the tracker. Found by CodeRabbit reviewing #271. Not fixed here yet. |
| 44 | `codonalign.build()` drops a row's codon positions when it pads the row against another row's frameshift insertion | `_align_shift_recs` builds `seq = CodonSeq(rf_table=recs[j].seq.rf_table)` and then `seq += ...` ([`Bio/codonalign/__init__.py:811-812`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/codonalign/__init__.py#L811), and again at L825-827, on `master` `1ef21971a`). `Seq.__add__` returns a new `CodonSeq` with the default `rf_table`, so a row with a frameshift deletion, padded against a frameshift insertion in another row, loses its codon positions, and `get_dn_ds_matrix` then raises "Sequence length is not a multiple of three" for every method. Run in this tree before the fix; `master` read, not run. Tracker checked 2026-11-10 ("rf_table", "codonalign build"): nothing; [#4999](https://github.com/biopython/biopython/issues/4999) is a different cause. **Fixed here in [#300](https://github.com/BioAGI-org/biop-ai/pull/300).** |
| 45 | `QCPSuperimposer` fits float32 coordinates in single precision, and from about 3,000 atoms returns an all-zero rotation | `set()` stores the arrays as given and `run()` only copies them ([`Bio/PDB/qcprot.py:294-331`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/PDB/qcprot.py#L294) on `master` `1ef21971a`), while `set_atoms()` converts to float64 (L275-276). `Atom.coord` is float32, so an array built from it is fitted in float32, the cofactors pass the float32 maximum, and the quaternion comes back zero. For all 6,315 atoms of 2XHE the fit leaves about 31 Å while about 1.7 Å is reported. Run in this tree before the fix; `master` read, not run. Tracker checked 2026-11-10 ("QCPSuperimposer", "qcprot", "QCPSuperimposer float32"): nothing on precision; [#5346](https://github.com/biopython/biopython/issues/5346) is the separate 180° defect. **Fixed here in [#280](https://github.com/BioAGI-org/biop-ai/pull/280)**, which converts to float64 in `run()` and `qcp()`. |
| 46 | `majority_consensus` stops at n − 2 clades whenever the consensus root has three children, so it drops clades every rooted tree has | The early exit at [`Bio/Phylo/Consensus.py:372`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Phylo/Consensus.py#L372) (`master` `1ef21971a`) was added on purpose in 2013 and is right for unrooted trees drawn with three clades at the root, as `nj` draws them. For rooted trees it is not: three copies of `(((A,B),C),(D,E))` give `(D,E,(C,(A,B)))`. Run on 1.88 and upstream `master` during [#283](https://github.com/BioAGI-org/biop-ai/pull/283). Tracker checked 2026-11-10 ("majority_consensus", "consensus"): nothing; the unmerged [#2956](https://github.com/biopython/biopython/pull/2956) removes the exit as part of a rewrite. **Fixed here in [#283](https://github.com/BioAGI-org/biop-ai/pull/283)**, which keeps the exit only when every input tree has three clades at its root. |
| 47 | The v3 CEL reader builds its DatHeader subfield-count error and never raises it | `ValueError("Subfields length isn't 12")` has no `raise` ([`Bio/Affy/CelFile.py:413`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Affy/CelFile.py#L413) on `master` `1ef21971a`), since [#3931](https://github.com/biopython/biopython/pull/3931) turned the assert into it in 2022. The reader then takes subfields by fixed index, so 7 or 11 subfields give a bare `IndexError`, an extra subfield in the middle shifts the fields into an `int()` `ValueError`, and an extra one at the end is accepted silently. Run in this tree before the fix; `master` read, not run. Tracker checked 2026-11-10 ("DatHeader", "CelFile subfields", "Affy CEL IndexError"): [#4144](https://github.com/biopython/biopython/issues/4144), [#3795](https://github.com/biopython/biopython/issues/3795) and [#4635](https://github.com/biopython/biopython/pull/4635) are other DatHeader defects. **Fixed here in [#318](https://github.com/BioAGI-org/biop-ai/pull/318).** |
| 48 | `Entity.center_of_mass()` raises `IndexError` for an entity with children but no atoms | The walk calls `entities.popleft()` on an empty deque ([`Bio/PDB/Entity.py:366-367`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/PDB/Entity.py#L366) on `master` `1ef21971a`), for example on a `Structure` whose only `Model` is empty, or a `Chain` whose one `Residue` has no atoms; an entity with no children already gets a `ValueError`. Copilot's review of [#5279](https://github.com/biopython/biopython/pull/5279) pointed this out, and nothing came of it; there is no issue. Run in this tree before the fix; `master` read, not run. Tracker checked 2026-11-10 ("center_of_mass", "radius of gyration"). Comment on #5279, or a short issue of its own, since the defect is in existing code. **Fixed here in [#306](https://github.com/BioAGI-org/biop-ai/pull/306)** (`ValueError`). |

Rows 24–29 come from the §1.1 re-scope in `IMPROVEMENTS.md`. They were checked
on 2026-10-06 against upstream `master` `372c71069`, and again on 2026-10-08
against `bb295f15a`. The code each row describes has not changed upstream, and
it behaved the same in this fork before rows 24, 25, 26 and 28 were fixed here
(rows 27 and 29 still behave the same). The GitHub tracker has nothing filed for any of
them. The nearest are
[#4147](https://github.com/biopython/biopython/issues/4147) (seeking a SeqIO
iterator) and [#358](https://github.com/biopython/biopython/issues/358) (`.`
against `-` in Stockholm), and neither is the same defect. Line numbers are
upstream's, at `bb295f15a`, because the fixes here moved or removed these
lines. This fork's differ in `Bio/Align/__init__.py`, which row 24
cites by name, and in `stockholm.py` and `phylip.py`, which PR #155 edited.

### From the October 2026 sprint, not yet duplicate-checked

Found by reviewers during the sprint and listed in `IMPROVEMENTS.md` under
"Follow-ups from the October 2026 sprint", where each entry has its file,
reproduction and PR. Each defect is still present on upstream `master`
(`fc0c15678`, checked with `git grep`, not run there); a few, such as the
`mktest` and `_get_pi` lines, are already fixed in this fork. Unless a line
says otherwise, none has been searched for in the tracker yet, so each needs
that check before it joins the table above. The C crashes are reached through
a caller's arguments, not a downloaded file, so this draft sends them to the
public tracker, as with row 30. **Dan to decide:** PR #196's reviewer
recommended the held route for the two `cluster.c` heap overflows, and `pca`
on data read by `Bio.Cluster.read` from a downloaded file is arguably
file-driven.

- `Bio.Cluster`: `clusterdistance(..., method="m")` with an index list longer than the matrix has rows (columns, with `transpose`) overflows the heap.
- `Bio.Cluster`: `pca` on data holding `inf`, or with a column sum that overflows, makes `svd` read `w[-1]`. So does finite data whose column sums are finite but which overflow inside `svd`, such as `[[1e308, 1], [-1e308, 2], [0, 3]]`, so checking the input alone does not fix it.
- `Bio.Cluster`: the C `pca` never checks the eigenvalues array's length, and its shape errors print the wrong array.
- `Bio.Cluster`: `PyTree_length` returns `int` but is installed as a `lenfunc`.
- `Bio.Cluster`: `somcluster` given a plain list raises `AttributeError`; `_cluster.Tree` and `_cluster.pca` ignore keyword arguments.
- `PairwiseAligner`: `del aligner.mode` or `del aligner.wildcard` segfaults; no setter checks for `NULL`.
- `CodonAligner`: `del aligner.wildcard` segfaults (the score setters raise `ValueError`/`TypeError`), and `align()` decrefs its `PathGenerator` twice when a trace-row allocation fails ([`_codonaligner.c:950`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/Align/_codonaligner.c#L950) then L1046; ASan reports use-after-free under `_testcapi.set_nomemory`).
- `CodonAligner.score()` and `.align()` test `isinstance(seqA, str)` where they mean `seqB`, so `CodonAligner().score(Seq("MSR"), "ATGTCTCGT")` raises "seqB must be a string, ..." ([`Bio/Align/__init__.py:5182`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/__init__.py#L5182) in `align()`, and L5105 in `score()`, on `master` `f68b89603`; from [#4479](https://github.com/biopython/biopython/pull/4479)). Tracker searched 2026-10-20: nothing filed.
- `PairwiseAligner`: `str(aligner)` overflows a 1,024-byte stack buffer when the substitution matrix's class name is long or the scores are huge (each 1e308 prints 300+ digits; two such scores, or `gap_score = -1e308`, which sets 12, are enough), and uses the name as a format string.
- `PairwiseAligner` maps a list of letters twice when the substitution matrix's alphabet is a string: the Python layer maps each item to its index in the alphabet, and the C code then maps that index as a code point through the same alphabet. So `PairwiseAligner(scoring="blastp").score()` and `.align()` on two lists or tuples of letters raise "sequence contains letters not in the alphabet", as do `blastn` and BLOSUM62; with a large custom alphabet, such as all 95 printable ASCII characters, they can instead return a wrong score without an error. `Alignment.counts()` repeats the Python mapping, so it raises or miscounts the same way. All of this run on upstream `master` (`c911b1210`); the `counts()` code is unchanged on `master` `f68b89603` (read, not run). Tracker checked 2026-10-18: [#4431](https://github.com/biopython/biopython/issues/4431) (open, filed against 1.81, when the error was "sequence has unexpected format") reports this case; it does not mention `counts()`, and a search on 2026-10-20 found no issue that does. There mdehoon suggested letting `Array` keep a list alphabet when the user asks for one; this fork instead makes a list of single letters score, align and count like the string of those letters. That is a behaviour choice, so comment on #4431 rather than open a new issue, and leave the choice upstream.
- `PairwiseAligner` without a substitution matrix ignores the wildcard in a list or tuple of letters, and never matches such a list against a `str`, `Seq`, `bytes` or code-point array: the Python layer numbers list items by their index among the items seen ([`Bio/Align/__init__.py:4511-4518`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/__init__.py#L4511-L4518), and the same in `counts()` at [L3948-L3956](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/__init__.py#L3948-L3956), `master` `f68b89603`), while the C code compares them with the wildcard's code point. `PairwiseAligner(mismatch_score=-1, wildcard="N").score("GANT", "GACT")` is 3.0, 2.0 for lists, and without the wildcard `list("GACT")` against `"GACT"` scores -4.0 (4.0 as strings). With 79 or more distinct items, the 79th is taken for `"N"`. Read on `master`, not run. Tracker checked 2026-10-20 ("PairwiseAligner wildcard", "wildcard list PairwiseAligner"): nothing filed. This fork maps letters to code points (the same choice as for #4431, so mention it there). **Fixed here in [#241](https://github.com/BioAGI-org/biop-ai/pull/241).**
- `PairwiseAligner`: a Gotoh global alignment with infinite or NaN scores, or an infinite epsilon, segfaults when iterated, or loops on one path while `len()` is 0, as `PairwiseAligner(open_gap_score=-2, extend_gap_score=-1, end_deletion_score=float("-inf")).align("AA", "A")[0]` and `epsilon=float("inf")` with `"ACGT", "AGT"` do. Its traces can lead into row 0 or column 0 where no path starts, and unlike `GOTOH_LOCAL_ALIGN` it never removes them ([`_pairwisealigner.c:5291`](https://github.com/biopython/biopython/blob/f68b896031088a8ad11ae49f5e903d49228fd61f/Bio/Align/_pairwisealigner.c#L5291)). Both crash on Biopython 1.88 from PyPI (run); upstream `master` (`f68b89603`) read, not run. Tracker searched 2026-10-20: nothing.
- `Bio/PDB/kdtrees.c`: `repr()` of `Point` and `Neighbor` raises `SystemError`, and `KDTree(coords, bucket_size=10)` ignores the keyword.
- `FastaBlastIterator`'s `self_line = line` typo gives every record the first record's title.
- `Seq.join()` given a generator returns an empty sequence, because the type check uses it up.
- The FASTA, QUAL, PIR, PHD and ACE index key comes from splitting raw bytes, so an id containing a control or Unicode space gives a key that is not `record.id`.
- The UniProt XML and IntelliGenetics indexers compare line lengths with `tell()`, which fails an `assert` on multi-block BGZF.
- `SffIO`: truncated header fields are not length-checked, a negative name padding reads the rest of the file, and `SffWriter.write_records` can be called twice.
- `SffIO`: a file cut short inside a read raises `struct.error` ("unpack requires a buffer of 16 bytes") rather than `ValueError`, and `get_raw` returns the short read ([`_sff_read_seq_record`](https://github.com/biopython/biopython/blob/fc0c156786e86862112a05a720b6f4f47023f36a/Bio/SeqIO/SffIO.py#L898)).
- `InsdcIO`: a short locus name with a sequence of 10^11 bp or more raises instead of using the long LOCUS line.
- `CifSeqresIterator`'s docstring names `_entity_poly_seq`; the code reads `_pdbx_poly_seq_scheme`.
- `Bio.SeqFeature`: `ExactPosition + float` truncates, one `LocationParserError` lacks its `f` prefix (fixed upstream in [#5393](https://github.com/biopython/biopython/pull/5393), merged here in the 2026-11-07 sync), `Reference() == None` raises, and `CompoundLocation(generator)` raises.
- `Bio/Align/emboss.py`: `raise ValueError("...%s") % line` raises `TypeError` instead.
- `Bio/Align/exonerate.py`: an error message never gets its `%` arguments.
- `Bio/Align/sam.py`: the writer writes topology as `PP:` instead of `TP:`.
- `calculate_dn_ds`: NG86, LWL85 and ML do not check that the rows have the same number of codons (read, not run).
- `calculate_dn_ds(method="YN00")` on identical sequences raises "dictionary changed size during iteration". **Fixed here in [#290](https://github.com/BioAGI-org/biop-ai/pull/290).**
- `_count_diff_YN00` counts one step twice and passes `k[1]` where `k[2]` looks intended; check against PAML's `yn00.c` first. **Fixed here in [#290](https://github.com/BioAGI-org/biop-ai/pull/290)**, checked against `yn00.c`; the same PR found and fixed a third defect, `_count_site_YN00` taking the base frequency of the codon changed to rather than the codon's own base.
- `codonalign.mktest`: `"-" in all_codon` never matches, so a partly gapped codon raises `KeyError: '---'`.
- `_get_pi`'s F1x4 counting raises `KeyError` for a base absent from both sequences (masked by row 4).
- `Bio.Align`'s Stockholm writer raises `AttributeError` for an `Alignment` with no `column_annotations`.
- `Bio.Align`'s Stockholm reader raises `IndexError` on a `#=GR` line with no value.
- `Bio/Align/bigbed.py`: the writer's zoom record count disagrees for a zero-length item at the end of a chromosome (an `assert` upstream).
- `Bio.Align`'s strict PHYLIP writer does not reject names that truncation to 10 characters makes equal.
- `AlignIO`'s `RelaxedPhylipWriter` rejects ASCII whitespace only, so it can write a file it cannot read back.
- `Bio/Align/tabular.py`: BLASTN or BLASTP output with a `query seq` column raises "Unknown program".
- `AlignIO`'s Clustal reader hangs on a file holding only the header line.
- `AlignIO`'s Mauve reader never formats one error message, and shares `_ids` between iterators.
- `AlignIO`'s FASTA reader raises `IndexError` on a bare `>>` line.
- `qblast(format_type="HTML")` can never succeed: the polling loop's `<!DOCTYPE html` branch matches the result page.
- `Bio.Blast`'s XML2 parser fails with unhelpful errors on a root with no namespace or no `schemaLocation`.
- `SearchIO`'s BLAST XML `get_raw` loops forever if the file was truncated after its `index_db` was built.
- `SearchIO.index_db` does not wrap a lone `pathlib.Path` in a list, as it does a `str`.
- `Bio/PDB/binary_cif.py`'s string decoder ignores the column mask: `IndexError`, or wrong values, on current RCSB files.
- `Bio/PDB/binary_cif.py` and `FragmentMapper` do not accept `os.PathLike`.
- `FastMMCIFParser`, `PDBMLParser` and the BinaryCIF parser do not upper-case `type_symbol`, so `Fe` raises.
- `StructureAlignment(..., aligner=...)` raises `UnboundLocalError`; the Tutorial shows this call.
- `PICIO._default_dihedron` has two unreachable branches that would crash (`list.add`).
- The `AtomKey` docstring lists a tuple form as accepted that raises.
- `MissingAtomError` is advertised and never raised.
- `Phylo.read(StringIO(...), "cdao")` raises on `handle.name`, and the RDFlib version error is swallowed.
- `Medline.parse`: records separated only by a line of six or more spaces merge (the trade-off of #4557/#4558).
- `Bio.SCOP`: a CLA sunid missing from the DES file gives a bare `KeyError`.
- `GenePop.FileParser`: `FileRecord.start_read()` again leaks the handle and doubles the loci; `remove_loci_by_position()` mutates the caller's list.
- `GenomeDiagram`: `Diagram.write` rejects a `pathlib.Path` for PDF output (ReportLab); SVG works.
- `GenomeDiagram`: `new_graph` silently accepts `centre=`, though `DEPRECATED.rst` says it went in 1.62; `new_track` and `new_set` silently accept misspelt keywords.
- `SVDSuperimposer.set` with 1-D arrays raises `IndexError` instead of its own error.

## Needs work before it can be filed

- **`translate()` drops an explicit `gap=` for `Seq` input.** Open PR
  [#4994](https://github.com/biopython/biopython/pull/4994) changes the
  module-level default and does not forward the argument, so this survives it.
  Any report must say so or it will be closed as a duplicate.
- **Infernal `*[ 0]*` hit lines.** `SearchIO`'s Infernal text parser's
  `assert gap_len > 0` probably rejects valid `*[ 0]*` output. Needs a real
  Infernal file showing it first.
- **`SeqIO.index()` masks a proxy's `TypeError`.** The only reproduction is a
  patched SFF header. Needs a natural malformed input first.
- **Location truthiness.** The dead Python 2 `__nonzero__` methods and the
  stale docstring are real, but the resolution is a behaviour choice, and
  either direction touches `Bio/SeqIO/SnapGeneIO.py:180`.
- **`cpairwise2.rint` parses `int` with `"l"`.** A genuine LP64 stack
  overwrite, but reachable only from the deprecated `Bio.pairwise2`.

- **TwoBit `__getitem__` mis-slices three ways.** Reproduced here and on
  stock Biopython 1.85: `_twoBitIO.c` truncates `(end - start) / step`
  instead of ceiling-dividing, so extended slices silently drop the final
  base (`seq[0::3]` loses the last codon position unless `len % 3 == 0`);
  negative-step slices compute a negative byte count and raise
  `RuntimeError` (`seq[::-1]`); and an integer index equal to the record
  length reads a base from the *next* record's packed data instead of
  raising `IndexError`. Needs a duplicate-check against the tracker and
  re-verification on `master` before filing.
- **`BgzfWriter` crashes on poorly-compressible input.** ~64 KB of
  high-entropy data deflates past the 65,536-byte block cap and raises
  `RuntimeError: TODO - Didn't compress enough`; compressed sizes
  65,511–65,536 slip the guard and fail later in `struct.pack("<H", ...)`.
  htslib caps per-block input instead. Fixed here in
  [#126](https://github.com/BioAGI-org/biop-ai/pull/126) (2026-10-06) by
  capping blocks at 0xFF00 bytes. Fixed upstream independently in
  [#5404](https://github.com/biopython/biopython/pull/5404) (2026-11-06, for
  [#5403](https://github.com/biopython/biopython/issues/5403), filed by a third
  party), which instead splits an oversized block in half after compressing it.
  **That leaves a defect worth reporting:** a virtual offset `tell()` returned
  before the split no longer seeks. With 65,536 random bytes written as 40,000
  then the rest, `tell()` gives `(0, 40000)` and `BgzfReader.seek()` then
  raises "Within offset 40000 but block size only 32768" on upstream `master`
  `1ef21971a` (run 2026-11-07); the fork reads the right bytes. The 2026-11-07
  sync kept the fork's fix and adapted #5404's test to it. Comment on #5403.

### Found while fixing known upstream defects, 9 and 10 November 2026

Found by the pull requests in
[Known upstream, fixed here](#known-upstream-fixed-here). Each is still on
upstream `master` `1ef21971a` (read, not run there, unless the line says
otherwise) and **is not fixed here either**, so each needs a fix, or at least
a reproduction, before it is worth filing. The tracker was searched on
2026-11-10 with the terms given, and had nothing on any of them.

- **LWL85 depends on the order of the two sequences.** `_diff_codon` classes
  each differing position by the degeneracy of the first sequence's codon only
  (`fold_dict[codon1]`,
  [`Bio/Align/analysis.py:367`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Align/analysis.py#L367));
  PAML's `yn00` counts half in each codon's class. On the EGFR test pair
  [#289](https://github.com/BioAGI-org/biop-ai/pull/289) got dN 0.0211 and dS 0.0170 one way round and 0.0199 and 0.0216 the
  other; averaging the two difference counts reproduces `yn00` exactly. Run
  here. Searched "LWL85": only [#5392](https://github.com/biopython/biopython/pull/5392) and [#1611](https://github.com/biopython/biopython/issues/1611).
- **NG86 and YN00 raise a bare `KeyError` on a stop codon.** A terminal `TAA`
  with table 1 gives `KeyError('TAA')` instead of a clear `ValueError`; `_count_site_NG86`
  looks the codon up in `forward_table`
  ([`analysis.py:175`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Align/analysis.py#L175)).
  Run on this fork's `main` during [#290](https://github.com/BioAGI-org/biop-ai/pull/290). Searched "stop codon
  calculate_dn_ds", "calculate_dn_ds stop codon KeyError".
- **`Primer3.parse` loops forever on empty input.** The loop that skips
  leading blank lines never checks for end of file: `readline()` keeps
  returning `""`
  ([`Bio/Emboss/Primer3.py:98-101`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Emboss/Primer3.py#L98)).
  Seen on 1.88 during [#319](https://github.com/BioAGI-org/biop-ai/pull/319). Searched "Primer3", "Primer3 parse empty",
  "Primer3 infinite loop".
- **`Bio.UniGene` sets `homol = True` for `HOMOL NO`.** Both branches assign
  `True`
  ([`Bio/UniGene/__init__.py:299-302`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/UniGene/__init__.py#L299)),
  which looks like a copy-paste slip. Read, not run. Searched "UniGene",
  "UniGene HOMOL": [#2853](https://github.com/biopython/biopython/issues/2853) quotes the code but does not mention this.
- **`BgzfWriter` rejects a `gzip.GzipFile` handle on Python 3.10 to 3.12.**
  There `GzipFile.mode` is an `int`, so `"b" not in getattr(fileobj, "mode",
  "wb")`
  ([`Bio/bgzf.py:805`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/bgzf.py#L805))
  raises `TypeError`. Read, not run. Searched "BgzfWriter", "GzipFile",
  "BgzfWriter GzipFile".
- **`get_support` raises `TypeError` for a tree with an unnamed tip among
  named ones.** It sorts the tip names with plain `sorted()`
  ([`Bio/Phylo/Consensus.py:514`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Phylo/Consensus.py#L514)),
  which cannot compare `None` with a string; [#283](https://github.com/BioAGI-org/biop-ai/pull/283) changed the consensus
  functions but not this one. Read, not run. Searched "get_support",
  "get_support TypeError": [#1887](https://github.com/biopython/biopython/issues/1887) and [#1537](https://github.com/biopython/biopython/issues/1537) are about other things.
- **`DistanceCalculator` with a scoring matrix raises `AttributeError` for a
  `Bio.Align.Alignment` holding a letter outside the matrix.** The intended
  `ValueError` message reads `seq1.id`
  ([`Bio/Phylo/TreeConstruction.py:534`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Phylo/TreeConstruction.py#L534)),
  but `Alignment` rows are plain strings. Seen on 1.88 and this fork's
  `main` during [#292](https://github.com/BioAGI-org/biop-ai/pull/292). Searched "DistanceCalculator", "DistanceCalculator
  AttributeError".
- **`FastMMCIFParser` ignores an indented `_atom_site.` header.** It tests
  `line.startswith("_atom_site.")`
  ([`Bio/PDB/MMCIFParser.py:441`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/PDB/MMCIFParser.py#L441)),
  so the file attached to [#2259](https://github.com/biopython/biopython/issues/2259) gives a `KeyError` there. Seen during [#313](https://github.com/BioAGI-org/biop-ai/pull/313).
  Searched "FastMMCIFParser", "FastMMCIFParser indented".
- **`Bio.Align`'s FASTA and A2M readers raise `IndexError` on a blank line
  before the first record.** A line that is not a header is appended to
  `lines[-1]`
  ([`Bio/Align/fasta.py:73`](https://github.com/biopython/biopython/blob/1ef21971aa828fd357897567044357e9f985abc8/Bio/Align/fasta.py#L73)),
  which does not exist yet. This fork's `Tests/test_Align_streams.py` pins
  the `IndexError` as a known bug; [#323](https://github.com/BioAGI-org/biop-ai/pull/323) fixed the same shape in Clustal and
  MSF only. Read, not run. Searched "Align fasta IndexError".
- **The Newick writer turns an unnamed clade's comment-borne support into a
  name.** `(C:0.02,D:0.02):0.1[99]`, read with `comments_are_confidence=True`,
  is written `(C:0.02,D:0.02)99.00:0.1[99]`, and reading that back the same
  way gives a clade named `99.00`. Seen on 1.88 and upstream `master` during
  [#297](https://github.com/BioAGI-org/biop-ai/pull/297), which fixed only named clades. Searched "comments_are_confidence",
  "Newick comments_are_confidence write": only the unrelated [#360](https://github.com/biopython/biopython/pull/360).

## Draft only, not sent

Material for Dan to rewrite in his own words and send, or not. **Nothing in
this section has been sent: no issue, comment or pull request has been opened
on `biopython/biopython` about it.**

### A format-registration hook and plugin entry points

- **What the fork is adding.** A public
  `Bio.SeqIO.register_format(name, iterator=None, writer=None, *, replace=False)`
  and `Bio.Align.register_format(name, module, *, replace=False)`, plus
  entry-point plugins in the groups `biop_ai.seqio.iterators`,
  `biop_ai.seqio.writers` and `biop_ai.align`. The contract is decided in
  `IMPROVEMENTS.md` §1.2. The private registry it builds on
  ([#149](https://github.com/BioAGI-org/biop-ai/pull/149)), the hook
  ([#185](https://github.com/BioAGI-org/biop-ai/pull/185)) and the plugins
  ([#198](https://github.com/BioAGI-org/biop-ai/pull/198)) have all landed.
- **Why upstream may care.** Upstream's format tables are private. It has no
  public registration hook and no entry points, and no open issue or pull
  request proposes one (searched 2026-10-06, `master` at `372c71069`,
  1.89.dev0). Downstream code assigns into the private
  `SeqIO._FormatToIterator` and `_FormatToWriter` instead: hybran, LMAT,
  recentrifuge and psico do, and tfbayes does the same to AlignIO's table.
- **The question worth asking is the group names.** The fork uses its own
  `biop_ai.*` groups so that it claims nothing in upstream's namespace. If
  upstream wanted a hook, shared `biopython.*` group names would let one plugin
  package serve both projects. The fork could read those groups as well
  without breaking its own plugins.
- **Related upstream threads.** #3614, #3926 and #894 are open; #4627 is
  merged.
  [#3614](https://github.com/biopython/biopython/issues/3614): the "Unknown
  format" error should list the formats.
  [#3926](https://github.com/biopython/biopython/pull/3926): a stalled 2022 WIP
  pull request for public readable and writable format lists. peterjc liked it
  as a single source of truth; mdehoon wanted one `formats`, as `Bio.Align`
  has. [#894](https://github.com/biopython/biopython/issues/894): delay
  SeqIO's imports. [#4627](https://github.com/biopython/biopython/pull/4627)
  (merged 2024) went the same way for `MultipleSeqAlignment`. Upstream's
  `import Bio.SeqIO` still imports every format module, `Bio.AlignIO`,
  `Bio.Align` and NumPy; this fork's
  [#73](https://github.com/BioAGI-org/biop-ai/pull/73) made it lazy.
- **Form.** Information, not a pull request, as in
  [How reports are made](#how-reports-are-made). Any code is offered as a
  description.

## Disclosed, in progress

Both of the memory-safety findings have now been passed on. Neither is fully
fixed upstream.

- **`Bio/PDB/bcifhelpermodule.c` writes past the output buffer.** Sent by email
  to the `Bio/PDB` owners on 2026-08-03 under their `.github/SECURITY.md`,
  rather than to the public tracker, because it is driven by file content. The
  patch is AI-written, so upstream is reimplementing it clean-room from the
  report; Will Tyler has been asked to take it but can't promise time. The CVE
  is going through a GitHub Security Advisory raised by upstream, not MITRE:
  on 2026-08-06 Peter Cock offered the same route he used for the NEXUS
  `eval()` bug (GHSA-cxcq-r753-xf25, still without a CVE weeks later), with
  this fork's maintainer credited as finder. On 2026-09-06 he asked Will and
  Michiel de Hoon to confirm the severity and log it. As of 2026-10-16 there
  has been no reply, no published advisory and no fix on upstream `master`;
  still none on 2026-11-07, and the file is unchanged there since 2025-08-12. The precedent
  before that, CVE-2025-68463, was assigned by MITRE at the reporter's request.
- **`Bio/Align/_aligncore.c` returns heap contents.** Filed publicly as #5272,
  since it needs a caller to pass malformed input rather than arriving in a
  downloaded file. Different module and different owner (@mdehoon) from the
  BinaryCIF one. Closed on 2026-08-31 with
  [#5298](https://github.com/biopython/biopython/pull/5298), a third party's
  hand-written fix for the unbounded `offset` only. The second half of the
  report is still open on `master` `1ef21971a`: the row-length check
  `buffer + m != s` cannot fire, so ragged rows are accepted and
  `parse_printed_alignment()` returns coordinates read from the heap.
  [#5297](https://github.com/biopython/biopython/pull/5297), an earlier
  AI-assisted version that fixed both halves, was closed unmerged.

## Known upstream, fixed here

These defects already have an upstream issue or pull request, so there is
nothing new to report. The fixes were merged here on 9 and 10 November 2026.
Upstream states were read with `gh` on 2026-11-10: each issue cited is still
open there unless the row says otherwise. A comment on the upstream issue
pointing at the fix is a possible follow-up, not a plan. Like any report it
would be written by a person, sent one at a time, and re-checked against
`master` on the day; see [How reports are made](#how-reports-are-made).
Nothing in this section has been sent.

"Adopted" means the upstream commits were cherry-picked with their authors
kept, as recorded in `ADOPTED.md`; "written here" means the fork's own code.
Where a fork correction found a fault in an open upstream pull request, the
fourth column says so; that is review material for the pull request, not a new
issue. Defects these pull requests found that upstream has no report of are
queued as rows 44–48, or listed under
[Needs work before it can be filed](#needs-work-before-it-can-be-filed) when
they are not fixed here either.

| upstream | defect | fixed here | how | still not fixed |
|---|---|---|---|---|
| [#902](https://github.com/biopython/biopython/issues/902), [#4205](https://github.com/biopython/biopython/issues/4205); PR [#905](https://github.com/biopython/biopython/pull/905) (open) | The INSDC writers write a value-less qualifier such as `/pseudo` as `/pseudo=""`, which INSDC forbids, so a RefSeq record with pseudogenes does not round-trip | [#273](https://github.com/BioAGI-org/biop-ai/pull/273) | Written here, on the writer side only. #905 is unlicensed and changes the parser's data shape | A repeated bare qualifier is still dropped by the parser |
| [#4644](https://github.com/biopython/biopython/issues/4644); PR [#5366](https://github.com/biopython/biopython/pull/5366) (open) | `Atom._assign_element` gives the unknown element `X` for four-character hydrogen names such as `HNZ1`, and for unpadded mmCIF names such as `C1` | [#274](https://github.com/BioAGI-org/biop-ai/pull/274) | Adopted from #5366, plus a test of unpadded names | Left-justified hydrogen names such as `"H1  "` still give `X` |
| [#1383](https://github.com/biopython/biopython/issues/1383); PR [#5391](https://github.com/biopython/biopython/pull/5391) (open) | `PDBIO` writes `TER` 81 columns wide, with the serial number the next chain's first atom then reuses | [#275](https://github.com/BioAGI-org/biop-ai/pull/275) | Adopted from #5391. Its `preserve_atom_numbering` docstring says stored serials leave room for each `TER`; they do not when a chain's waters follow its polymer (1LCD, 2XHE, 7DDO), and the fork corrects the text | With `preserve_atom_numbering=True` a `TER` can still share a serial; [#5381](https://github.com/biopython/biopython/pull/5381) (`write_ter`) not taken |
| [#3456](https://github.com/biopython/biopython/issues/3456) | `Bio.Nexus.Trees` strips only the first `[&...]` comment, so a MrBayes consensus tree gives branch lengths as strings such as `'0.1[&l=2]'` | [#276](https://github.com/BioAGI-org/biop-ai/pull/276) | Written here | The reporter's own file still fails in `Bio.Nexus`, because `_parse` counts parentheses inside comments; a labelled internal node with a comment loses its label |
| [#4521](https://github.com/biopython/biopython/issues/4521) | SwissProt CC text wrapped after a hyphen is joined with a space: `GPI- anchor`, `HAMAP- Rule` | [#277](https://github.com/BioAGI-org/biop-ai/pull/277) | Written here; changes comment strings in 62% of Swiss-Prot entries | — |
| [#745](https://github.com/biopython/biopython/issues/745) | Rerooting moves branch lengths but leaves support values on the wrong branches | [#278](https://github.com/BioAGI-org/biop-ai/pull/278) | Written here | PR [#5390](https://github.com/biopython/biopython/pull/5390) (branch lengths inflated when the outgroup is a child of the root) is not fixed here; comments are not moved |
| [#5346](https://github.com/biopython/biopython/issues/5346); PR [#5351](https://github.com/biopython/biopython/pull/5351) (closed) | `QCPSuperimposer` returns a wrong rotation for an exact 180° turn or a mirror image while reporting the right RMSD | [#280](https://github.com/BioAGI-org/biop-ai/pull/280) | Written here, with the diagnosis and the eigenvector check from #5351; its `eigh` fallback is not taken. The same PR fixes a float32 defect that has no report, queued as row 45 | — |
| [#1208](https://github.com/biopython/biopython/issues/1208); PRs [#1213](https://github.com/biopython/biopython/pull/1213) (closed), [#1506](https://github.com/biopython/biopython/pull/1506) (open) | `DistanceTreeConstructor.upgma()` computes WPGMA | [#281](https://github.com/BioAGI-org/biop-ai/pull/281) | Written here to the same maths as #1506, which no longer applies; the old method stays as `wpgma()` | — |
| [#1277](https://github.com/biopython/biopython/issues/1277) | `Tm_NN` charges the initiation terms on a mismatched or dangling end, so one duplex melts at a different Tm depending on which strand is given | [#282](https://github.com/BioAGI-org/biop-ai/pull/282) | Written here | `salt_correction` still takes the whole of `seq`, so methods 5 to 7 still differ between strands |
| [#3345](https://github.com/biopython/biopython/issues/3345); PR [#2956](https://github.com/biopython/biopython/pull/2956) (open) | Strict, majority-rule and bootstrap consensus trees get wrong clades and support values when the trees list their tips in different orders | [#283](https://github.com/BioAGI-org/biop-ai/pull/283) | Written here, with the sorted-name fix the reporter proposed; #2956 is unlicensed. The same PR fixes an early stop that has no report, queued as row 46 | `get_support` still sorts tip names with plain `sorted()`; see Needs work |
| [#5345](https://github.com/biopython/biopython/issues/5345); PR [#5353](https://github.com/biopython/biopython/pull/5353) (open) | `CodonAdaptationIndex.calculate()` always leaves out ATG and TGG and skips the standard stop codons, whatever genetic code the index was built with | [#284](https://github.com/BioAGI-org/biop-ai/pull/284) | Written here, with the same single-codon rule as #5353, which also changes illegal-codon and stop-codon behaviour | The `ZeroDivisionError` when nothing is left to count, [#5306](https://github.com/biopython/biopython/issues/5306) and PR [#5311](https://github.com/biopython/biopython/pull/5311) (open) |
| [#2112](https://github.com/biopython/biopython/issues/2112) | GenBank and EMBL join a wrapped `/anticodon` or `/transl_except` with a space, and the writers break these values at column 80 anywhere | [#285](https://github.com/BioAGI-org/biop-ai/pull/285) | Written here | — |
| [#5113](https://github.com/biopython/biopython/issues/5113) | `internal_coords` string keys such as `"N:CA:C:1N"` return `None` when the neighbouring residue's atom has altlocs or an occupancy other than 1 | [#286](https://github.com/BioAGI-org/biop-ai/pull/286) | Written here | — |
| [#5368](https://github.com/biopython/biopython/issues/5368) (closed upstream as solved by #5317, which it is not); PRs [#5405](https://github.com/biopython/biopython/pull/5405) (closed), [#5373](https://github.com/biopython/biopython/pull/5373) (closed), [#5317](https://github.com/biopython/biopython/pull/5317) (merged) | `QCPSuperimposer` reports a wrong RMSD, and prints "did not converge", for points exactly on a line | [#287](https://github.com/BioAGI-org/biop-ai/pull/287) | Adopted from #5405, and the test from #5373; the RMSD is now measured on the fit on every path, which also corrects sets spread over a few thousandths of an Å. #5405 as upstream has it gives silent wrong rotations without #5346's check. A comment on #5368 should say all of this; the defect is still on `master` `1ef21971a` | Points close to a line but not on it can still get a slightly worse fit |
| [#5312](https://github.com/biopython/biopython/issues/5312); PRs [#5356](https://github.com/biopython/biopython/pull/5356) (closed), [#5324](https://github.com/biopython/biopython/pull/5324) (closed) | `IsoelectricPoint.pi()` searches only pH 4.05 to 12 and silently returns an end of that range | [#288](https://github.com/BioAGI-org/biop-ai/pull/288) | Adopted from #5356, one of its four commits | — |
| PR [#5392](https://github.com/biopython/biopython/pull/5392) (open) | LWL85 weights the twofold-site distances by the fourfold site count | [#289](https://github.com/BioAGI-org/biop-ai/pull/289) | Adopted from #5392 | LWL85 still depends on the order of the two sequences; see Needs work |
| [#1611](https://github.com/biopython/biopython/issues/1611), the YN00 half | `calculate_dn_ds(method="YN00")` returns `None` or raises, and otherwise differs from PAML's `yn00` (queued row 42) | [#290](https://github.com/BioAGI-org/biop-ai/pull/290) | Written here, as a port of PAML 4.10.10's `yn00.c` | The NG86 and LWL85 half of #1611 (saturation); NG86 and YN00 raise a bare `KeyError` on a stop codon, see Needs work |
| [#5340](https://github.com/biopython/biopython/issues/5340); PR [#5363](https://github.com/biopython/biopython/pull/5363) (open) | `ProteinAnalysis.protein_scale` divides by zero for window 1 and double-counts a residue for even windows | [#291](https://github.com/BioAGI-org/biop-ai/pull/291) | Adopted from #5363, plus `ValueError` for a window below 1 and for window 2 with edge 0 | — |
| [#5332](https://github.com/biopython/biopython/issues/5332); PR [#5365](https://github.com/biopython/biopython/pull/5365) (open) | `DistanceCalculator`'s identity model with `skip_letters` divides by the full length, so skipped columns inflate the distance | [#292](https://github.com/BioAGI-org/biop-ai/pull/292) | Adopted from #5365, plus a test with `Bio.Align.Alignment`, which #5365's test never passes | A scoring-matrix model given an `Alignment` with a letter outside the matrix raises `AttributeError`; see Needs work |
| PR [#5389](https://github.com/biopython/biopython/pull/5389) (open) | `PDBMLParser` reads no resolution from any file (`if element:` on a leaf) | [#293](https://github.com/BioAGI-org/biop-ai/pull/293) | Adopted from #5389. Its `arctan2` line turns Z onto −Z for a zero-length vector whose z is −0.0; the fork adds `0.0` to z | — |
| [#2814](https://github.com/biopython/biopython/issues/2814) | `Bio.Nexus` cannot read `SYMBOLS = "0 1 2"`, and `AlignIO`, `SeqIO` and `Bio.Align` cannot read any `datatype=standard` matrix | [#294](https://github.com/BioAGI-org/biop-ai/pull/294) | Written here | An interleaved standard matrix raises `TypeError`; such an alignment cannot be written as NEXUS; Mesquite files with no `NTAX`, or with `BLOCKID`, fail |
| [#2449](https://github.com/biopython/biopython/issues/2449) | The GenBank writer turns a WGS, TSA or TLS master record's `rc` into `bp` and drops its accession-range lines | [#295](https://github.com/BioAGI-org/biop-ai/pull/295) | Written here | `Bio.GenBank.Record.__str__` has the same two gaps |
| [#4969](https://github.com/biopython/biopython/issues/4969); PR [#5206](https://github.com/biopython/biopython/pull/5206) (closed) | An Ensembl insertion written `7..6` on a linear record is read as spanning the origin, the location becomes `None`, and writing the record crashes | [#296](https://github.com/BioAGI-org/biop-ai/pull/296) | Written here, with #5206's condition apart from `1..0`; #5206 has no licence statement | — |
| PR [#3016](https://github.com/biopython/biopython/pull/3016) (open) | The Newick writer renames a named internal clade that has a confidence: `N1` with 0.9 becomes `N10.90` | [#297](https://github.com/BioAGI-org/biop-ai/pull/297) | Written here. #3016's guard is always true, so it would also rewrite every unnamed confidence as a branch length | An unnamed clade read with `comments_are_confidence=True` reads back as a clade named `99.00`; see Needs work |
| [#4878](https://github.com/biopython/biopython/issues/4878); PR [#5350](https://github.com/biopython/biopython/pull/5350) (open) | The GenBank writer puts the `ss-`/`ds-` strand prefix in the molecule-type columns of the LOCUS line | [#298](https://github.com/BioAGI-org/biop-ai/pull/298) | Tests adopted from #5350; writer fix written here. #5350's own writer change still pads a type to 7 columns after the prefix columns, so `snRNA`, `snoRNA` and `circRNA` give LOCUS lines of 80 to 82 characters that `SeqIO` cannot read back | Reading a long-header LOCUS line keeps only `DNA` from `ss-DNA`; `Bio.GenBank.Record` writes the old layout |
| [#2477](https://github.com/biopython/biopython/issues/2477) | The GenBank writer puts structured comments before the free-text comment | [#299](https://github.com/BioAGI-org/biop-ai/pull/299) | Written here | A RefSeq line of free text after the tables, such as `COMPLETENESS: full length.`, is now written before them; the EMBL writer drops structured comments |
| [#4999](https://github.com/biopython/biopython/issues/4999) | `cal_dn_ds` on codons that gaps put out of frame raises misleading errors or returns wrong numbers | [#300](https://github.com/BioAGI-org/biop-ai/pull/300) | Written here. The same PR fixes `build()` losing `rf_table`, which has no report, queued as row 44 | `get_full_rf_table` and `full_translate` still raise `IndexError` on such rows; column slices of a `build()` alignment lose `rf_table` |
| [#3751](https://github.com/biopython/biopython/issues/3751) | `BgzfWriter` accepts a read-only or text handle and fails only at the first full block, `flush()` or `close()` | [#301](https://github.com/BioAGI-org/biop-ai/pull/301) | Written here, with the `write(b"")` probe suggested on the issue | A `gzip.GzipFile` handle raises `TypeError` on Python 3.10 to 3.12; see Needs work |
| [#3532](https://github.com/biopython/biopython/issues/3532) | MEME XML from a run on more than 1000 sequences, which has no sites, cannot be read | [#302](https://github.com/BioAGI-org/biop-ai/pull/302) | Written here | Those motifs have no counts, so consensus, PWM and the like fail on them |
| [#4988](https://github.com/biopython/biopython/issues/4988); PRs [#4990](https://github.com/biopython/biopython/pull/4990) (open), [#5213](https://github.com/biopython/biopython/pull/5213) (closed) | `degenerate_consensus` gives `V` (or `D`, `B`) for a column with no A, C, G or T counted | [#303](https://github.com/BioAGI-org/biop-ai/pull/303) | Written here | A column of only `Y` gives `N`, not the `Y` the reporter asked for; a design choice |
| [#555](https://github.com/biopython/biopython/issues/555) | `dist_pearson` raises `UnboundLocalError`, or silently gives a wrong offset, when a matrix has `-inf` | [#304](https://github.com/BioAGI-org/biop-ai/pull/304) | Written here (`ValueError`) | `dist_pearson_at` raises `ZeroDivisionError` on an overlap with no variance |
| PR [#3022](https://github.com/biopython/biopython/pull/3022) (open) | The Newick writer raises `RecursionError` from a depth of about 496 | [#305](https://github.com/BioAGI-org/biop-ai/pull/305) | Adopted from #3022 | `ladderize` and NEXUS output still recurse, from a depth of about 1000 |
| [#5257](https://github.com/biopython/biopython/issues/5257); PR [#5279](https://github.com/biopython/biopython/pull/5279) (open) | A request: radius of gyration for `Structure`, `Model`, `Chain` and `Residue` | [#306](https://github.com/BioAGI-org/biop-ai/pull/306) | Adopted from #5279. As it stands #5279 sums in float32, so the answer depends on where the structure sits, and its docstring example uses `PDBParser` on an mmCIF file. The same PR fixes a `center_of_mass` crash, queued as row 48 | `center_of_mass` still averages in float32 |
| [#1022](https://github.com/biopython/biopython/issues/1022); PR [#5278](https://github.com/biopython/biopython/pull/5278) (open) | An unquoted hyphenated taxon name in a NEXUS `TRANSLATE` table raises `NexusError` | [#307](https://github.com/BioAGI-org/biop-ai/pull/307) | Adopted from #5278 | Hyphenated names in `MATRIX` rows are still split at the hyphen |
| [#5329](https://github.com/biopython/biopython/issues/5329); PRs [#5362](https://github.com/biopython/biopython/pull/5362) (closed), [#5330](https://github.com/biopython/biopython/pull/5330) (closed) | `Seq({}, length=N)` builds a partly defined sequence with no segments, which crashes on `+` and `*` | [#308](https://github.com/BioAGI-org/biop-ai/pull/308) | Adopted from #5362 | `Seq({0: ""}, length=4)` still builds zero-length segments; it does not crash |
| [#1016](https://github.com/biopython/biopython/issues/1016); PR [#5318](https://github.com/biopython/biopython/pull/5318) (open) | The EMBL reader ignores `DT` lines, so a record written as GenBank gets `01-JAN-1980` | [#309](https://github.com/BioAGI-org/biop-ai/pull/309) | Adopted from #5318. As it stands #5318 raises `IndexError` on a `DT` line with no date; the fork warns and skips it | The EMBL writer writes no `DT` lines |
| [#4291](https://github.com/biopython/biopython/issues/4291); PR [#5233](https://github.com/biopython/biopython/pull/5233) (open) | BioSQL with `driver="pymysql"` cannot load the schema, and MySQL 8 rejects `"rank"` | [#310](https://github.com/BioAGI-org/biop-ai/pull/310) | Adopted from #5233, minus its `db`/`passwd` hunk, which makes PyMySQL warn 440 times | — |
| [#729](https://github.com/biopython/biopython/issues/729); PR [#5255](https://github.com/biopython/biopython/pull/5255) (open) | `kegg_conv` rejects the documented entry conversions and `chebi`; several calls no longer match the KEGG REST API | [#311](https://github.com/BioAGI-org/biop-ai/pull/311) | Four of #5255's commits adopted; its `conv` and `list` commits were not, and those are written here | — |
| PR [#5194](https://github.com/biopython/biopython/pull/5194) (open) | A minimal MEME file with no background-frequency section raises `ValueError` | [#312](https://github.com/BioAGI-org/biop-ai/pull/312) | Adopted from #5194. As it stands #5194 raises `UnboundLocalError` on a file that ends after its `ALPHABET` line; the fork fixes that | — |
| [#2259](https://github.com/biopython/biopython/issues/2259) | `MMCIFParser` crashes with a bare `int()` `ValueError` on an `auth_seq_id` such as `20A` | [#313](https://github.com/BioAGI-org/biop-ai/pull/313) | Written here (`PDBConstructionException` naming the value) | `PDBParser` and `PDBMLParser` have the same bare error; `FastMMCIFParser` misses indented `_atom_site.` headers, see Needs work |
| [#3987](https://github.com/biopython/biopython/issues/3987); PRs [#5157](https://github.com/biopython/biopython/pull/5157) (open), [#5268](https://github.com/biopython/biopython/pull/5268) (open) | `PDBList` adds `/pub/pdb` to every server, so the EBI mirror gives 404 | [#314](https://github.com/BioAGI-org/biop-ai/pull/314) | Written here; neither PR agrees to the licence | — |
| [#1663](https://github.com/biopython/biopython/issues/1663); PR [#1664](https://github.com/biopython/biopython/pull/1664) (merged, writer only) | A GenBank or EMBL record whose description is `<unknown description>` fails to parse | [#315](https://github.com/BioAGI-org/biop-ai/pull/315) | Written here | The writers still write that description as empty |
| [#1624](https://github.com/biopython/biopython/issues/1624); PR [#1625](https://github.com/biopython/biopython/pull/1625) (closed) | The FASTA, FASTQ and QUAL writers raise `IndexError` on a whitespace-only description | [#316](https://github.com/BioAGI-org/biop-ai/pull/316) | Written here | — |
| [#2853](https://github.com/biopython/biopython/issues/2853) | `Bio.UniGene` fails on the blank lines between fields in NCBI's Xenopus files | [#317](https://github.com/BioAGI-org/biop-ai/pull/317) | Written here | `read()` still rejects a blank line after the final `//`; `HOMOL NO` sets `True`, see Needs work |
| [#4144](https://github.com/biopython/biopython/issues/4144); PRs [#4138](https://github.com/biopython/biopython/pull/4138) (merged), [#3931](https://github.com/biopython/biopython/pull/3931) (merged) | A v3 CEL file whose DatHeader scanner ID has no scanner type cannot be read; fixed upstream by #4138, without the test asked for on the issue | [#318](https://github.com/BioAGI-org/biop-ai/pull/318) | Written here: the test. The PR's main fix, a subfield-count check that is built and never raised, has no report and is queued as row 47 | — |
| [#4926](https://github.com/biopython/biopython/issues/4926); PR [#3031](https://github.com/biopython/biopython/pull/3031) (merged, other parsers) | Ten more parsers fail on a `tempfile.NamedTemporaryFile` handle, which has no `__next__` | [#319](https://github.com/BioAGI-org/biop-ai/pull/319) | Written here | `Primer3.parse` loops forever on empty input; see Needs work |
| [#969](https://github.com/biopython/biopython/issues/969); PR [#970](https://github.com/biopython/biopython/pull/970) (closed) | `Polypeptide.get_ca_list`, `get_tau_list` and `get_theta_list` raise `KeyError` for a residue with no CA | [#320](https://github.com/BioAGI-org/biop-ai/pull/320) | Written here; #970 has no licence answer | — |
| [#4652](https://github.com/biopython/biopython/issues/4652); PR [#5380](https://github.com/biopython/biopython/pull/5380) (open) | The Tutorial has no example of finding residues near a residue with `NeighborSearch` | [#321](https://github.com/BioAGI-org/biop-ai/pull/321) | Adopted from #5380 | — |
| [#2828](https://github.com/biopython/biopython/issues/2828); PR [#5386](https://github.com/biopython/biopython/pull/5386) (open) | The ruff settings live only in the pre-commit hook, so a bare `ruff check` disagrees with CI | [#322](https://github.com/BioAGI-org/biop-ai/pull/322) | Three of #5386's four commits adopted | — |
| PR [#5374](https://github.com/biopython/biopython/pull/5374) (closed) | The Clustal and MSF readers raise `IndexError` on a blank first line | [#323](https://github.com/BioAGI-org/biop-ai/pull/323) | Written here, in all four readers; #5374 covered `ClustalIO` only, without a licence | `Bio.Align`'s FASTA and A2M readers raise `IndexError` on blank lines; see Needs work |
| [#5292](https://github.com/biopython/biopython/issues/5292); PRs [#5293](https://github.com/biopython/biopython/pull/5293) (closed), [#5300](https://github.com/biopython/biopython/pull/5300) (closed) | `Seq` ordering against an unsupported type passes the internal `bytes` on, so the error names `bytes` and a reflected method receives `bytes` | [#324](https://github.com/BioAGI-org/biop-ai/pull/324) | Written here, the approach of #5293. The maintainer preferred only rewording the error; returning `NotImplemented` is a behaviour choice | — |
| [#2839](https://github.com/biopython/biopython/issues/2839) | The PHD reader crashes on a TRIM line with a split exponent, `-1 -1 -1.000000e 000` | [#325](https://github.com/BioAGI-org/biop-ai/pull/325) | Written here, with the `BiopythonParserWarning` a maintainer said was acceptable | — |

## Checked and not reported

- **`nt_search(seq, "")` never terminates** (queued row 11). Filed upstream
  by a third party as
  [#5406](https://github.com/biopython/biopython/issues/5406) on 2026-11-07,
  so there is nothing to send. Fixed here; see `NEWS.rst`.
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
- **The second commit of [#5127](https://github.com/biopython/biopython/pull/5127)
  breaks `coord_space`.** Already raised. peterjc posted mypy's output on
  2025-12-27, the day after the pull request opened, naming
  `sc[2][0]` at lines 572 and 619 — "Value of type `float` is not indexable" —
  and asking the author to avoid the double indexing. The thread did not move
  after that; the pull request was approved by JoaoRodrigues, and mdehoon
  closed it unmerged, without comment, on 2026-10-23.

  This fork adopted the first commit only, and confirmed the second fails at
  runtime with `IndexError: invalid index to scalar variable` against upstream
  `master` `e136be720`, whose `test_PDB_vectors` still passes flat arrays.

  Two things neither the linter nor the thread has, kept here in case it
  revives. First, the root cause: each changed line sits under a commented-out
  reference line using `[0][0]`, written when the callers passed column arrays;
  the callers switched to flat ones and the comments were never updated, so the
  commit is faithfully following stale documentation. Second, **the suggested
  fix does not work either** — `get_spherical_coordinates` returns a *tuple*,
  so `sc[2, 0]` raises `TypeError: tuple indices must be integers or slices,
  not tuple`. The correct expression is the original `sc[2]`.

- **The BioSQL tests shared from `common_BioSQL` run with the wrong driver.**
  In one pytest session they read module globals that each driver module's
  `check_config` overwrites when it is imported, so all of them run with the
  last driver imported, `sqlite3`, and the MySQL shared tests never reach
  MySQL (found in [#310](https://github.com/BioAGI-org/biop-ai/pull/310)). Not fixed here. Not an upstream defect:
  upstream's `run_tests.py` imports and runs one test module at a time
  (read, not run), so the globals are right for each; the fault comes from
  this fork running the suite under pytest.

## Limits of the survey

The duplicate check covered the GitHub issue and pull request tracker only. It
did not cover the `biopython-dev` mailing list or the pre-GitHub Redmine and
Bugzilla archives, so "unreported" here means "not in the GitHub tracker".
