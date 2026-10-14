# BioPAIthon Improvements Plan

A prioritised, evidence-backed plan for improving this codebase, produced by a
comprehensive review of the tree at the point BioPAIthon forked from Biopython
(`5d6fe8d22`, version `1.88.dev0`).

This is a living document. Each item cites the file and line it was found at so
it can be checked rather than believed. Where a claim was reproduced by running
code, it is marked **[reproduced]**; where it rests on reading the source, it is
marked **[from source]**. Nothing had been fixed when this was first written;
items fixed since carry an inline **Status** note citing the fork pull request,
with the original text kept below for the record.

**Effort** is S (hours), M (days), L (weeks). **Impact** is this fork's judgement
of user-visible value.

---

## Tier 0 — Correctness bugs

These produce wrong answers or destroy diagnostic information. They are cheap to
fix and should go first. Several are silent, which is what makes them serious:
users get a plausible number rather than an error.

### 0.1 `ProteinAnalysis.flexibility()` ignores the residue at the centre of its window **[reproduced — FIXED]**

> **Status: fixed** in PR #19. `Bio/SeqUtils/ProtParam.py` now uses
> `middle = subsequence[window_size // 2]` and
> `range(self.length - window_size + 1)`, the docstring documents the
> behaviour change, and the fixture was regenerated.

`Bio/SeqUtils/ProtParam.py:168-180`. The sliding window is nine residues wide,
but the "middle" index is computed as `window_size // 2 + 1`, which is offset 5,
not 4. Offset 4 — the actual centre — is never read, and offset 5 is added twice.

Reproduced on a poly-alanine background by substituting a single tryptophan at
each window offset:

| substitution at window offset | change in score |
|---|---|
| 3 | −0.012381 |
| **4 (the centre)** | **0.000000** |
| 5 | −0.027619 |

The centre residue has no effect on its own flexibility score, and its
neighbour has 2.23× the weight it should — matching the ratio 1.8125 / 0.8125
predicted by reading the code. Separately, the loop is `range(self.length -
window_size)`, so the final window is dropped: a 9-mer, which contains exactly
one complete window, returns `[]`.

Every Vihinen-1994 flexibility profile this library has produced is affected.
`Tests/test_ProtParam.py:97-146` hard-codes the current output as a 143-element
literal, so the test suite pins the bug rather than catching it.

**Fix:** `middle = subsequence[window_size // 2]` and `range(self.length -
window_size + 1)`. Regenerate the fixture from an independent reference (ExPASy
ProtScale), not from the code. Add a symmetry test asserting that substituting
at offset *k* and at *window_size−1−k* moves the score identically.
**Effort S · Impact high**

### 0.2 `Bio.Seq.translate()` silently discards `gap=` for `Seq` inputs **[reproduced — FIXED]**

> **Status: fixed** in PR #1. All three branches of the module-level
> `translate()` (`Bio/Seq.py:3009-3015`) now forward `gap=gap`. The two
> defaults still differ (`None` at module level, `"-"` on the method), which
> that PR deliberately left alone — see the note on upstream #4994 under
> "Do not adopt".

`Bio/Seq.py:3006-3013` forwards to `sequence.translate(table, stop_symbol,
to_stop, cds)` without passing `gap`. The module-level function defaults to
`gap=None` (`Bio/Seq.py:2913`) while the method defaults to `gap="-"`
(`Bio/Seq.py:1521`), so the two disagree even before the argument is dropped.

Reproduced:

```
translate("ATG+++TAA", gap="+")       -> "M+*"
translate(Seq("ATG+++TAA"), gap="+")  -> TranslationError: Codon '+++' is invalid
```

Same input, same argument, different answer depending on whether the caller
wrapped the string in a `Seq`. A user translating a gapped codon alignment is
told the gap character they just declared is an invalid codon.

**Fix:** pass `gap=gap` through both branches and reconcile the two defaults.
Add a regression test asserting `translate(s, gap=g) == translate(Seq(s),
gap=g)` across a matrix of gap characters.
**Effort S · Impact high**

### 0.3 Derived `SeqRecord`s share mutable annotations with their parent **[reproduced — annotations FIXED]**

> **Status: the annotations half is fixed.** A `_copy_annotations` helper now
> deep-copies the annotations dictionary at all six derive sites, and the
> `print`-then-reraise in `__add__`'s per-letter concatenation handler is now a
> `TypeError` naming the key and the two types. `DerivedRecordIsolation` in
> `Tests/test_SeqRecord.py` covers all six methods and fails against the old
> shallow copy. **The features half is deliberately not done here:** `_flip`
> and `_shift` already return new features, so `reverse_complement` and the
> shifted half of `__add__` do not share, but `features[:]` in `upper`, `lower`
> and the other `__add__` paths does. Deep-copying features is far more
> expensive than annotations (measured ~675 µs for 41 features versus ~34 µs
> for a whole annotations dict, and it scales with feature count on a genome
> record), so it belongs in its own change with its own benchmark rather than
> bundled into a metadata-correctness fix.

`Bio/SeqRecord.py:994-1002` (`__add__`), and identically at `:1083`
(`__radd__`), `:1136` (`upper`), `:1187` (`lower`), `:1418`
(`reverse_complement`), `:1527` (`translate`). Each uses
`annotations=self.annotations.copy()` — a shallow copy. Annotation *values* in
parsed records are routinely lists (`keywords`, `accessions`, `references`,
`taxonomy`). `features[:]` likewise copies the list but shares the `SeqFeature`
objects.

Reproduced:

```python
r = SeqRecord(Seq("ACGT"), annotations={"keywords": ["a", "b"]})
r2 = r + "TT"
r2.annotations["keywords"].append("EVIL")
r.annotations["keywords"]  # -> ['a', 'b', 'EVIL']
```

The whole point of these methods returning a new record is that the original is
untouched. A pipeline that derives edited copies, adjusts their keywords or
feature qualifiers, then writes both out silently corrupts the originals — and
the mutation site is far from the corrupted object, so it is very hard to debug.

While here: `Bio/SeqRecord.py:1042` does `print("Failed while try to
concatenate letter annotations")` from inside a library exception handler.
Libraries should not write to stdout, and the message has a typo.

**Fix:** deep-copy annotation values (or a documented `_copy_annotations()`
helper that deep-copies containers and shares immutables); copy features
explicitly. Delete the `print` and re-raise with the offending key and types.
**Effort M · Impact high**

### 0.4 `Seq.search()` never matches at the final position **[reproduced — bug FIXED]**

> **Status: the off-by-one is fixed** in PR #7 — the loop at `Bio/Seq.py:997`
> is now `range(len(self) + 1)`. The separate rewrite suggested below
> (`bytes.find` cursors, deterministic ordering) has not been done: the loop
> still allocates a slice per (position, length) pair and yield order across
> differing lengths still depends on dict insertion order.

`Bio/Seq.py:997` iterates `range(len(self) - 1)`, one short.

```
Seq("ACGTA").search(["A"])  -> [(0, 'A')]      # the 'A' at index 4 is dropped
```

Substrings of length ≥ 2 are unaffected, which is why this has survived.
`search()` is the API that `find`/`rfind` docstrings point users to for scanning
multiple motifs, so this is a silently wrong answer on the recommended path.

The same six lines are also O(n·k): they allocate a fresh `bytes` slice per
(position, length) pair instead of using C-level `bytes.find`, and yield order
across differing lengths depends on dict insertion order.

**Fix:** `range(len(self))`. Separately, reimplement over `bytes.find` cursors
merged by position so results are ordered deterministically, and document the
ordering. Test a match at index `len-1` for each substring length.
**Effort S (bug) / M (rewrite) · Impact high**

### 0.5 `as_handle()` swallows `TypeError` raised by its caller's block **[from source, with repro — FIXED]**

> **Status: fixed** in PR #6. The `try` in `Bio/File.py` now wraps only the
> `open()` call; the `yield` sits outside it, so a `TypeError` from the
> caller's block propagates unchanged.

`Bio/File.py:71-75`. The `try` is meant to catch `TypeError` from `open()`
(meaning "this is already a handle"), but it also wraps the `yield`, so *any*
`TypeError` raised inside the caller's `with` body is caught and the generator
yields a second time:

```
with File.as_handle("/etc/hostname") as fp:
    raise TypeError("a real parser error")
# -> RuntimeError: generator didn't stop after throw()
```

`as_handle` is used by 16 modules including `Bio/SeqIO/QualityIO.py`,
`Bio/AlignIO/__init__.py`, `Bio/SearchIO/__init__.py`, `Bio/GenBank/Scanner.py`,
`Bio/Nexus/Nexus.py`, `Bio/Phylo/_io.py` and six `Bio/PDB/*` modules. Any
`TypeError` a parser raises on malformed input becomes a message with zero
diagnostic value, with the original traceback destroyed.

**Fix:** narrow the guard to the `open()` call only — `try: fp = open(...)
except TypeError: yield handleish else: with fp: yield fp`. Add a regression
test that a `TypeError` from inside the block propagates unchanged.
**Effort S · Impact high**

### 0.6 `SeqRecord.reverse_complement()` crashes on features with `UnknownPosition` **[from source — FIXED]**

> **Status: fixed** in PR #5, with exactly the sentinel proposed below: the
> sort key in `Bio/SeqRecord.py` (now ~:1446) returns `(0, int(start))` for
> known positions and `(1, 0)` for unknown ones.

`Bio/SeqRecord.py:1406-1413`. The sort key catches `TypeError` "expected for
UnknownPosition" and returns `None` — but a key list containing `None` cannot be
ordered, so `list.sort` raises `TypeError: '<' not supported between instances
of 'NoneType' and 'NoneType'`. The `except` clause defeats its own purpose.

`UnknownPosition` is produced by the GenBank/EMBL parsers for `?` locations, so
this fires on real files with no obvious connection to the user's action.

**Fix:** return a total-order sentinel — `(0, int(start))` for known and `(1, 0)`
for unknown — so unknown positions sort stably to the end.
**Effort S · Impact medium**

### 0.7 `SeqIO.index()` and `SeqIO.parse()` disagree on the same file **[from source — FIXED]**

> **Status: this case is fixed** in PR #8. `Bio/SeqIO/_index.py` (~:220) now
> yields an empty ID for a bare title line in the fasta and pir formats,
> matching the parsers. The structural duplication that caused it — the index
> layer re-deriving IDs independently of each parser — has since been removed
> as well: PR #88 moved the marker and id rules onto the parser classes
> (§1.4), so this whole category of disagreement is closed, not just this
> instance.

`Bio/SeqIO/_index.py:219` re-derives the record ID independently of the parser
(`line[marker_offset:].strip().split(None, 1)[0]`), while
`Bio/SeqIO/FastaIO.py:245-249` explicitly handles a bare `>` line. On a FASTA
record with an empty description, `parse` yields `id=''` and `index` raises a
bare `IndexError: list index out of range` with no filename or line number.

`Bio/File.py:243-244` (`raise ValueError(f"Key did not match ({key} vs
{key2})")`) exists solely to police this duplication after the fact. See §1.4
for the structural fix.
**Effort S (this case) · Impact medium**

### 0.8 343 `assert` statements validate untrusted file content **[sweep FIXED, lint hook beyond the parsers open]**

> **Status: the sweep is done, and a lint hook guards the parser packages.**
> PR #69 and the adopted upstream #5175 converted the first twelve sites.
> The sweep then converted the asserts that check input to `ValueError`
> (or `TypeError` for a wrong type), and listed each kept assert with the
> reason it holds: #155 `Bio.Align`, #156 `Bio.SearchIO`, #157 `Bio.Blast`,
> #159 `Bio.SeqIO`, #160 `Bio.AlignIO`, and #165 the rest of `Bio`,
> including `Bio/GenBank/Scanner.py:479`. A second `ruff-check` hook,
> `ruff S101 (parser packages)`, now fails any new `assert` in
> `Bio/{Align,AlignIO,Blast,SearchIO,SeqIO}/`. It skips
> `Bio/Align/substitution_matrices/`, `Bio/Align/analysis.py` and the
> deprecated `Bio/Blast/NCBI{XML,WWW}.py`. The 132 kept asserts there carry
> `# noqa: S101`.
>
> Open:
>
> - **GenBank is not guarded yet.** Extending the hook's `files` to
>   `Bio/GenBank/` needs `# noqa: S101` on its 15 kept asserts (#165 lists
>   them); the rest of `Bio/` could follow.
> - **One kept assert is not an invariant.** #156 kept
>   `assert gap_len > 0` in `Bio/SearchIO/InfernalIO/infernal_text.py` on
>   purpose: Infernal may print `*[ 0]*` on the hit line, which is valid.
>   Relaxing it to allow 0 needs a real Infernal file to test against.
> - **Offset asserts in `Bio/SeqIO/_index.py` fail on BGZF.** The
>   IntelliGenetics and UniProt XML indexers check offsets built by adding
>   up line lengths against `handle.tell()`. On a BGZF file of more than one
>   block, `tell()` is a virtual offset, so `SeqIO.index` raises
>   `AssertionError` (`:555`, `:462`). Under `python -O` the UniProt index
>   works and the IntelliGenetics one raises a misleading BGZF error.

`grep -c "^\s*assert " Bio/{SeqIO,AlignIO,Align,GenBank}/*.py` → 343. These are
not internal invariants; they check parsed input — `Bio/SeqIO/_index.py:671`
(`assert line[0:1] == b"+"`, the FASTQ separator), `:443` (offset arithmetic),
`:395`, `:231`, `Bio/GenBank/Scanner.py:479`, `Bio/SeqIO/FastaIO.py:247`.

Under `python -O` — common in containerised pipelines — every one vanishes. A
malformed FASTQ then produces a *silently wrong index*: `get_raw()` returns the
wrong bytes rather than raising. Even with asserts enabled,
`AssertionError: b'garbage\n'` is a poor error with no file or offset.

**Fix:** convert content-checking asserts to `raise ValueError(...)` with
filename and offset. Keep asserts only for provably internal invariants. Add
ruff `S101` scoped to parser modules to prevent regression.
**Effort M · Impact medium**

### 0.9 Malformed GenBank is a warning where malformed EMBL is an error **[length check FIXED, `strict=` open]**

> **Status: the length-mismatch half is fixed.** PR #71 made a
> declared-vs-actual sequence length mismatch a `ValueError` naming the
> record and both lengths, for `Bio.SeqIO` (genbank, embl, imgt) and
> `Bio.GenBank.FeatureParser` — a truncated download no longer yields a
> silently short `SeqRecord`. Two malformed-but-recoverable shapes still
> parse with a warning as before: complete sequence data missing only the
> final `//`, and a record with no sequence data at all. The broader
> `strict=` option promoting parser warnings to errors remains unbuilt.

Two branches of the same file. EMBL (`Bio/GenBank/Scanner.py:656-662`) raises
`ValueError("Premature end of file in sequence data")`. GenBank
(`:1213-1223`) warns and substitutes `line = "//"`, i.e. silently truncates. The
declared-vs-actual length check is also only a warning
(`Bio/GenBank/__init__.py:834-843`).

A truncated GenBank download yields a `SeqRecord` whose `.seq` is shorter than
its LOCUS length, with a warning that is trivially lost in a pipeline. There are
88 `BiopythonParserWarning` sites across `Bio/`, 23 in `Scanner.py` alone, and
no `strict=` option anywhere.

**Fix:** reconcile the two paths; make length mismatch an error. Add
`SeqIO.parse(..., strict=True)` promoting parser warnings to `ValueError`, and
include record ID and stream position in every message.
**Effort M · Impact high**

### 0.10 Two more error handlers report the wrong cause **[from source — FIXED]**

> **Status: both halves are fixed.** PR #13 made `Bio/SeqIO/__init__.py`
> (~:974) validate the filename type up front and construct the proxy outside
> any `try`, so construction errors propagate. PR #11 hoisted the writer
> *construction* in `Bio/Align/__init__.py` (~:2409) out of the guard, so the
> `except AttributeError` now covers only the bare `module.AlignmentWriter`
> lookup.

`Bio/SeqIO/__init__.py:936-941` wraps proxy *construction* and reports any
`TypeError` as "Need a string or path-like object for the filename (not a
handle)", with `from None` erasing the cause — so a corrupt SFF header is
reported as a filename-type mistake. `Bio/Align/__init__.py:2384-2389` wraps
writer construction and reports any `AttributeError` as "Formatting alignments
has not yet been implemented for the {fmt} format" — so a bad keyword argument
is reported as an unsupported format. The same check is done correctly on the
bare attribute lookup at `Bio/Align/__init__.py:4871-4876`.

**Fix:** validate the filename type up front and let construction errors
propagate; hoist the writer lookup out of the `try`.
**Effort S · Impact medium**

### 0.11 Heap buffer overflow in `_bcif_helper` driven by untrusted file metadata **[reproduced by review — FIXED]**

> **Status: fixed** in PR #3. Every unpack function in
> `Bio/PDB/bcifhelpermodule.c` now reads `out_view->shape[0]` and bounds
> `out_index` against it, failing cleanly when the two file fields disagree
> in either direction.

**This is the most serious item in this document and should be fixed first.**

`Bio/PDB/bcifhelpermodule.c:6-33` (and the identical `_u16`/`_i8`/`_i16`
variants at `:36,66,96`) bound their loop by the *input* buffer only. The output
index is never checked against the destination:

```c
Py_ssize_t in_size = in_view->shape[0];
uint32_t *out_data = out_view->buf;
while (in_index < in_size) {
    ...
    out_data[out_index] = sum;   /* out_view->shape/len never consulted */
    in_index += 1; out_index += 1;
}
```

`integer_unpack` (`:126-181`) validates only `ndim` and the input `format`. The
caller sizes the destination from the file itself —
`Bio/PDB/binary_cif.py:110-123` does `src_size = encoding["srcSize"]` then
`np.empty((src_size,), dtype)`. But `srcSize` and the packed-data length are two
independent fields of a `.bcif` file, so a corrupt or hostile structure file
makes them disagree and the extension writes past the array. Confirmed: a
200,000-byte input with a 1-element output **segfaults the interpreter**
(exit 139, core dumped).

Related in the same file: passing a 2-D input sets `ValueError` then falls
through the `exit:` label at `:176-180`, which returns `Py_None` regardless →
`SystemError: returned a result with an exception set`.

**Fix:** pass `out_view->shape[0]` into each unpack function and bound
`out_index`; verify the output's itemsize and format; raise `ValueError` on
truncation. Change `exit:` to return `NULL` when `PyErr_Occurred()`. Add a
regression test with a deliberately wrong `srcSize`.
**Effort S · Impact high (memory safety)**

### 0.12 Out-of-bounds read and a dead validation check in `_aligncore` **[reproduced by review — FIXED]**

> **Status: all three defects are fixed** in PR #2. `offset` is validated
> against `PyBytes_GET_SIZE(line)` and raises `ValueError` when outside it;
> the tautological length check is now `else if (m != self->m)` and raises;
> and the error paths return `NULL` with `PyErr_NoMemory()` set on allocation
> failure.

`Bio/Align/_aligncore.c:141-143` uses a caller-supplied `offset` unvalidated:
`buffer = PyBytes_AS_STRING(line) + offset`. Confirmed:
`PrintedAlignmentParser().feed(b'ACGT', 100000)` returns unrelated heap memory
as the "sequence", and `feed(b'ACGT', 1<<40)` **segfaults**.

Second, at `:187-194`, the length check is a tautology:

```c
m = s - buffer;
if (n == 0) self->m = m;
else if (buffer + m != s) {   /* always false: m was just defined as s - buffer */
```

The intended comparison is `m != self->m`. Confirmed: feeding `b'ACGT-'` then
`b'AC'` is accepted silently and yields `shape == (2, 4)`. `_aligncore` is on the
default path for `Bio.Align` printed-alignment parsing
(`Bio/Align/__init__.py:1057`), so ragged alignment blocks produce silently
wrong coordinate arrays instead of an error.

Third, `Parser_fill` (`:273,279-286,322-327`) reaches `Py_RETURN_NONE` on both
allocation failure and its `ValueError` path.

**Fix:** bounds-check `offset` against `PyBytes_GET_SIZE(line)`; fix the length
comparison; return `NULL` from every error path and set `PyErr_NoMemory()`.
**Effort S · Impact high (memory safety)**

### 0.13 `ccealign` leaks ~116 KB per call and segfaults on tuple coordinates **[reproduced by review — FIXED]**

> **Status: now fully fixed.** PR #68 closed everything the note below left
> open: the `CEAlignment` struct-sequence type is created once at module load
> and exposed as `Bio.PDB.ccealign.CEAlignment`, so all results share one
> type and pickle; every allocation is checked and raises `MemoryError`
> instead of crashing; coordinate lengths the extension's arithmetic cannot
> hold are rejected up front; and `PyInit_ccealign` carries `PyMODINIT_FUNC`.
> The history of the leak and validation fixes is below.

> **Status of the earlier fixes: partly fixed, in two PRs.** PR #16 removed the three `Py_INCREF`s
> and frees every surviving `pathBuffer[]` entry before returning (~:690), so
> both the object and byte leaks are gone — the measurements are in the
> update blocks below. PR #39 added argument validation: coordinates go
> through `PySequence_Size`/`PySequence_GetItem` with every lookup checked
> (so tuple coordinates and short inner lists raise instead of segfaulting —
> the two later-found crashes below), `fragmentSize` and `gapMax` are
> range-checked (so `run_cealign(c, c, 0, 30)` raises), and the
> `PyArg_ParseTuple` return is checked. **Still open:** the fresh
> `PyStructSequence_NewType` built per path inside the result loop (~:677,
> so `type(r[0]) is type(r[1])` is still `False`), the unchecked
> `PyMem_RawMalloc`s in `calcDM`, `calcS` and `findPath` (:108, :111, :194,
> :197, :472), and the missing `PyMODINIT_FUNC` on `PyInit_ccealign`
> (~:727).

`Bio/PDB/ccealignmodule.c:610-611,617-618,636-637` contain three `Py_INCREF`s on
references that are immediately stolen (by `Py_BuildValue("[NN]", ...)` and
`PyStructSequence_SetItem`). Measured on 200-residue chains over 100 calls:
**116.5 KB and 261 tracked objects leaked per `run_cealign()` call**. An
all-against-all `CEAligner` run over 10,000 pairs leaks roughly 1.1 GB.

> **Correction.** An earlier version of this entry claimed that removing those
> three `Py_INCREF`s reduced the leak to "0 objects, ~0 KB per call, nothing
> else changed". **That was wrong**, and a later adversarial review of PR #16
> disproved it by measuring three builds side by side:
>
> | build | objects/call | RSS KB/call | tracemalloc KB/call |
> |---|---|---|---|
> | unpatched | 261.1 | 230.1 | 107.9 |
> | `Py_INCREF`s removed | 0.1 | 70.4 | 31.7 |
> | plus `pathBuffer` freed | 0.1 | 11.7 | 0.5 |
>
> The object leak is fixed; **most of the byte leak is not**. The residual
> 31.7 KB/call is exactly `MAX_PATHS(20) × sizeof(afp)(8) × n(200)` = 31.25 KiB:
> `curPath` is `PyMem_RawMalloc`'d at `:451` and stored into `pathBuffer[]` at
> `:574`, but freed at `:578` **only** in the `bufferSize == MAX_PATHS` overflow
> case. The surviving entries are never freed before `return result;` at `:672`.
> The complete fix needs, before that return:
> `for (int i = 0; i < bufferSize; i++) PyMem_RawFree(pathBuffer[i]);`
>
> The original claim came from re-measuring only tracked-object counts and
> reporting bytes as though they had been measured too. Treat single-number
> "leak eliminated" claims in this document with suspicion unless a
> before/after byte measurement is shown.
>
> **Update: the byte leak is now fixed** on PR #16, which adds exactly the loop
> above. Re-measured across three builds in one session — `main`, the branch
> with the new free deleted, and the branch as submitted:
>
> | build | objects/call | RSS KB/call | tracemalloc KB/call |
> |---|---|---|---|
> | `main` | 261.00 | 211.0 | 149.28 |
> | branch minus the new free | 0.00 | 60.5 | 31.70 |
> | branch as submitted | 0.00 | 25.8 | **0.45** |
>
> Confirmed independently by LeakSanitizer, which attributes 576,000 bytes to
> `findPath` `:451` on `main` and **0** on the branch, and by an asymptotic
> scaling run showing 0 B/call retained at 1,280 calls. The residual RSS is
> transient and plateaus; it is not retention.
>
> One methodological correction to the note above: `tracemalloc` **does** track
> `PyMem_RawMalloc` on CPython 3.12 — it hooks `PYMEM_DOMAIN_RAW`. An earlier
> claim here that it does not was wrong.
>
> The other five defects listed below remain unfixed, and two further crashes
> were found later that are not listed at all: `run_cealign(c, c, 0, 30)` and
> a coordinate list whose inner lists are too short both segfault.

Also: `:660-661` builds a fresh heap type inside the loop (20 distinct
`CEAlignment` types per call, so `type(r[0]) is type(r[1])` is `False` and the
results will not pickle); `:381-386` calls `PyList_GetItem` with no NULL check,
so passing tuples instead of lists **segfaults**; `:685` discards the
`PyArg_ParseTuple` return; every `PyMem_RawMalloc` at `:107,110,193,196,451` is
unchecked; `:751` is missing `PyMODINIT_FUNC`.

**Fix:** delete the three `Py_INCREF`s; hoist the struct-sequence type to a
module-level singleton created in `PyInit_ccealign`; NULL-check everything.
**Effort S (leak) / M (hardening) · Impact high**

### 0.14 `cpairwise2.rint` writes 8 bytes into a 4-byte stack slot **[reproduced by review — FIXED]**

> **Status: fixed** in PR #9. `Bio/cpairwise2module.c:429` now parses with
> `"d|i"`, and no `"l"` format remains anywhere in the file.

`Bio/cpairwise2module.c:415-425` declares `int precision` but parses it with
format `"l"`, which requires `long *`. On LP64 that is a 4-byte out-of-bounds
stack write. Confirmed indirectly: `cpairwise2.rint(2.5, 2**62)` is accepted and
returns `0`, where an `"i"` conversion would have raised `OverflowError` — so the
parser is demonstrably storing 64 bits.

**Fix:** change the format to `"d|i"`. Then sweep every `PyArg_Parse*` call site
in the tree for width mismatches between format character and declared C type.
Note this is reachable only from the deprecated `Bio.pairwise2` (§2.3), so
deleting that module resolves it too.
**Effort S · Impact medium**

---

## Tier 1 — Structural work

Large, high-leverage changes. Each removes a whole category of future bug rather
than one instance.

### 1.1 `Bio.AlignIO` and `Bio.Align` are two complete parser stacks for the same formats **[re-scoped — shim declined, new scope FIXED]**

> **Status: done, on the new scope.** All five `Bio.Align` reader PRs
> listed below have merged: PR #182 (non-seekable streams), PR #183
> (`#=GR` by column), PR #197 (interleaved Stockholm), PR #184 (several
> PHYLIP alignments in one file) and PR #195 (`phylip-relaxed`, the first
> `Bio.Align` format name that is not a module name). **Decided:** the
> `AlignIO` shim was dropped in favour of closing `Bio.Align`'s gaps
> (PR #142); the reasons are below. The notes below are kept as written
> on 2026-10-06, apart from the outcome column of the table.
>
> **Status: re-scoped on 2026-10-06. The shim is declined. The new goal is
> that `Bio.Align` reads the files `AlignIO` reads, closing the gaps
> measured on real inputs.** `Bio.AlignIO` is not rewritten; it takes only
> small fixes, such as PR #116's. This fork will not rebuild it as a shim
> over `Bio.Align`, delete its format modules, port `MafIndex`, add
> `Bio.Align.convert()`, or migrate the internal consumers. This note
> replaces the **Plan** in the original text, which is kept below for the
> record. All numbers below were measured on 2026-10-06 on `main` and on
> upstream `master` `372c71069`, which is also the merge-base. They were
> re-checked on 2026-10-08, on `main` after PR #155 and on upstream `master`
> `bb295f15a`; the two upstream commit counts are from that day, and line
> numbers are `main`'s. The commands that reproduce them are in PR #142.
>
> **Why the shim is declined:**
>
> - **It fails its own acceptance gate.** A prototype backed `AlignIO`'s ten
>   reader names with `Bio.Align`, with adaptors for every mechanical
>   difference (name and description defaults, MAF and Mauve start and strand
>   annotations, EMBOSS annotation keys, empty input). It fails 34 of 50 tests
>   in `test_AlignIO` and 17 of 31 in `test_AlignIO_MsfIO`. That comes from a
>   throwaway harness and is a lower bound: the Clustal, EMBOSS and Mauve test
>   modules call the iterator classes directly, so they were not exercised.
> - **Two differences destroy information, so no adaptor can recover them.**
>   An `Alignment` cannot hold an all-gap column:
>   `Alignment(*Alignment.parse_printed_alignment([b"A-C", b"A-C"])).shape` is
>   `(2, 2)`. So `Tests/Nexus/test_Nexus_input.nex` is 9×48 in `AlignIO`
>   (pinned by `test_reading_alignments_nexus1`) and 9×46 in `Bio.Align`. And in
>   MAF, `Bio.Align` reads `.` as a gap where `AlignIO` puts in the reference
>   letter (`Bio/AlignIO/MafIO.py:181-194`, deliberate since `9470216c4` in
>   2011, pinned by `test_reading_alignments_maf1`).
> - **The writers have to stay anyway.** Only nexus output is byte-identical
>   between the stacks. The eight `Bio/AlignIO` format modules hold 3,326
>   lines: 634 are writers and 614 are `MafIndex`. At most the 1,556 lines of
>   readers could go, and they encode exactly the behaviour the tests pin.
> - **`Bio.SeqIO` goes through `AlignIO`.** It reads 11 alignment format
>   names and writes 8 with `AlignIO`'s classes
>   (`Bio/SeqIO/__init__.py:473-558`), so a shim would change SeqIO output
>   too.
> - **The sync cost runs the wrong way.** `Bio/AlignIO` has had 1 upstream
>   commit since the merge-base (`git rev-list --count
>   372c71069..upstream/master -- Bio/AlignIO`), and it was a sweep: ruff
>   import sorting, which edits `NexusIO.py`
>   ([#5360](https://github.com/biopython/biopython/pull/5360), merged
>   2026-10-08). More sweeps that touch it are open
>   ([#5219](https://github.com/biopython/biopython/pull/5219) and
>   [#5361](https://github.com/biopython/biopython/pull/5361)), and
>   upstream still fixes its bugs
>   ([#5355](https://github.com/biopython/biopython/issues/5355)). Rewritten
>   modules would turn each of those into a conflict. Closing the shim's
>   gaps would also mean editing `Bio/Align/__init__.py`, upstream's
>   most-edited file: 67 upstream commits since 2023 (`git rev-list --count
>   --since=2023-01-01 upstream/master -- Bio/Align/__init__.py`).
> - **Deleting the modules would break this fork's own rule** — "Nothing has
>   been removed here that upstream still ships" (see "Where this fork
>   differs from upstream" below). The submodules are also in public use.
>   GitHub code search, with the method of "How much is any of this used?":
>   `from Bio.AlignIO.Interfaces import` 143 files, `MafIO.MafIndex` 116,
>   `from Bio.AlignIO.MafIO import` 28, `from Bio.AlignIO.PhylipIO import`
>   23, `from Bio.AlignIO.ClustalIO import` 14.
> - **Upstream keeps both stacks.** Its Tutorial calls
>   `MultipleSeqAlignment` and `AlignIO` "older"
>   (`Doc/Tutorial/chapter_msa.rst:6`), and mdehoon called
>   `MultipleSeqAlignment` "essentially obsolete"
>   ([#5051](https://github.com/biopython/biopython/pull/5051)). But nothing
>   is deprecated, and no issue proposes rebuilding `AlignIO` on
>   `Bio.Align`. The maintainers also disagree on how strict the parsers
>   should be ([#3747](https://github.com/biopython/biopython/issues/3747));
>   a shim would have to settle that alone.
>
> **The rest of the old plan:**
>
> - **`MafIndex` stays in `AlignIO`**, where PR #116 fixed its `close()`.
>   `AlignIO` keeps its own copy, so a port would add a third copy rather
>   than remove one. `Bio.Align`'s route to indexed MAF is bigMaf with
>   `search()`: `Tests/MAF/ucsc_mm9_chr10.maf` converts to a 50 KB bigMaf,
>   and searching mm9.chr10 3,014,000–3,015,000 returns 7 alignments. But
>   it is not a drop-in replacement. The bigMaf writer needs every
>   alignment in memory and `.targets` set by hand (it ignores its
>   `targets=` argument, `Bio/Align/bigmaf.py:135`), so it does not suit
>   the genome-scale MAF files that `MafIndex` is for. `MafIndex` reads
>   plain MAF only: `MafIndex.__init__` opens the file with the built-in
>   `open(..., "rb")`, so "bgzip-aware" in the original text below is
>   wrong.
> - **No `Bio.Align.convert()`.** `Align.write(Align.parse(a, f), b, g)` is
>   the whole function, and it would have to live in upstream's most-edited
>   file.
> - **No consumer migration.** "21 internal consumers" was stale.
>   `git grep -n AlignIO -- 'Bio/*.py' ':!Bio/AlignIO'` finds 2 functional
>   consumers: SeqIO's alignment-format bridge
>   (`Bio/SeqIO/__init__.py:473-558`) and `MultipleSeqAlignment.__format__`
>   (`Bio/Align/__init__.py:383-386`). Both are `AlignIO` by definition. The
>   rest is 6 doc examples (`Bio/Align/__init__.py:94,528`,
>   `Bio/Phylo/TreeConstruction.py:384,605,1177`, `Bio/Phylo/_utils.py:306`)
>   and prose.
>
> **New scope: `Bio.Align` reads what `AlignIO` reads.** Users who take
> upstream's advice and move to `Bio.Align` should not lose files. Only gaps
> measured on real inputs get closed, and §1.1 is done when the five PRs
> below have landed. Some files will still read only in `AlignIO`: the rare
> inputs the table marks out of scope, FASTA `-m 10` output (`fasta-m10`,
> a format only `AlignIO` has), and files that `Bio.Align`'s stricter
> checks reject (`funny.sth`, `simple_short.xmfa`). All five have merged
> (see the status note above):
>
> 1. `align-parse-nonseekable` — `Bio.Align.parse` reads stdin, pipes and
>    other non-seekable streams.
> 2. `align-stockholm-gr` — per-residue `#=GR` annotations drop each row's
>    gap columns, instead of every `.`.
> 3. `align-stockholm-blocks` — interleaved (multi-block) Stockholm. After 2.
> 4. `align-phylip-multi` — PHYLIP files holding several alignments
>    (seqboot output).
> 5. `align-phylip-relaxed` — `phylip-relaxed` in `Bio.Align`: a reader and
>    writer for names of any length, separated from the sequence by
>    whitespace. Today's `phylip` cuts names at 10 characters
>    (`Bio/Align/phylip.py:39,98-99`). A name with a hyphen used to be
>    impossible, because `Bio.Align` imported `Bio.Align.<name>`. §1.2's
>    registry core (PR #149) removed that: `Bio.Align._load` now looks names
>    up in a registry (`Bio/Align/__init__.py:5319`), so a format name no
>    longer has to be a module name, and naming the new format takes one
>    registry entry. It waited on §1.2's public `register_format`
>    contract, which PR #146 decided.
>
> **Rules the code PRs follow:**
>
> - Seekable streams keep today's behaviour, including handles positioned
>   part-way into a file and `next()` followed by `for`.
> - `len()` of a non-seekable stream raises `TypeError`, so `list()` falls
>   back to plain iteration.
> - Each `#=GR` annotation must span the alignment, as each `#=GC`
>   annotation already must. In an interleaved file this is checked after
>   the blocks are joined, as `AlignIO` does.
> - Repeated Stockholm names must form whole blocks: each block lists every
>   sequence once. Anything else raises `ValueError`.
> - `Bio.Align` gets no `phylip-sequential` name. Its reader already detects
>   the layout, and its writer, one line per sequence, already writes valid
>   sequential PHYLIP.
>
> **Where the two stacks differ.** "Upstream row" means a row in the Queued
> table of `UPSTREAM.md`.
>
> | format | difference | outcome |
> |---|---|---|
> | all | Record defaults. `AlignIO` copies the id into `name`, `description` or both; `Bio.Align` leaves the defaults. MSF is the same in both. | kept |
> | Stockholm, EMBOSS, MAF, Mauve | Annotation key names. Stockholm `secondary_structure` and `GS:`/`GC:` prefixes against `secondary structure` and no prefix; EMBOSS `identity` against `Identity`. For MAF and Mauve, `AlignIO` keeps start, size and strand in each record's `annotations`, `Bio.Align` in `coordinates`. | kept |
> | Nexus, Stockholm | All-gap columns. `Alignment` cannot hold one, so `Bio.Align` drops them: `test_Nexus_input.nex` has 48 columns in `AlignIO` and 46 in `Bio.Align`. | kept |
> | MAF | `.` in an `s` line. `AlignIO` reads it as the reference letter, `Bio.Align` as a gap. `AlignIO`'s reading contradicts the line's size: in `Tests/MAF/humor.maf`, mm3 has size 3,424 but `AlignIO`'s row has 5,166 residues (`Bio.Align`'s has 3,424). | kept; reported upstream as a behaviour choice (upstream row 27) |
> | MAF | Size check. `Bio.Align` raises when an `s` line's size does not match its sequence (`Bio/Align/maf.py:472-476`); `AlignIO` does not check. | kept |
> | MAF | `Bio.Align` needs the `##maf` header line (`Bio/Align/maf.py:305-306`). `Tests/MAF/ucsc_mm9_chr10_big.maf` has none: 983 alignments in `AlignIO`, `ValueError` in `Bio.Align`. | out of scope (rare input) |
> | Mauve | Different id model. `AlignIO` ids are `1/0-5670` (sequence number and range); `Bio.Align` ids are `0`, `1`, with the file name in the description. `Bio.Align` also checks each sequence against its stated range, so `simple_short.xmfa`, whose sequences are cut short, raises `ValueError` there and reads in `AlignIO`. | kept |
> | EMBOSS | Reverse strand (`-sreverse`). `AlignIO` raises on `Tests/Emboss/water_reverse1-4.txt` ([#1376](https://github.com/biopython/biopython/issues/1376), open since 2017); `Bio.Align` reads them. | out of scope (an `AlignIO` fix) |
> | MSF | Strictness and error messages differ. `Tests/msf/DOA_prot.msf`, whose header says length 62 where the sequences have 250, reads only in `Bio.Align` (with a warning). | kept |
> | Clustal | Several alignments in one file. `AlignIO` reads them; `Bio.Align` raises at the second header. | out of scope (rare input) |
> | PHYLIP | Several alignments in one file, such as seqboot output or `Align.write` given several alignments. `AlignIO` reads them; `Bio.Align` raises. | fixed in PR #184; upstream row 26 |
> | PHYLIP | Relaxed PHYLIP. `Bio.Align` has no `phylip-relaxed`. | fixed in PR #195 |
> | PHYLIP | `AlignIO` ignores the header's length (`Bio/AlignIO/PhylipIO.py:206`), so `phylip` and `phylip-relaxed` can misread silently: relaxed short names read under `phylip` give rows of 3 residues with ids `s1 AAAAAAA`, and names repeated in a later block end up inside the sequence. `Bio.Align` raises on both. Every silent misread found was `AlignIO`'s, so `Bio.Align`'s auto-detection is not changed to match. | reported upstream, not fixed here (upstream row 29) |
> | PHYLIP | Interleaved without blank lines between blocks. `AlignIO` reads it; `Bio.Align` raises (it does not misread). | out of scope |
> | Stockholm | Interleaved blocks. `Bio.Align` makes one row per sequence line: `simple.sth` raises, an equal-width two-block file silently reads as 4×5 instead of 2×10, and duplicate names become extra rows. | fixed in PR #197; upstream row 25 |
> | Stockholm | `#=GR`. `Bio.Align` deletes every `.` from the annotation, not the row's gap columns (`Bio/Align/stockholm.py:299`). The annotation shifts (`example_nonstandardannotations.sth`), and `Bio.Align` cannot read back its own output when an annotation has a `.` at a residue (`TypeError`). | fixed in PR #183; upstream row 28 |
> | Stockholm | `#=GF SQ` must match the number of sequences in `Bio.Align` (`Bio/Align/stockholm.py:233`), so `funny.sth` raises; `AlignIO` ignores it. Deliberate ([#3747](https://github.com/biopython/biopython/issues/3747)). | kept |
> | all | Non-seekable streams. `Bio.Align.parse` seeks, so stdin and pipes raise `io.UnsupportedOperation`; `AlignIO` reads them. | fixed in PR #182; upstream row 24 |
> | FASTA `-m 10` | Only `AlignIO` has `fasta-m10` (`Bio/AlignIO/__init__.py:165`): `Tests/Fasta/output001.m10` reads as 4 alignments there. `Bio.Align` has no such format. Its `tabular` reads FASTA's `-m 8CB` and `-m 8CC` output, not `-m 10`, and raises `ValueError` on that file. | out of scope (not one of the eight formats both stacks implement) |
> | writers | Only nexus output is byte-identical. Clustal, PHYLIP, Stockholm, MAF and Mauve output differ, and `AlignIO` cannot write EMBOSS or MSF. | kept: both sets of writers stay |

Eight formats implemented twice — ~3,300 lines in `Bio/AlignIO/` duplicating
~2,450 in `Bio/Align/`: clustal, emboss, msf, nexus, phylip, stockholm, maf,
mauve. Neither is deprecated, yet `DEPRECATED.rst:175-177` steers users to
`Bio.Align`.

The migration target is a feature *subset*: `Bio.Align` has no
`phylip-relaxed`/`phylip-sequential` (`Bio/AlignIO/PhylipIO.py:291,353`), no
`MafIndex` (`Bio/AlignIO/MafIO.py:257` — a bgzip-aware SQLite interval index),
and no `convert()`, `index()` or `to_dict()`. The two PHYLIP parsers also behave
differently: `Bio/Align/phylip.py:134-143` auto-detects interleaved vs
sequential, `AlignIO` requires you to name the variant — so the same file can
parse differently under the two APIs.

**Plan:** close the gaps in `Bio.Align` (port `MafIndex`, add the PHYLIP
variants, add `convert`); migrate the 21 internal `Bio/` consumers off
`AlignIO`; then reimplement `Bio/AlignIO/__init__.py` as a thin shim wrapping
`Alignment` → `MultipleSeqAlignment` and delete the duplicate format modules.
Do **not** remove `AlignIO.read`/`parse` — they are among the most-used entry
points in the ecosystem. Existing `Tests/test_AlignIO_*.py` passing unchanged is
the acceptance gate.
**Effort L · Impact high**

### 1.2 Five incompatible format-registration mechanisms, none extensible **[FIXED]**

> **Status: done, as decided below.** PR #149 added a private, lazy
> `FormatRegistry` (`Bio/_io_registry.py`). SeqIO and `Bio.Align` look format
> names up through it, and `Bio.Phylo` has too since PR #174, so a
> `Bio.Align` format name no longer has to be a module name. PR #185 added
> `register_format()` to SeqIO and `Bio.Align`, and PR #198 the entry-point
> plugins, both following the "Decided API" list below, which PR #146
> settled. The entry-point groups are `biopaithon.*`. AlignIO and SearchIO
> are left as they are, by decision. The NumPy consequence this section
> first reported went with PR #73: `import Bio.SeqIO` loads no format module,
> nor NumPy, `Bio.Align`, `Bio.AlignIO`, `urllib.request` or `xml.sax`, and
> FASTA parses with NumPy uninstalled (re-checked on 2026-10-08).
> `test_SeqIO.LazyFormatRegistries.test_import_seqio_is_lazy` and
> `test_io_registry.BuiltinTables.test_import_seqio_stays_light` guard it,
> and `test_register_format.StaysWithoutImportlibMetadata` guards that the
> built-in formats never look for entry points.

`Bio/SeqIO/__init__.py:575-675` holds lazy `"Module.Class"` strings,
`Bio/Align/__init__.py:5327-5334` builds its table from the `formats` tuple, and
`Bio/Phylo/_io.py:21-38` holds `"package.module"` specs. All three are
`FormatRegistry` tables, resolved on first use (PRs #149 and #174).
`Bio/AlignIO/__init__.py:162-185` uses dicts of eagerly imported classes.
`Bio/SearchIO/__init__.py:210-261` holds lazy
`(module, class)` string tuples, which `get_processor`
(`Bio/SearchIO/_utils.py:35-64`) imports on use.

Consequences: until PR #149, `Bio.Align` format names had to be valid Python
module names. That is *why* it could not offer `phylip-relaxed` (§1.1), and why
it says `tabular` where SeqIO says `fasta-m10`. Case handling differs —
`Align.read(f, "FASTA")` works, `SeqIO.parse(f, "FASTA")` raises. Until
PR #185 there was no registration hook at all, so a downstream package had to
mutate private dicts. hybran, LMAT, recentrifuge and psico assign into
`SeqIO._FormatToIterator`, and tfbayes into AlignIO's table.
Until PR #73 the eager imports also made NumPy a requirement for parsing FASTA
(see the status note). `import Bio.AlignIO` still pulls in `Bio.Align` and
NumPy.

**Plan:** one private, lazy `FormatRegistry` (built in PR #149), a `dict`
subclass generalising SeqIO's old `_LazyFormatRegistry`. It holds
`"package.module"` or `"package.module:attr"` specs, the form entry points use,
keyed on an explicit name so subtype names stay free-form. SeqIO and `Bio.Align`
resolve through it, and so does `Bio.Phylo` (PR #174). Each keeps its existing
table name, case rules and error messages. AlignIO and SearchIO are
deliberately left as they are. Lazy AlignIO tables would save little: most of
`import Bio.AlignIO` is `Bio.Align` and NumPy, pulled in by its module-level
`MultipleSeqAlignment` import, and its format modules add a few milliseconds.
SearchIO is already lazy, with no known downstream use of its maps. Add a
public `register_format()` for SeqIO and `Bio.Align`, plus `importlib.metadata`
entry-point groups so plugins work.

**Decided API:**

1. **Entry-point groups:** `biopaithon.seqio.iterators`,
   `biopaithon.seqio.writers` and `biopaithon.align`, after this fork's
   distribution name. They claim nothing in upstream's namespace. Adding
   `biopython.*` groups later would break nothing, while starting with them
   and then retreating would break plugins.
2. **Signatures:** `Bio.SeqIO.register_format(name, iterator=None, writer=None,
   *, replace=False)` and `Bio.Align.register_format(name, module, *,
   replace=False)`. A SeqIO `iterator` or `writer` is a callable or a
   `"package.module:attr"` string, and at least one must be given. The Align
   `module` is a module, an object with an `AlignmentIterator` attribute
   (`AlignmentWriter` is optional), or a `"package.module"` or
   `"package.module:attr"` string. AlignIO, SearchIO and `Bio.Phylo` get no
   public hook.
3. **Name rules follow each package's own lookup.** A SeqIO name must pass
   the name checks `SeqIO.parse` makes before its table lookup (a non-empty
   string for which `str.islower()` is true), and fails with the same
   `TypeError` or `ValueError` message. A `Bio.Align` name is any non-empty
   string, stored lowercased because `Bio.Align` lowercases on lookup; a
   non-string name raises `TypeError`. Entry-point names follow the same
   rules: an invalid SeqIO name is skipped with a warning, and an Align name
   is lowercased.
4. **Precedence comes from where a name came from, not when:** built-in, then
   `register_format`, then entry point. `replace=True` is needed to replace a
   built-in or another `register_format` call's entry. Registering the same
   value again is a no-op. A SeqIO call that sets both iterator and writer is
   all-or-nothing. `register_format` overrides an entry-point plugin without
   `replace=True`. A plugin naming a built-in is skipped with a warning. If two
   installed distributions offer one name with different targets, neither is
   registered and one warning names both.
5. **What `replace=True` does in SeqIO:**
   - It drops the `SeqIO.convert` fast paths that use the replaced role, so
     `convert` cannot silently ignore the override.
   - For the same reason, replacing an iterator drops the `sff`, `sff-trim`
     and `uniprot-xml` index proxies, which parse with the built-in parser
     whatever the table holds. `SeqIO.index` and `index_db` then refuse those
     names. Replacing only a writer leaves every index proxy as it is.
   - It keeps every other built-in index proxy, as assigning into the
     private dict does today. The proxy finds record boundaries and keys, and
     the replacement parses each record from a text handle. So the
     replacement must read text, and give each record the id the proxy's key
     rule gives, or a lookup raises `ValueError` (`StreamModeError` or "Key
     did not match"). For ace, fasta, phd, pir and qual that rule is the
     replacement's `parse_id_from_header`; every other kept proxy uses its
     own built-in rule.
   - The proxies that find records by a start marker (ace, embl, fasta, gb,
     genbank, imgt, phd, pir, qual, swiss) read `record_start_marker` from
     the replacement. A subclass of the built-in iterator inherits it;
     without it, `SeqIO.index` raises.
   - An `index_db` database keeps the keys it was built with. Rebuild one
     made before a replacement that changes ids.
6. **`SeqIO.index` for a new name** works when the iterator is a
   `SequenceIterator` subclass that has `"t"` in its `modes`, sets
   `record_start_marker`, and overrides `parse_id_from_header` (a public hook
   since PR #133). Built-in names with no index proxy stay unsupported, and
   the proxy classes stay private.
7. **`Bio.Align.formats` lists built-in names only,** and stays unchanged.
   Registered and plugin names never appear in it.
8. **Discovery is lazy and runs once per process.** Entry points are scanned
   on a table's first miss or first listing (iteration, `keys()`, `len()`,
   `values()`, `items()`), never on a hit, so `import Bio.SeqIO` still does not
   import `importlib.metadata`. There is no rescan after a package is installed
   into a running process. A failed scan warns once and means "no plugins"; a
   plugin's own import error is not caught.
9. **Registration is per process.** Register at module top level, or use an
   entry point. Multiprocessing's spawn and forkserver workers re-run module
   top-level code and rediscover entry points, but they never see a
   registration made inside `if __name__ == "__main__":` or a notebook cell.

**Effort L · Impact high**

### 1.3 The typing story is worse than having no types at all **[harness and the sequence spine FIXED, I/O and alignment open]**

> **Status: about half done.** PR #152 added the harness: the
> `disallow_untyped_defs` allowlist in `.mypy.ini` proposed below, a test
> of the config's shape, a full-tree mypy hook and a `mypy --strict`
> check of downstream code (`Tests/downstream_typing/`). Then these
> modules were fully annotated and joined the allowlist: `Bio.Seq`
> (PR #176), `Bio.File` (PR #177, its index classes now generic in the
> record type), `Bio.SeqRecord` (PR #187), `Bio.SeqFeature` (PR #188) and
> `Bio.Cluster` (PR #192). `from Bio.Restriction import EcoRI` no longer
> fails type checking (PR #179): the enzymes type as `Any` until they are
> built lazily (§1.7). The C extension stubs are §1.8. Four planned steps
> were not written, and are the rest of this item:
> `annotate-seqio-alignio` (the `Bio.SeqIO` and `Bio.AlignIO` entry
> points), `annotate-align-msa-io` (`MultipleSeqAlignment`,
> `Bio.Align.read`/`parse`/`write` and `Bio/Align/interfaces.py`),
> `annotate-alignment` (`Alignment` and its containers) and
> `annotate-aligners` (`PairwiseAligner`, `CodonAligner` and the
> convenience functions, after which `Bio.Align` joins the allowlist).
> Docstring styles are converted only in modules being annotated anyway
> (see addendum B and `AGENTS.md`).

`Bio/py.typed` and `BioSQL/py.typed` ship (`MANIFEST.in:66-68`), which is a PEP
561 promise to type checkers. An AST census of all 298 modules found **296 of
4,787 definitions fully annotated (6.2%)**. `Bio/SeqFeature.py` 0/96,
`Bio/Align/__init__.py` 2/92, `Bio/Seq.py` 4/136; `Bio/SearchIO` 0/352,
`Bio/Graphics` 0/229, `Bio/KEGG` 0/156, `Bio/motifs` 0/130.

Measured downstream consequence: `r: int = s.translate()` where `s = Seq("ACGT")`
produces **no mypy error**, and `reveal_type(s.reverse_complement())` is `Any`.
Shipping `py.typed` actively *suppresses* the errors a checker would otherwise
raise against an unstubbed package — users get false confidence.

**Plan:** annotate the spine in user-reach order — `Bio/Seq.py`,
`Bio/SeqRecord.py` (already 31/37, finish it), `Bio/SeqFeature.py`,
`Bio/SeqIO/__init__.py`, `Bio/AlignIO/__init__.py`, `Bio/Align/__init__.py` —
with a per-module `disallow_untyped_defs` allowlist in `.mypy.ini` that grows as
each is finished, locking in the gains.
**Effort L overall, S–M per module · Impact high**

### 1.4 The indexing layer re-implements each parser's record-boundary logic **[FIXED]**

> **Status: done, essentially as planned.** PR #88 moved the record markers
> and id rules onto the format's own parser classes, so `Bio.SeqIO.index`
> and `index_db` derive each key from the same code that sets `record.id`,
> and the corpus-wide test proposed below exists: index keys are asserted
> equal to parsed ids for every indexable format across the fixtures. Two
> formats changed behaviour, both converging on their parser: "pir" ids
> keep internal spaces, and pretty-printed UniProt XML with indented
> `<entry>` elements indexes instead of silently returning an empty index.
> The three unrelated random-access implementations named below still exist
> (`Bio/File.py`, `MafIO`, `bigbed`) — that consolidation is §1.1/§1.2
> territory.

`Bio/SeqIO/_index.py:187-201` hardcodes byte markers for *other modules'*
formats inside a base class, and `:219` re-derives record IDs independently —
the direct cause of §0.7. Only 19 of ~40 formats are indexable, and there are
three further unrelated random-access implementations
(`Bio/File.py:150,267`, `Bio/AlignIO/MafIO.py:257`, `Bio/Align/bigbed.py:1487`).

**Plan:** move the marker and ID-extraction rules onto each format's own parser
class (a `record_start_marker` attribute and a `parse_id_from_header`
classmethod on `SequenceIterator`), and have the random-access base read them
from `_FormatToIterator` — deleting the marker dict and making index support
automatic for any sequential format. Add a corpus-wide test asserting
`list(index(f, fmt)) == [r.id for r in parse(f, fmt)]`.
**Effort M · Impact high**

### 1.5 Three parallel BLAST stacks, and the tutorial teaches the superseded one **[a and b FIXED, c not done]**

> **Status: (a) and (b) are done.** PR #90 rewrote the Tutorial's BLAST
> chapter against `Bio.Blast`, including a "Migrating from the older BLAST
> modules" section mapping the old record attributes onto the new classes,
> and PR #91 added `BiopythonDeprecationWarning`s to `NCBIWWW` and `NCBIXML`
> with the `DEPRECATED.rst` entry — deprecation only, the modules still
> work, and removal waits at least a release as the plan requires.
> (c) — whether SearchIO's `blast-xml` becomes an adapter over `Bio.Blast`
> — was not done in the October 2026 sprint. It stays recorded here as an
> option, with no plan to do it.

New: `Bio/Blast/__init__.py` + `_parser.py` + `_writers.py` (~4,100 lines).
Old and carrying **no deprecation warning**: `Bio/Blast/NCBIXML.py` (1,331) and
`NCBIWWW.py` (373). Third: `Bio/SearchIO/BlastIO/`, a fourth independent XML
implementation. `Doc/Tutorial/chapter_blast.rst:1536,1574,1825` still teaches
`NCBIWWW.qblast` and `NCBIXML.read`.

**Plan, in order:** (a) rewrite the tutorial against `Bio.Blast` — high value,
breaks nothing; (b) add deprecation warnings to `NCBIXML`/`NCBIWWW` with a
`DEPRECATED.rst` entry; (c) decide separately whether SearchIO's `blast-xml`
becomes an adapter. Do not remove `NCBIXML` in the same release as the warning.
**Effort L (M if scoped to a+b) · Impact high**

### 1.6 `check_untyped_defs` is off, hiding 2,172 real errors **[FIXED — ratchet in place]**

> **Status: enabled, with the ratchet.** PR #96 turned
> `check_untyped_defs = True` on globally, with a per-module
> `check_untyped_defs = False` baseline in `.mypy.ini` for the modules that
> do not yet pass — shrink-only: a module leaves the baseline when it is
> fixed and nothing may join it. The approach paid for itself immediately:
> working through the baseline surfaced `Bio.Pathway.System.stochiometry`,
> broken since its introduction in 2001 (`[] * len(reactions)` is always
> empty), fixed in PR #97 with its first-ever regression tests, and its
> module left the baseline.

`.mypy.ini:6` has `#check_untyped_defs = True` commented out, so mypy skips
essentially every function body. Enabling it: **2,172 errors in 186 files**.
By category: `attr-defined` 791, `assignment` 331, `union-attr` 329,
`var-annotated` 162, `index` 147. Worst: `Bio/Align/bigbed.py` 148,
`Bio/Seq.py` 120, `Bio/Blast/_parser.py` 118, `Bio/GenBank/Scanner.py` 55.

The 329 `union-attr` and 791 `attr-defined` hits are exactly the shape of latent
`AttributeError`/`None`-dereference bugs in parser paths that only fire on
unusual input — the category this library gets issues about.

**Plan:** enable globally, then add a mechanically generated per-module
`ignore_errors` baseline for the 186 failing modules as a ratchet, so new and
touched code is checked from day one and the list only shrinks. Triage
`Bio/GenBank/Scanner.py` and `Bio/Blast/_parser.py` first.
**Effort M · Impact high**

### 1.7 Import cost: ~190 ms warm for `import Bio.SeqIO`, and 77 ms of it is one module body **[mostly FIXED]**

> **Status: the `CodonTable` and `SeqIO` pieces are done.** PR #53 made the
> ambiguous codon list expansion lazy — expanded on first use instead of at
> import — and PR #73 made the per-format `SeqIO` parser imports lazy:
> `import Bio.SeqIO` is roughly ten times faster (about 107 ms down to
> 11 ms, median of five `python -X importtime` runs on 3.12), submodules
> stay reachable as `Bio.SeqIO.FastaIO` and friends, and parsing FASTA now
> works on a machine without NumPy. The `Restriction` PEP 562 `__getattr__`
> is deliberately deferred on the usage evidence (the survey later in this
> document found `Bio.Restriction` far less imported than assumed);
> revisit only if import-cost complaints actually name it. That deferral
> was confirmed in October 2026. Its typing half no longer waits on it:
> PR #179 stopped `from Bio.Restriction import EcoRI` failing under mypy,
> with the enzymes typed as `Any` until they are built lazily.

Measured cumulative `-X importtime`: `Bio.Seq` 45.3 ms, `Bio.Align` 155.9 ms,
`Bio.Restriction` 239.6 ms, `Bio.SeqIO` 547.8 ms cold / 188.7 ms warm.

Three specific causes:
- `Bio/Restriction/Restriction.py:2600-2637` synthesises ~1,000 enzyme classes
  in a module-level loop — **77.5 ms of self-time**, paid even to use one enzyme.
  It is also invisible to type checkers: `from Bio.Restriction import EcoRI`
  raises `attr-defined` under plain, non-strict mypy today.
- `Bio/Data/CodonTable.py:603-1300` builds ~34 codon tables at import — a ~29 ms
  module body, imported from `Bio/Seq.py:30`, so *every* entry point pays it.
- `Bio.SeqIO` eagerly pulls `numpy` (101 ms), `urllib.request` (19 ms) and
  `xml.sax.saxutils` (22 ms).

This is the most-felt performance characteristic of the library — startup
latency on every CLI tool, notebook and workflow task rule — and it is unrelated
to parser speed.

**Plan:** lazy table construction in `CodonTable` behind a module `__getattr__`;
PEP 562 `__getattr__` in `Restriction` building enzymes on first access (plus a
generated `.pyi` declaring the names); lazy per-format imports in `SeqIO`. Add a
CI budget check on `python -X importtime -c "import Bio.SeqIO"`.
**Effort M · Impact high**

### 1.8 No `__all__` on re-exporting packages, and no stubs for 13 C extensions **[FIXED]**

> **Status: both halves are done.** PR #178 added `.pyi` stubs for the
> `Bio.Align` C extensions, with `stubtest` in CI, so `PairwiseAligner`
> is no longer `Any` to a type checker. PR #192 stubbed `Bio.Cluster._cluster`
> and PR #191 the rest, and made a stub compulsory for every extension.
> The `__all__` note below is from when only that half was done.
>
> **Status: the `__all__` half is done; the stubs half is not.** Measured
> rather than estimated: a file importing every public re-exported name of
> every package, checked with `mypy --strict --follow-imports=silent`, gave
> **175** "does not explicitly export" errors across **30** packages — more
> than the 13 guessed at below. All 30 now declare `__all__` and the count is
> 0. `Tests/test_public_exports.py` asserts each `__all__` resolves and still
> covers its package's public namespace, so a new re-export cannot quietly
> reintroduce the problem.
>
> Two corrections to the text below. **`Bio.SeqIO` was only partly affected** —
> `parse`, `read` and `write` are defined in `Bio/SeqIO/__init__.py` rather
> than re-exported, so they never had the problem; 9 of its names did.
> **`Bio.Restriction` is not fixable this way and is excluded**: its ~1,000
> enzyme classes are synthesised at import, so mypy cannot see them at all and
> reports *zero* re-export errors for it. It needs the generated stubs, which
> is the other half of this item.
>
> Note what this does *not* buy. `from Bio.PDB import PDBParser` now resolves,
> but `PDBParser()` still fails `--strict` with `Call to untyped function
> "PDBParser" in typed context`. That is §1.3, and it is the larger job.

Only 21 of 298 files define `__all__`, 13 of them inside `Bio/SearchIO`.
`Bio/PDB/__init__.py` re-exports 36 names with none. Since `.mypy.ini:8` sets
`no_implicit_reexport = True`, the project has accepted these semantics
internally but never declared its exports — so `from Bio.PDB import PDBParser`,
the exact import the tutorial teaches, fails under downstream `mypy --strict`.

Separately, `find Bio BioSQL -name "*.pyi"` returns nothing for 13 C extensions,
each imported with a blanket `# type: ignore`. Because
`Bio/Align/__init__.py:4158` is `class PairwiseAligner(_pairwisealigner.PairwiseAligner)`,
the base class is `Any` and the library's headline API is fully opaque:
`a.mode = 12345` type-checks clean despite runtime validation at
`Bio/Align/__init__.py:4402-4420` that rejects it.

**Plan:** add `__all__` to the 13 re-exporting `__init__.py` files, starting
with `Bio/PDB`, `Bio/SeqIO`, `Bio/AlignIO`; add a test asserting `__all__`
matches the documented API. Hand-write `.pyi` stubs starting with
`Bio/Align/_pairwisealigner.pyi` and `Bio/Cluster/_cluster.pyi`.
**Effort S (`__all__`) / M (stubs) · Impact high**

### 1.9 No C extension ever releases the GIL, so threads give zero speed-up **[FIXED, except Ctrl-C in `Bio.Cluster` and `ccealign`]**

> **Status: the `cluster.c` slice is done too.** PR #196 releases the GIL
> around all nine `Bio.Cluster` kernel calls; `cluster.c` now allocates
> with `PyMem_Raw*` so it can run without the GIL. What is left is the
> plan's Ctrl-C half for those kernels: the planned step
> `cluster-interrupt` (a periodic signal check in the clustering loops)
> was not written, so Ctrl-C still waits for a long clustering call to
> finish. `ccealign`'s kernels still hold the GIL. The note below is from
> before PR #196.
>
> **Status: the two hottest extensions now release the GIL.** PR #75 wrapped
> the NW/SW/Gotoh pairwise alignment kernels and PR #74 the `kdtrees` build
> and search kernels (measured there: two threads on separate
> 100,000-point datasets run at about 1.75× serial speed, previously 1.0×).
> PR #74 also made `kdtrees` reentrant — searching one shared tree from
> several threads is safe, and `neighbor_simple_search` no longer corrupts
> the tree by re-sorting the shared point list. The `cluster.c` kernels are
> the remaining slice; PR #89's removal of the file-scope RNG state was its
> prerequisite. `ccealign`'s kernels stay GIL-bound for now. On the
> free-threaded build of Python 3.14 there is no GIL to release: every
> extension now declares it can run without one (§1.10), so threads run
> all of them in parallel there, `Bio.Cluster` and `ccealign` included.

`grep -rn "Py_BEGIN_ALLOW_THREADS" --include=*.c .` returns **nothing** across
all 13 extensions. Measured on a 14-core machine,
`PairwiseAligner(scoring="blastn").score()` on two random 20,000 nt sequences:

| | time |
|---|---|
| 1 alignment | 1.41 s |
| 4 sequential | 5.64 s |
| **4 threads** | **5.72 s (speed-up 0.98×)** |

`SIGINT` delivered 0.5 s into a 1.4 s `score()` call raised `KeyboardInterrupt`
only at 1.54 s — Ctrl-C is dead for the entire duration of any C call, which for
a long alignment means hours.

`ThreadPoolExecutor` over alignments, clustering or KD-tree searches therefore
gains literally nothing, forcing users to `multiprocessing` and re-pickling.

**Plan:** wrap the pure-C kernels that touch no Python objects in
`Py_BEGIN_ALLOW_THREADS`/`Py_END_ALLOW_THREADS` — the `*_align`/`*_score` macro
bodies in `_pairwisealigner.c:4593-5300`, `findPath`/`calcS`/`calcDM` in
`ccealignmodule.c` after `getCoords` returns, `KDTree_build_tree` and
`KDTree_neighbor_search` in `kdtrees.c`, and the `cluster.c` kernels. Add
periodic `PyErr_CheckSignals()` in each outer loop so Ctrl-C works.
**Effort M · Impact high**

### 1.10 Free-threaded CPython is unsupported and blocked by mutable file-scope state **[FIXED — without multi-phase init]**

> **Status: the named static-state blockers are gone.** PR #74 threads the
> kdtrees sort dimension through per-call state instead of
> `DataPoint_current_dim`, and PR #89 replaced `Bio.Cluster`'s
> `srand(time(0))`-seeded file-scope RNG with caller-owned xoshiro256++
> state — which, exactly as predicted below, also delivered a reproducible
> `rng_seed=` keyword for `kcluster`/`kmedoids`/`somcluster` and stopped
> the library silently resetting the process-wide `rand()` stream (a
> regression test now checks `Bio.Cluster` never touches `rand()`).
>
> **Free threading is supported.** Every extension now declares, in its
> single-phase `PyInit_*`, that it can run without the GIL
> (`PyUnstable_Module_SetGIL` through `Bio/_freethreading.h`, which compiles
> to nothing on GIL builds and PyPy): PR #171 for the six that needed no
> other change, PR #189 for `Bio.Cluster` and PR #190 for the five
> `Bio.Align` ones, which lock the per-object state other calls could free.
> PRs #173 and #181 first fixed, on every build, the bugs that made a
> shared aligner unsafe: each call now snapshots the aligner's settings and
> checks private copies of the sequence indices. The 3.14t CI job runs the
> suite with no `PYTHON_GIL` override, and the release workflow builds
> `cp314t` wheels, whose test fails if importing the extensions turns the
> GIL on (PR #201).
>
> **The multi-phase init migration planned below was rejected.** Free
> threading does not need it: CPython's free-threading HOWTO gives
> `PyUnstable_Module_SetGIL` as the route for single-phase modules. The
> file-scope state left (type pointers and static type objects) is written
> once at import, under the import lock, and only read after, so it is no
> data race. The real bugs were per-object, and module state fixes none of
> them. Moving `_pairwisealigner`'s `Array_Type` into module state would
> also mean heap types for the aligner, rewriting upstream's most-edited C
> file and changing behaviour users can see. Module state is for isolating
> subinterpreters; if that is ever wanted, it is an item of its own.
>
> One step planned under this item was not written: `cluster-interrupt`,
> which would let Ctrl-C stop a long `Bio.Cluster` kernel (see §1.9).
> PR #196 (`cluster-release-gil`) and PR #201 (`ft-release`) were the
> last steps that merged.

No `Py_mod_gil`, `Py_MOD_GIL_NOT_USED` or `Py_GIL_DISABLED` anywhere; all 13
modules use single-phase init with `m_size = -1`. A single-phase module without a
`Py_mod_gil` slot causes the free-threaded interpreter to **re-enable the GIL
process-wide at import** — so one `import Bio` silently disables free-threading
for the whole application.

Concrete blockers: `Bio/PDB/kdtrees.c:20` `static int DataPoint_current_dim`,
written by `DataPoint_sort` and read by the `qsort` comparator, so two concurrent
`KDTree()` builds corrupt each other's sort key; `Bio/Cluster/cluster.c:64`
`static int TEMP_SWAP_INT` used as the temporary inside the `swap_int` macro in
the sort hot path; `cluster.c:73` and `:2030-2038` static RNG state seeded with
`srand(time(0))`; `Bio/Align/_pairwisealigner.c:40` `static PyTypeObject *Array_Type`.

`pyproject.toml` already advertises 3.13/3.14, so users will reach for `3.14t`,
and the fix window is now — before wheels start shipping `cp314t`.

**Plan:** migrate to multi-phase init with module state; make `TEMP_SWAP_INT` a
local; thread the sort dimension through a `qsort_r` context; replace the
`srand`-based RNG with caller-supplied state (which also finally makes
`Bio.Cluster.kcluster` reproducible and stops the library perturbing the
process-wide `rand()` stream); then add the `Py_mod_gil` slot and a `3.14t` CI job.
**Effort L · Impact high**

### 1.11 FASTQ parsing costs 4.9× the time and 3.9× the memory of the raw iterator **[FIXED]**

> **Status: part (b) is done too.** PR #167 gave the three FASTQ parsers
> a keyword-only `compact=False`. With `compact=True` each read's
> qualities are an `array.array("b")`: about 5× smaller for the
> qualities alone, and 514 MB down to 312 MB for 200,000 reads held in
> a list, as the record objects themselves do not shrink. The default
> output is unchanged. The note below is from before PR #167.
>
> **Status: part (a) is done.** PR #17 moved quality decoding into an
> overridable `_decode_quality` method: the base class returns
> `list(byte_scores)` and only `FastqSolexaIterator` keeps the signed
> `array.array("b", ...).tolist()` round-trip.
> `Tests/test_SeqIO_QualityIO.py` covers the full valid byte range for all
> three variants. Part (b) — opt-in `array('b')`/`bytes` backing for
> `letter_annotations` — remains open.

Measured on 200,000 × 150 bp reads (a 65 MB file):

| | time | peak RSS |
|---|---|---|
| `FastqGeneralIterator` | 0.23 s | 40 MB |
| `SeqIO.parse(..., "fastq")` | 1.12 s | 40 MB |
| `list(FastqGeneralIterator(...))` | 0.24 s | 140 MB |
| `list(SeqIO.parse(..., "fastq"))` | 2.36 s | **538 MB** |

A quarter of the time is one line — `Bio/SeqIO/QualityIO.py:1115`,
`array.array("b", byte_scores).tolist()`. Micro-benchmarked at this size:
`array('b', b).tolist()` takes 0.284 s versus `list(b)` at 0.119 s (**2.4×
faster**), and the resulting list is 5.2× the size of the array.

The signed round-trip exists only for Solexa. `FastqPhredIterator.q_mapping` and
`FastqIlluminaIterator.q_mapping` (`:1405-1412`) map every valid byte to 0–93 —
all non-negative — so `list(byte_scores)` is exactly equivalent there.

**Plan:** (a) move quality decoding into an overridable method and use
`list(byte_scores)` in the Phred and Illumina subclasses, keeping the array
round-trip only for Solexa — ~15% off `SeqIO.parse` with no API change; (b) offer
an opt-in `array('b')`/`bytes` backing for `letter_annotations`, which
`SeqRecord` already only requires to be a sized sequence, for a 5× memory
reduction.
**Effort S (a) / M (b) · Impact medium**

### 1.12 `PairwiseAligner.align()` allocates the full O(n·m) traceback matrix **[guards and Needleman-Wunsch FIXED, Gotoh and Smith-Waterman open]**

> **Status: linear space is in for Needleman-Wunsch.** PR #186 added a
> checkpointed linear-space traceback that yields the same paths in the
> same order as the full matrix (not Hirschberg's, which would pick a
> different co-optimal path). PR #199 counts the alignments in linear
> space too, so `len()`, `repr()` and `bool()` no longer build the
> matrix, and switched it on: Needleman-Wunsch alignments whose
> traceback matrix would exceed 512 MiB, about 23 kb × 23 kb and up,
> now take the linear path. At 20 kb, `align()[0]` takes about 1.5× the
> time in 67 MB instead of 415 MB. Three planned steps were not written,
> so Gotoh (affine gaps, which the `blastn`, `megablast` and `blastp`
> presets use), Smith-Waterman and the other modes still build the full
> matrix: `pairwise-linear-gotoh-global`, `pairwise-linear-sw` and
> `pairwise-linear-gotoh-local`.
>
> **Status: the short-term guards are in.** PR #67 computes the matrix size
> in checked `size_t` arithmetic before allocating: a size that cannot fit
> raises `MemoryError` naming the predicted size and both sequence lengths,
> a failed allocation reports the same numbers, and lengths above the
> aligner's `int` storage (`INT_MAX - 1`) raise `ValueError` — previously
> exactly `INT_MAX` slipped past into undefined behaviour, and FOGSAA
> computed its cell count in C `int`, so ~65,536-letter sequences wrapped
> the multiplication and corrupted the heap from pure Python. There is
> deliberately no cap below what the platform can address. The
> Hirschberg/Myers-Miller linear-space traceback remains the long-term
> item (see the status note above for what landed).

`Bio/Align/_pairwisealigner.c:69-80`. Fresh-process measurements with
`scoring="blastn"`:

| N | `score()` | `align()` |
|---|---|---|
| 10,000 | 0.36 s, 29.3 MB | 1.70 s, **220.1 MB** |
| 20,000 | 1.38 s, 29.5 MB | 6.63 s, **792.6 MB** |

`score()` is O(n) in memory; `align()` is exactly quadratic (2× N gives 4× RSS).
A 100 kb × 100 kb alignment would need roughly 20 GB. `nA`/`nB` are also stored
as `int` (`:72-73`), so lengths above `INT_MAX` are silently truncated rather
than rejected. This is the hard ceiling on `Bio.Align` for anything longer than a
gene, and it is invisible until the process is OOM-killed.

**Plan:** short term, validate the lengths and raise a clear error naming the
predicted allocation size. Longer term, add a Hirschberg/Myers-Miller
linear-space traceback for the single-best-path case (`align(...)[0]`, by far the
common use), keeping the full matrix only when enumerating all optimal
alignments. A banded mode is a cheaper intermediate win for near-identical
sequences.
**Effort S (validation) / L (Hirschberg) · Impact medium**

---

## Tier 2 — Cruft removal

A fork has more latitude here than upstream, but users are real. Where an item
is widely copy-pasted, prefer an `ImportError` stub with a migration pointer over
silent deletion.

### 2.1 Already-broken or already-empty — no user impact **[FIXED]**

> **Status: all four bullets are resolved.** `Scripts/xbbtools` is removed
> (PR #86 — dead at import since 1.86; `nextorf.py`, which shared the
> directory but not the broken import chain, moved to `Scripts/nextorf.py`).
> `import Bio.HMM` now raises an informative `ImportError` naming the
> removed modules and pointing at hmmlearn (PR #92, taking the
> `Bio.Alphabet` stub route argued for below; the package deliberately
> stays in `pyproject.toml` so the message ships). The unreachable
> `run_tests.py` doctest-exclusion block was already gone, and the two
> `FIXME remove this after 1.87` sites are resolved by PR #83, which
> removed the transitional `warn_defaults_changed()` (see
> `DEPRECATED.rst`).

- `Scripts/xbbtools/xbb_blastbg.py:20-24` and `xbb_blast.py:23,200` import
  `Bio.Blast.Applications`, **which no longer exists**, and `MANIFEST.in:14`
  ships it in every sdist. The reach is worse than "on first use" as written
  here originally: `xbb_widget.py:21` does `from xbb_blast import BlastIt` at
  module scope, and `xbbtools.py:17` imports `xbb_widget`, so `import xbbtools`
  raises `ModuleNotFoundError` before a window is ever drawn. The whole GUI is
  unreachable, not just its BLAST feature.
- `Bio/HMM/__init__.py` is a 5-line docstring-only file — all four submodules
  were removed in 1.86 — yet `pyproject.toml` still ships `Bio.HMM` as a package.
  Note `import Bio.HMM` currently *succeeds* and yields nothing, so removing it
  is a public-API change needing a `DEPRECATED.rst` entry, not a free deletion.
  `Bio/Alphabet/__init__.py` is the precedent for the alternative: keep the
  package and raise `ImportError` with a pointer.
- `Tests/run_tests.py:69-73` still excludes doctests for five modules removed in
  1.86, and `:96` for `Bio.PDB.Vector`, removed in **1.74**. The whole 45-entry
  list is guarded by `if np is None:` (`:61`) — unreachable, since numpy is a
  hard dependency. **[FIXED]** — the block and the `try: import numpy` that fed
  it are gone. A full offline run reports the same 520 modules before and after,
  which is what "unreachable" predicts.
- `Bio/Align/__init__.py:4449,4507`: `# FIXME remove this after 1.87 is out`.
  1.87 is out.

> **Correction.** This section originally also listed "`pyproject.toml:58`
> still ships `Bio.Alphabet` six years after removal". That is wrong, and it is
> the kind of wrong this document is meant to prevent. `Bio/Alphabet/__init__.py`
> is a deliberate stub whose body is `raise ImportError(...)` pointing at
> https://biopython.org/wiki/Alphabet. It has to be *packaged* for that error to
> reach anyone; dropping it from `pyproject.toml` would replace a helpful
> message with a bare `ModuleNotFoundError`. `Tests/run_tests.py` excludes it
> from doctests for the same reason. Leave it alone.

**Effort S · Impact medium** (nothing breaks — it is already broken)

### 2.2 `Bio.PDB.mmtf` — the format's server no longer resolves in DNS **[deprecation SHIPPED, removal scheduled]**

> **Status: step one is done.** PR #72 deprecated `Bio.PDB.mmtf` with a
> warning pointing at `Bio.PDB.binary_cif`, and the `DEPRECATED.rst` entry
> is explicit that this deprecation is BioPAIthon's own. Relatedly,
> `PDBList` no longer offers MMTF downloads and gained BinaryCIF (adopted
> upstream #4938), and the new `[structure]` extra deliberately omits
> `mmtf-python`. The removal itself — subpackage, test files, packaging
> entry, CI dependency — waits the customary release. **Decided:** that
> second step does not happen before the BioPAIthon 1.88 release, so
> 1.88 ships the deprecation and removal comes after it.

RCSB decommissioned MMTF in July 2024; `mmtf.rcsb.org` does not resolve. The
dependency `mmtf-python` last released 2022-07-06. `Tests/test_mmtf_online.py:28`
calls `get_structure_from_url("4ZHL")` and can never pass again. 563 lines plus a
CI dependency (`ci-dependencies.txt:18`).

**Plan:** one release with a deprecation warning, then remove the subpackage,
both test files, the `pyproject.toml` entry and the CI dependency; repoint
`Tests/test_PDB_internal_coords.py:63-66` at an existing mmCIF fixture.
`Bio.PDB.binary_cif` is the successor.
**Effort S · Impact medium-high**

### 2.3 `Bio.pairwise2` — deprecated eight releases ago, still built and shipped **[FIXED]**

> **Status: removed, exactly as planned.** PR #79 deleted the ~2,900 lines
> and ships `Bio/pairwise2.py` as a stub whose import raises `ImportError`
> naming `Bio.Align.PairwiseAligner` — and warning that migration is not
> mechanical: `pairwise2`'s default gap score is 0 where `PairwiseAligner`'s
> has been -1 since 1.86, so scores should be set explicitly when porting.
> This also removes the last consumer of `cpairwise2module.c` (§0.14's
> module) and its per-platform CI compile cost. Upstream Biopython 1.88
> still ships the deprecated module; this is a deliberate difference.

Deprecated in 1.80 (`Bio/pairwise2.py:274-284`); current version 1.88.dev0.
~2,900 lines: the module (1,441), `Bio/cpairwise2module.c` (479), three test
files, and a whole tutorial chapter. Zero internal consumers. The C extension
must compile on every platform/Python/PyPy combination in CI.

**Plan:** delete the implementation but leave a small `Bio/pairwise2.py` raising
`ImportError` with a `PairwiseAligner` migration pointer — mirroring the
`Bio/Alphabet/__init__.py:20` pattern. It is heavily copy-pasted in tutorials and
StackOverflow answers, and semantics differ from `PairwiseAligner` (notably the
1.86 gap-score default change), so migration is not mechanical and a silent
deletion would be hostile.
**Effort M · Impact high**

### 2.4 The 1.86 deprecation cohort is ripe for one batched removal **[FIXED]**

> **Status: removed in one batch** in PR #70 — the seven `as_*` writer
> helpers and `SummaryInfo`, whose module `Bio.Align.AlignInfo` (widely
> imported in public code) ships as an `ImportError` stub naming
> `msa.alignment` as the replacement. The `PairwiseAligner` alias table was
> held one more release, exactly as the plan below says.

All emit `BiopythonDeprecationWarning`, are documented in `DEPRECATED.rst`, and
have zero internal callers: `as_fasta`/`as_fasta_2line`
(`Bio/SeqIO/FastaIO.py:645,664`), `as_fastq`/`as_qual`/`as_fastq_solexa`/
`as_fastq_illumina` (`Bio/SeqIO/QualityIO.py:1652,1818,1920,1998`), `as_tab`
(`Bio/SeqIO/TabIO.py:136`), the `SummaryInfo` class
(`Bio/Align/AlignInfo.py:32-47`), and `PairwiseAligner.__setattr__`/`__getattr__`
with its 20-entry alias table (`Bio/Align/__init__.py:4392-4444`) — which sits on
the hot path of the most-used alignment class.

**Plan:** remove in one batch. Hold the `PairwiseAligner` alias table one more
release; the rename was only in 1.86 and those callers are more common.
**Effort S · Impact medium**

### 2.5 `Bio.codonalign` has warned "experimental" for twelve years **[FIXED]**

> **Status: decided and done.** The warning is dropped and the API is
> committed (PR #141): `Bio.codonalign` is now under the normal
> backwards-compatibility policy, so future changes go through
> `DEPRECATED.rst`. PR #154 then made `cal_dn_ds` and `mktest` delegate
> to `Bio.Align.analysis` (addendum C), removing about 1,230 lines. The
> alignment-building half was kept.

`Bio/codonalign/__init__.py:21-25` warns on every import; introduced 2014.
2,644 lines whose last substantive commit was a formatting pass. An import-time
warning that has fired unchanged for twelve years trains users to filter *all*
Biopython warnings, which then masks the real deprecations above.

**Plan:** decide and record. Either drop the warning and commit to the API, or
migrate alignment-building onto `Bio.Align.CodonAligner` and keep only
`cal_dn_ds` and `mktest`, which are what people actually use.
**Effort M · Impact medium**

### 2.6 Python-2-era shims, one of which silently breaks a documented promise **[FIXED, `colour` removal scheduled]**

> **Status: the `colour` aliases now warn.** PR #168 made each of the 15
> aliases (the 101 text matches below come down to 15) raise
> `BiopythonDeprecationWarning`, pointing at the caller's line. They are
> to be removed in the release after 1.88. The note below is from before
> PR #168.
>
> **Status: the `__nonzero__` half is fixed** in PR #12, which took the
> "real `__bool__`" branch of the decision below: both methods are now
> `__bool__` returning `True` per the documented promise, the docstrings are
> corrected, and truthiness of zero-length locations is tested. The
> `iteritems()` remnant was deprecated in PR #66 (deprecated rather than
> deleted, since upstream still ships it as public API). The 101 unwarned
> `colour` aliases remain fully open — add warnings, then remove a release
> later (tracked in `TODO.md`).

`Bio/SeqFeature.py:1117,1515` define `__nonzero__` — the **Python 2** name.
Neither class defines `__bool__`, but both define `__len__`, so zero-length
locations are falsy. The docstring at `:1118-1120` promises "Return True
regardless of the length of the feature… for backwards compatibility" — a
contract silently broken since Python 3 became the only runtime. Also
`Bio/SearchIO/_model/query.py:233` `iteritems()`, and 101 UK-spelling `colour`
aliases in `Bio/Graphics/GenomeDiagram/` deprecated in **1.55 (2010)** that emit
no warning at all.

**Plan:** delete both `__nonzero__` methods and consciously decide the
truthiness question — either add a real `__bool__` returning `True` per the
documented intent, or update the docstrings to state that zero-length locations
are falsy. The current state is wrong either way, so the decision must be
explicit. Add warnings to the `colour` aliases before removing them.
**Effort S · Impact low-medium**

### 2.7 Unhashable value objects **[FIXED]**

> **Status: fixed** in PR #48. All four classes now define `__hash__` in
> `Bio/SeqFeature.py`, with `SeqFeature`'s excluding the mutable
> `qualifiers` dict.

`SimpleLocation`, `CompoundLocation`, `SeqFeature` and `Reference` define
`__eq__` without `__hash__` (`Bio/SeqFeature.py:1197,1538,227,599`), so Python
sets `__hash__ = None`. `set(record.features)`, deduplicating features across
annotation sources, and using a location as a dict key all fail with an opaque
error. This is drift, not policy: `Bio/Seq.py:2163` deliberately restores
`Seq.__hash__` for exactly this reason, and `UnknownPosition` defines one.

**Fix:** add `__hash__` to all four, excluding the mutable `qualifiers` dict from
`SeqFeature`'s.
**Effort S · Impact medium**

---

## Tier 3 — Build, packaging and CI

### 3.1 Supply-chain exposure and a shell-injection sink **[mostly FIXED]**

> **Status: the pinning and the injection sink are fixed** in PR #18. Every
> action in `ci.yml` — including `tj-actions/changed-files` — is SHA-pinned
> with a version comment, and the changed-files output now reaches
> pre-commit through an `env:` variable consumed as JSON by a Python step,
> not shell interpolation. **Still open:** `.github/dependabot.yml` covers
> only `github-actions`, with no `pip` ecosystem, and the `test_*` jobs
> still carry `secrets.CODECOV_TOKEN`.

`.github/workflows/ci.yml:30` uses `tj-actions/changed-files@v47` — a floating
tag on the action compromised in CVE-2025-30066 — and `:56` interpolates its
output unquoted into a shell body:
`pre-commit run --files ${{ steps.changed-files.outputs.all_changed_files }}`,
with attacker-controlled filenames from a PR branch. Every action is tag-pinned
rather than SHA-pinned, and `test_*` jobs carry `secrets.CODECOV_TOKEN`
(`:218,264,310`). `.github/dependabot.yml` covers only `github-actions`, not pip.

**Fix:** SHA-pin every third-party action with a version comment; move the file
list into an `env:` var and quote it, or drop `tj-actions/changed-files` and run
`pre-commit run --all-files`. Add a `pip` ecosystem to Dependabot.
**Effort S · Impact high**

### 3.2 CI builds 15 wheels per run; none is installable, tested, or published **[FIXED]**

> **Status: fixed.** PR #4 (`96aba7bda`) removed `build_wheels` and the dead
> `cleanup_wheels` job. `release.yml` now builds real wheels with cibuildwheel,
> import-tests all 13 compiled extensions in each, and publishes on a `v*` tag
> via Trusted Publishing. `biopaithon 1.88.dev0` shipped 35 wheels plus an
> sdist. Job runs per push: 29 → 19.

`ci.yml:114-136` runs a 3-OS × 5-Python matrix whose only step is
`python -m build --wheel`. The upload step is commented out (`:138-149`), so
every wheel is discarded — while `cleanup_wheels` (`:396-403`) deletes an
artifact that is never created. There is no `cibuildwheel`, `auditwheel`,
`delocate` or `twine` anywhere, and no release workflow at all. Linux wheels
built this way are tagged `linux_x86_64` and rejected by PyPI regardless.

This matters more for a fork: `pyproject.toml` now says `biopaithon` while the
import package is still `Bio`, so installing both distributions would silently
clobber the same directory. That collision is now documented in `README.rst`.

**Fix:** replace with `pypa/cibuildwheel` including a `CIBW_TEST_COMMAND` that
actually imports the compiled extensions; upload artifacts; delete the dead
cleanup job; add a tag-triggered release workflow using Trusted Publishing.
**Effort M · Impact high**

### 3.3 Caching defeats the "test against latest dependencies" intent **[FIXED]**

> **Status: fixed** in PR #63. `ci.yml` now has a weekly `schedule:` trigger
> (Mondays 03:00 UTC) on which the cache steps are skipped entirely, so the
> weekly run installs everything fresh and a new dependency release can
> actually break CI within a week; a bumpable `CACHE_EPOCH` env var is part
> of every cache key for manual invalidation.

`ci.yml:76-89` caches the whole `${{ env.pythonLocation }}` keyed only on
`pyproject.toml` + `ci-dependencies.txt`, and gates installation on a cache miss.
Since those files rarely change, `--upgrade-strategy eager` (`:89,194`) almost
never runs, so a new numpy release cannot break CI until an unrelated commit
touches those files — the exact failure mode that motivated `numpy!=2.1.0`
(`ci-dependencies.txt:16`). Also `ci.yml:372` builds a cache key from
`matrix.python-version` in a job that **has no `strategy:` block**, so the
expression is empty.

**Fix:** add a weekly `schedule:` trigger running uncached; add a bumpable
`CACHE_EPOCH` to the key; fix the docs job key.
**Effort S · Impact high**

### 3.4 Two of five supported Pythons are tested; MySQL is started but unusable **[MOSTLY FIXED]**

> **Status: both named halves are fixed.** PR #4 extended `test_linux` and
> `linux_prep` to 3.10–3.14, so every advertised version is now compiled *and*
> tested rather than compiled and discarded, and PR #87 made one Linux CI
> job actually run the BioSQL suite against the MySQL server that was
> previously started and never used. macOS and Windows still run only 3.10
> and 3.14.

Wheels are built for 3.10–3.14 (`ci.yml:119`) and `pyproject.toml:23-32`
advertises all five, but every test matrix is `["3.10","3.14"]`
(`:169,226,272`). Meanwhile `:173-175` runs `sudo /etc/init.d/mysql start`
before the Linux tests, while both drivers are commented out in
`ci-dependencies.txt:19-21` — so `Tests/common_BioSQL.py:93-100` skips, and the
BioSQL layer is effectively untested outside the fragile AppVeyor job.

**Fix:** extend the test matrix to all five, or drop the untested ones from the
wheel matrix so build and test agree. Either restore the MySQL driver and add a
`biosql.ini` step, or delete the pointless server start.
**Effort S · Impact high**

### 3.5 Three CI systems that have already drifted **[FIXED]**

> **Status: fixed** (commit `b793f8145`). `.appveyor.yml` and `.circleci/`
> are gone entirely — CircleCI was removed rather than kept for the docs
> deploy — and `requirements-sphinx.txt` now lives in `Doc/`, consumed by
> the GitHub Actions docs job.

GitHub Actions, CircleCI and AppVeyor all run the same offline suite with three
different dependency sets. `.appveyor.yml:4` pins the long-deprecated
`Visual Studio 2015` image and installs whatever Python Miniforge ships that
day, making failures unreproducible. `.circleci/config.yml:50` has a copy-paste
bug: `--source Bio,BioSQL --source Bio,BioSQL`. `.circleci/requirements-sphinx.txt`
is consumed by the *GitHub Actions* docs job (`ci.yml:381`).

**Fix:** delete AppVeyor (superseded by the `windows-2022` job), keep CircleCI
only for the docs deploy, move `requirements-sphinx.txt` to `Doc/`, fix the
duplicated `--source`.
**Effort M · Impact medium**

### 3.6 Tooling configuration is split and already inconsistent **[mypy half FIXED]**

> **Status: the mypy-cannot-see-numpy half is fixed** in PR #20. The hook in
> `.pre-commit-config.yaml` now has `additional_dependencies: [numpy]`, the
> `[mypy-numpy.*]` override is gone from `.mypy.ini`, and the five
> `Bio/PDB/Atom.py` errors are resolved (a `_isclose` helper handling
> `float | None`). **Still open:** `ci-dependencies.txt:10` still pins
> `black==22.12.0` against pre-commit's `24.10.0`, ruff configuration still
> lives only in pre-commit CLI args with no `[tool.ruff]` section,
> `ci-dependencies.txt` has not been replaced by `[dependency-groups]`, and
> there is no scheduled full-tree `pre-commit run --all-files`.

`ci-dependencies.txt:10` pins `black==22.12.0`; `.pre-commit-config.yaml:33` uses
`24.10.0`. All ruff configuration lives in pre-commit CLI args, with no
`[tool.ruff]` section anywhere — which is why running `ruff check` bare gives
different results from CI. `.pre-commit-config.yaml:112-120` skips flake8,
rstcheck, doc8 and codespell on pre-commit.ci, and GHA only lints changed files,
so those four linters never see the whole tree.

Critically for this fork: the mypy hook (`.pre-commit-config.yaml:44-49`) has
**no `additional_dependencies`**, so mypy runs without numpy and `.mypy.ini:66-67`
turns every numpy symbol into `Any`. Measured: without numpy, `Success: no issues
found in 298 source files`; with numpy 2.5.1, **5 real errors**, all in
`Bio/PDB/Atom.py:290-292` where `bfactor`/`occupancy` are `float | None` and are
passed to `np.isclose`. CI reports a clean run that is not clean.

**Fix:** add `additional_dependencies: [numpy]` and delete `[mypy-numpy.*]`; fix
the five `Atom.py` errors; move ruff config into `pyproject.toml`; drop the stale
black pin; replace `ci-dependencies.txt` with `[dependency-groups]`; add a weekly
full-tree `pre-commit run --all-files`.
**Effort M · Impact high**

### 3.7 No optional-dependency extras; unbounded `numpy` **[FIXED]**

> **Status: fixed** in PR #64. `pip install biopaithon[graphics]`,
> `[biosql]`, `[structure]`, `[phylo]`, `[all]` and `[test]` now work;
> `numpy` is bounded to `>=1.24,<3`; and `ci-dependencies.txt` just
> installs the `test` extra, so the two cannot drift. `[structure]`
> deliberately omits `mmtf-python` (§2.2).

`pyproject.toml` declares exactly `dependencies = ["numpy"]`, unbounded, and no
`[project.optional-dependencies]` — despite `.mypy.ini:43-77` and
`ci-dependencies.txt:17-27` documenting ten optional packages. Users have no
supported way to ask for the graphics, BioSQL or structure stacks.

**Fix:** add `graphics`, `biosql`, `structure`, `phylo` and `all` extras; have CI
install `.[all]` so extras are exercised; set `numpy>=1.24,<3`. Drop `wheel` from
`build-system.requires` (unnecessary since setuptools 70.1).
**Effort S · Impact medium**

### 3.9 The declared setuptools floor is too low for the license syntax in use **[reproduced — FIXED]**

> **Status: fixed in this fork.** `build-system.requires` is now
> `["setuptools>=77"]` and `wheel` has been dropped. The CI job that would test
> the floor is still outstanding.

`pyproject.toml:2` declared `requires = ["setuptools>=74.1", "wheel"]`, but
`license = "LicenseRef-Biopython-License-Agreement"` and the `license-files` list
are PEP 639 syntax, which setuptools only accepts from **77.0.0**. Any build that
resolves setuptools to 74.1–76.x fails while *reading* the config:

```
ValueError: invalid pyproject.toml config: `project.license`.
configuration error: `project.license` must be valid exactly by one definition
```

Reproduced against a system setuptools older than 77. Build isolation normally
hides this because pip fetches the newest setuptools, so it surfaces only for
`--no-build-isolation` builds, distro-packaged setuptools, and pinned or offline
environments — where it looks like a corrupt `pyproject.toml` rather than a
version floor.

**Fix:** raise the floor to `setuptools>=77` (which also covers the `>=74.1`
needed for the `ext-modules` table) and drop `wheel`. Add a CI job that builds
with `--no-build-isolation` against exactly the declared minimum, so the floor
is tested rather than assumed.
**Effort S · Impact medium**

### 3.8 The sdist ships ~108 MB of test data **[FIXED]**

> **Status: fixed** in PR #84, taking the "stay self-testing" branch of the
> decision below: `Tests/` and its data still ship in full, but unpacked
> size dropped 32% (124.1 MB to 84.3 MB). The five orphaned fixture
> directories (§4.8) are deleted, the 40 MB `Tests/PDB/6WG6.xml` ships
> gzipped and its test reads it through `gzip.open`, and CI fails if the
> compressed or unpacked sdist grows past 40 MB / 170 MB respectively —
> about twice today's sizes — so the payload cannot silently regrow.

`MANIFEST.in:15` is `recursive-include Tests *` with no filter. Tracked `Tests/`
totals 108 MB across 1,607 files; `Tests/PDB/6WG6.xml` alone is 39 MB — roughly
triple the entire `Bio/` package.

**Fix:** decide explicitly whether the sdist must be self-testing. If yes, gzip
the multi-MB fixtures (the suite already reads `.gz` fixtures) and cap the rest;
if no, use targeted includes. Either way add a CI assertion on sdist size so it
cannot silently regrow.
**Effort M · Impact medium**

---

## Tier 4 — Test suite and testing infrastructure

The `unittest.TestCase` tests themselves are fine and should stay — rewriting
them would be pure churn. The problems are in the bespoke *runner* and in how
the suite depends on ambient process state.

### 4.1 Five test modules fail if `PYTHONWARNINGS` is set **[reproduced by review — FIXED]**

> **Status: fixed** in PR #10. Every recorded-warning block now establishes
> its own filter — `assertWarnsRegex` where one warning is the contract, an
> explicit `simplefilter("always")` where a recorded list is inspected — and
> the five modules above pass under `PYTHONWARNINGS=ignore`.

39 sites use `warnings.catch_warnings(record=True)`; seven never call
`simplefilter("always")` inside the block, so they inherit whatever global filter
is active — `catch_warnings` saves and restores filters but does not reset them.

```
$ cd Tests && PYTHONWARNINGS=ignore python run_tests.py --offline \
    test_Entrez test_Align_msf test_PDB_PDBParser test_AlignIO test_SeqIO_features
FAILED (failures = 5)
```

Every failure is `AssertionError: 0 != 1` (`Tests/test_Entrez.py:144`,
`test_Align_msf.py:221`, `test_PDB_PDBParser.py:627`), which gives no hint of the
real cause. `PYTHONWARNINGS` is set by many CI images, editors and wrapper
scripts. The guard is unreliable in both directions: under an `ignore` filter a
genuinely regressed warning would still pass.

**Fix:** replace the `catch_warnings(record=True)` + `assertEqual(len(w), 1)`
pattern with `assertWarns`/`assertWarnsRegex`, which resets filters correctly and
is already used 184 times elsewhere in the suite. Failing that, add
`simplefilter("always")` as the first statement in every block.
**Effort S · Impact high**

### 4.2 The suite only runs from inside `Tests/` **[FIXED]**

> **Status: done.** Four sweep PRs moved every remaining module onto
> `support.DATA`: PR #162 (Ace, Affy and `Align*`), PR #164 (B to O),
> PR #163 (`[Pp]*`) and PR #166 (Q to Z), including the
> `os.chdir("Tests")` in `test_PDB_StructureAlignment.py`. Then PR #200
> dropped the session `chdir`: tests run in whatever directory pytest
> started in, and only the docstring examples change into `Tests/`. CI
> now starts the suite outside the checkout, so a new cwd-relative path
> fails there instead of hiding. The parallel runs this unlocked are
> §4.3. The note below is from before PR #162.
>
> **Status: the mechanism exists and five modules are migrated.** PR #93
> added `Tests/support.py` with `DATA = pathlib.Path(__file__).parent` —
> including guidance for tests that deliberately exercise relative paths —
> and freed the first five modules from the working directory. 108 of 217
> test modules remain cwd-relative; the migration is mechanical from here
> and tracked in `TODO.md`.

Data paths are relative, so `python Tests/test_GenBank.py` from the repo root
gives 83 `FileNotFoundError`s. The only thing making the suite work is
`os.chdir(self.testdir)` at `Tests/run_tests.py:237`. Individual modules then
fight over cwd — `test_SeqIO_index.py` calls `os.chdir` at twelve separate
lines, while `test_PDB_StructureAlignment.py:225` does `os.chdir("Tests")`, which
only works from the repo root, the exact opposite convention.

Every IDE runner, coverage tool and coding agent defaults to invoking from the
repo root and sees what looks like a broken checkout. It also forbids parallel
execution, since cwd is per-process global state. (Notably, `pytest` passes the
suite unmodified once cwd is right: `cd Tests && python -m pytest test_GenBank.py
test_Seq_objs.py -q` → `162 passed`.)

**Fix:** introduce `Tests/support.py` with `DATA = pathlib.Path(__file__).parent`
and migrate `open("GenBank/x.gb")` → `open(DATA / "GenBank/x.gb")` module by
module, removing the cwd dependency entirely and unlocking parallelism. Delete
the `os.chdir("Tests")` in `test_PDB_StructureAlignment.py:225`.
**Effort M · Impact high**

### 4.3 The runner reports 501 "tests" that are really 3,050 test cases, and discards their output **[FIXED]**

> **Status: pytest is the runner.** PR #175 replaced the bespoke runner
> with pytest 9, keeping the `unittest.TestCase` tests. Results are per
> test, a test's output is shown when it fails, and `-k`, `-x`, `--lf`,
> `--durations` and `--junitxml` all work. `run_tests.py` is
> now a shim that translates its old command line, so documented
> commands stay valid. The doctests are not collected with
> `--doctest-modules`, as the fix below proposed, but by name through
> `Tests/test_docstrings.py`, one test per docstring, with failures
> shown inline. From `Tests/`, `--doctest-modules` cannot see
> `Tests/conftest.py`, and a root conftest would import the source
> tree ahead of the installed package CI tests. PR #193 documented
> pytest in the Tutorial's testing chapter, `README.rst` and
> `CONTRIBUTING.rst`. Still open: parallel runs with pytest-xdist, CI
> calling pytest directly, and dropping the session `chdir` into
> `Tests/`.
>
> **Status, later: all three are done.** PR #194 runs the suite in
> parallel with pytest-xdist (`-n auto --dist loadfile`), and CI calls
> pytest directly, with coverage from pytest-cov and `--check-skips` on
> Linux. PR #200 dropped the session `chdir` (§4.2). One CI gap is left:
> the Windows test step is not a gate yet (`continue-on-error`), because
> 22 to 33 tests fail there on `main`; see "Follow-ups from the October
> 2026 sprint" below.

`Ran 501 tests` counts *modules* (209 unit-test + 292 doctest). Loading them
individually gives 3,050 actual test cases. Worse,
`output = StringIO()` (`run_tests.py:231`) is installed as `sys.stdout` at `:240`
and **never read anywhere** — grep returns only those two lines. Anything a test
prints is destroyed. This bites hardest in `test_Tutorial.py:255-266`, where a
`doctest.DocTestRunner` writes its got/expected diff to `sys.stdout` and the test
then raises only the failing doctest *names*.

The runner also offers no `-k` filtering, no per-test timing, no JUnit XML, no
parallelism, no skip accounting and no `--pdb`, in 250 lines of hand-rolled
result plumbing.

**Fix:** adopt pytest as the *runner* while keeping the `unittest.TestCase`
tests — add `[tool.pytest.ini_options]` to `pyproject.toml` and keep
`run_tests.py` as a thin shim so documented commands stay valid. Move the doctest
phase to `--doctest-modules` so failures print inline.
**Effort M · Impact high**

### 4.4 The previously documented test command did not work **[FIXED]**

> **Status: both halves are fixed.** PR #14 removed the module-scope
> `from setuptools import find_packages` from `Tests/run_tests.py`, so the
> runner works in a clean venv, and PR #64 added the `test` extra with
> `ci-dependencies.txt` reduced to installing it, so the two cannot drift.

`python setup.py test` does not work, yet the old `CLAUDE.md` documented it as
the primary way to run tests. `setup.py` is now a 19-line deprecation shim with
no `test` command. The exact failure depends on the setuptools version: with a
recent one it is `error: invalid command 'test'`; with one older than 77 it fails
earlier still, on the `license` field (see §3.9). **This fork has already fixed
the documentation side**: `AGENTS.md` documents `cd Tests && python
run_tests.py --offline`, which works.

Two related problems remain in the code: `run_tests.py:36` does `from setuptools
import find_packages` at module scope — needed only for doctest discovery — so
the runner fails with `ModuleNotFoundError: No module named 'setuptools'` in a
clean Python 3.12+ venv. And there is no `[project.optional-dependencies] test`;
the real test dependencies live in a bare `ci-dependencies.txt` that no
`pip install` or `uv sync` invocation references.

**Fix:** move the `find_packages` import inside the doctest branch or replace it
with `pkgutil.walk_packages`. Add a `test` extra mirroring
`ci-dependencies.txt`, and reduce that file to `-e .[test]` so the two cannot
drift. See also §3.7.
**Effort S · Impact high**

### 4.5 Skips are silent and unbounded, and a live API key is committed **[MOSTLY FIXED]**

> **Status: both main halves are fixed.** PR #15 removed the committed API
> key: `Tests/test_Entrez.py` uses the literal `"offline-test-api-key"`,
> and `Tests/test_Entrez_online.py` reads `NCBI_API_KEY` from the
> environment and skips when it is unset. PR #98 fixed the accounting: the
> summary now reads `Ran N modules (M cases), S skipped, F failed`, skip
> detection is a real `except MissingExternalDependencyError` around the
> import rather than substring matching (the same exception raised later,
> from code under test, is a failure again), and `Tests/expected_skips.txt`
> is a manifest the runner enforces — an unexpected skip fails the run.
> Still open: the macOS and Windows jobs install only numpy, so optional
> stacks remain untested on two of three platforms.

A fully-populated local run skips **61 of 501 modules (12%) and still exits 0**,
with no skip count in the summary. The macOS and Windows CI jobs install only
numpy, so on two of three platforms every reportlab, mmtf, msgpack, scipy,
networkx, igraph, rdflib and BioSQL module is skipped with no record. Skip
detection is substring matching on traceback text (`run_tests.py:250-256`), so a
`MissingPythonDependencyError` raised from *library code under test* also
converts to a green "skipping." line.

Separately: a live NCBI API key is committed at `Tests/test_Entrez.py:25` and
`Tests/test_Entrez_online.py:31`
(the literal value is in those files; it is deliberately not repeated here),
shared by every fork and CI run — a rate-limit and revocation hazard, and a
credential in version control regardless.

**Fix:** report `Ran N modules (M cases), S skipped, F failed`, and add an
expected-skip manifest so the Linux job fails if the count rises. Replace the
substring match with a real `except MissingPythonDependencyError` around the
import. Install the full dependency set on at least one macOS and one Windows
job. Move the API key to `os.environ.get("NCBI_API_KEY")` with a skip when unset.
**Effort S–M · Impact high**

### 4.6 Mega test methods with no `subTest` — one failure hides thousands of assertions

3,050 test cases carry roughly 110,000 assertions (95,611 `assertEqual`, 8,699
`assertAlmostEqual`), and the distribution is extreme:
`test_Align_tabular.py` is 43,592 lines with 18 test methods, one of which
(`test_2228_tblastx_001`) is 32,486 lines. `grep -rl subTest Tests/test_*.py`
returns **zero** — `subTest` is used nowhere.

A one-character parser regression fails at assertion #40 of a 30,000-assertion
method and the rest never execute, so you learn nothing about the blast radius.

**Fix:** convert the repetitive record-by-record loops in the worst offenders
(`test_Align_tabular`, `test_Align_hhr`, `test_pairwise_aligner`,
`test_Blast_parser`) to `with self.subTest(...)`. Move large inlined expected
output into fixture files so regenerating expectations is a data diff rather than
a 40,000-line source diff.
**Effort M–L · Impact medium-high**

### 4.7 Dead and unsound runner machinery **[FIXED]**

> **Status: the three remaining bullets are fixed** in PR #95, exactly as
> the fix line proposes: `--offline` now patches `socket.socket.connect`
> (allowing AF_UNIX and loopback), so `http.client`, `requests` and raw
> sockets are covered too; `requires_internet.check()` is lazy and
> time-boxed instead of a live `getaddrinfo` at import; and
> `test_PAML_tools.py` uses `shutil.which`.

- `run_tests.py:62-104` — a 43-line `if np is None:` block excluding 40 modules
  from doctests. numpy is a hard dependency, so this is unreachable (also §2.1).
  **[FIXED]**
- `run_tests.py:167-172` — `--offline` monkeypatches only
  `urllib.request.urlopen`. `http.client`, `socket` and `requests` are untouched,
  so the guarantee is not actually enforced.
- `Tests/requires_internet.py:24-31` does a live `socket.getaddrinfo` at *import*
  time of the 20 modules that use it, so a firewalled network stalls the run
  rather than skipping.
- `Tests/test_PAML_tools.py:26-47` — a hand-rolled `which()` whose entire
  Windows search path (`os_path`, lines 26-43) is built and then never used.
  `shutil.which` has existed since Python 3.3.

**Fix:** delete the `np is None` block; patch `socket.socket.connect` instead of
`urlopen`; make `requires_internet.check()` lazy and time-boxed; replace the
hand-rolled `which`.
**Effort S · Impact medium**

### 4.8 Test fixtures dominate the repository **[MOSTLY FIXED]**

> **Status: the actionable parts landed with §3.8.** PR #84 verified and
> deleted the five orphaned fixture directories and gzipped
> `Tests/PDB/6WG6.xml` (its test reads through `gzip.open`). The fixtures
> still dominate the repository by file count — that is the nature of a
> parser library's test suite and no further action is planned.

`Tests/` is 1,607 of 2,414 tracked files (67%) and 108 MB; the packed git repo is
77 MB. Beyond the sdist consequence in §3.8, roughly 1.1 MB is fully orphaned —
`Tests/FSSP/`, `Tests/CodonUsage/`, `Tests/MetaTool/`, `Tests/SubsMat/` and
`Tests/NeuralNetwork/` have no corresponding `Bio/` module and no test
references. The 40 MB `Tests/PDB/6WG6.xml` backs exactly two assertions, and its
`.bcif.gz` sibling is only 946 KB.

**Fix:** verify against `git log`, then delete the five orphaned directories and
gzip or truncate `6WG6.xml` (the PDBML parser can read through `gzip.open`).
**Effort M · Impact medium**

---

## How much is any of this used?

Several items above propose removing or deprecating a module, and until now the
argument for each rested on how dead the code *looks*. This section replaces
that with a measurement, so a removal decision can cite a number.

**Method.** GitHub code search, via `gh api -X GET search/code`, for the exact
import line with `language:python`. For example:

```bash
gh api -X GET search/code -f q='"from Bio import pairwise2" language:python' \
  -F per_page=1 --jq '.total_count'
```

Taken 2026-08-06.

**What the numbers are not.** They count *files* on public GitHub, not
projects, downloads or people, and notebooks and course material inflate the
popular entries. They are a current cross-section, not a trend: a real time
series would need PyPI BigQuery (package-level only, so no per-module
resolution), Software Heritage, or the Stack Overflow dumps. For deprecation
decisions the level matters more than the trend, so that gap is tolerable.

**One trap, found the hard way.** Searching for a bare name rather than an
import line counts Biopython's own tree. `"xbbtools"` returns 204, but
`"xbb_blastbg"` — a filename that exists nowhere else — returns 33, and
excluding `repo:biopython/biopython` moves 204 only to 202. Those are forks and
vendored copies. Only search for something a *user* would write.

### Baselines

| query | files |
|---|---|
| `from Bio import SeqIO` | 72,192 |
| `from Bio.PDB import` | 17,808 |
| `from Bio import Entrez` | 9,648 |
| `from Bio import Align` | 8,352 |
| `from Bio import AlignIO` | 6,400 |
| `from Bio import Phylo` | 4,488 |

### Modules with a pending decision

**"State" was this fork and upstream `master` on the same day**, when every
row was *identical* in both. That is no longer true: the fork has since acted
on most of these (the "item" column links the section recording what was
done), so the state column now notes where the two have diverged. The counts
remain the usage evidence the decisions rested on.

| module | files | state (here = upstream) | item | what the number says |
|---|---|---|---|---|
| `Bio.pairwise2` | **4,176** | **removed here** (PR #79, `ImportError` stub); deprecated upstream | §2.3 delete | Deprecated eight releases ago and still on the scale of `Bio.Phylo`. §2.3's "leave an `ImportError` stub" is not a courtesy, it is the whole change. |
| `Bio.Blast.NCBIXML` | 1,792 | **deprecated here** (PR #91); shipped, no warning, upstream | §1.5 deprecate | Real users on the superseded stack. Rewriting the tutorial first, as §1.5 sequences it, is right. |
| `Bio.Blast.NCBIWWW` | 1,040 | **deprecated here** (PR #91); shipped, no warning, upstream | §1.5 deprecate | As above. |
| `Bio.SearchIO` | 1,064 | shipped in both | §1.5 | Smaller than either old BLAST module, which bears on whether `blast-xml` is worth adapting. |
| `Bio.PDB.mmtf` | 154 | **deprecated here** (PR #72); shipped upstream | §2.2 remove | Low, and the format's server no longer resolves. Removal is defensible; a stub still costs nothing. |
| `Bio.Restriction` | 141 | shipped in both | §1.7 lazy imports | **Much lower than assumed.** Its 77 ms import cost is paid by few people, and no module in `Bio/` imports it, which is why the lazy-loading half of §1.7 was not done with the `CodonTable` half. |
| `Bio.HMM` | 35 | **`ImportError` stub here** (PR #92); empty package upstream | §2.1 | Imports of a package whose contents were removed in 1.86 — already broken for all 35. |
| `Bio.codonalign` | 27 | **warning dropped here** (PR #141); experimental warning upstream | §2.5 decide | Twelve years of an import-time warning, and this is the audience. Whatever is decided, it is not urgent. |
| `Scripts/xbbtools` | — | **removed here** (PR #86); shipped and broken upstream | §2.1 | Not measurable this way; see the trap above. Decide on the code, not a count. |

### The command line wrappers removed upstream in 1.86

The sharpest result, because these are not a prediction — the removal already
happened. **Upstream removed them in release 1.86 and this fork inherited that**;
neither tree ships them today. The fork's only difference is what an import
*says*, and only once PR #56 lands.

| module | files still importing it | here | upstream |
|---|---|---|---|
| `Bio.Align.Applications` | **1,784** | stub (PR #56) | absent |
| `Bio.Application` | **1,506** | stub (PR #56) | absent |
| `Bio.Blast.Applications` | **1,136** | stub (PR #56) | absent |
| `Bio.Emboss.Applications` | 284 | stub (PR #56) | absent |
| `Bio.Sequencing.Applications` | 214 | stub (PR #56) | absent |
| `Bio.Phylo.Applications` | 147 | stub (PR #56) | absent |

Roughly 5,000 files, against 27 for `Bio.codonalign`, which was never removed.
`Bio.Align.Applications` alone outweighs several modules still shipped. Every
one of those imports raises `ModuleNotFoundError`, which names a module and
says nothing about what replaced it.

That is the evidence for §2.3's position, applied retrospectively: silent
deletion is the expensive part, not the deletion.

## Where this fork differs from upstream

**Nothing has been removed here that upstream still ships, and nothing is
deprecated here that upstream has not deprecated.** `git diff --diff-filter=D
upstream/master main -- Bio/ BioSQL/` lists no files: this fork has deleted
none. Every module in the survey above is in the same state in both trees.

Nor did adopting upstream pull requests move anything: all seven adoptions in
`ADOPTED.md` are bug fixes, none is a removal, so there is no case where this
fork took a change ahead of upstream deciding it.

Measured against `upstream/master` on 2026-08-10, after Biopython 1.88 was
merged in, this fork is **102 commits ahead and 1 behind**, and **61 files
under `Bio/` and `BioSQL/` differ** — all by modification, none by deletion.
The single commit behind is upstream's post-release version bump, which this
fork deliberately does not take: it keeps `1.88.dev0`, because its 1.88 line
tracks the release's content plus its own changes rather than being that
release. Almost all of
that is bug fixes that upstream would presumably also want, and which
`UPSTREAM.md` tracks reporting back.

Two differences are deliberate policy rather than fixes, and are the ones a
user could notice:

- **Removed modules keep a signpost.** `Bio.Application`,
  `Bio.Align.Applications`, `Bio.Blast.Applications`,
  `Bio.Emboss.Applications`, `Bio.Phylo.Applications` and
  `Bio.Sequencing.Applications` ship as stubs raising `ImportError` with a
  migration message. Upstream ships nothing under those names. Nothing becomes
  usable — the stubs can only be read — and `except ImportError` behaves the
  same either way, `ModuleNotFoundError` being a subclass.
- **The contribution policy**, which is the reason the fork exists. See
  `AGENTS.md`.

Everything else in `DEPRECATED.rst` follows upstream. Where this fork adds an
entry of its own, it is marked in that file.

## Verified non-issues

Recorded so nobody spends time re-investigating:

- **`Bio.TogoWS`, `Bio.SCOP` and `Bio.ExPASy` are not dead.** Despite their age
  and plaintext `http://` URLs, `togows.dbcls.jp`,
  `scop.mrc-lmb.cam.ac.uk/legacy/` and the ExPASy CGI endpoints all returned
  HTTP 200 when probed. Do not remove them as dead-service cruft.
- **`Bio/PDB/PSEA.py:16`** points at `ftp://ftp.lmcp.jussieu.fr/`, which does
  refuse connections — but it is a docstring URL, and the parser still works
  against a locally installed `p-sea`. Fix the URL; keep the module.
- **`unittest.TestCase` is not the problem.** The tests themselves are fine, and
  pytest runs them today unmodified once the working directory is right. It is
  the bespoke *runner* that lacks filtering, parallelism and honest reporting.
  Do not rewrite the tests; replace the runner. See §4.3.
- **Fixture duplication is minor.** The bulk of `Tests/` is size, not
  redundancy — mostly small SCOP `.id` files and deliberate FASTQ round-trip
  pairs. Deduplication is not the win; §4.8 and §3.8 are.
- **`README.rst:239`** trips a `rstcheck` "enumerated list start value not
  ordinal-1" info message. This predates the fork and is intentional formatting.

---

## Harvested from upstream's stalled pull requests

142 pull requests are open upstream, some since 2021. Surveying them found a
seam: several are technically finished and maintainer-approved, blocked only on
provenance or on a reviewer who never came back. Contributors dual-licence
under BSD 3-Clause when they open a PR upstream, so these can be adopted here
with attribution — **check the box was actually ticked on each one before
taking it**.

Ranked with the same bias as the rest of this document: a silent wrong answer
outranks a crash, a crash outranks a missing feature, and anything reproduced
outranks anything inferred. Every item below was checked against this fork's
tree, not just read about.

> **Status: fourteen of the fifteen have been adopted** — #3897, #4866,
> #4450, #3812, #5127 (first commit only), #5181 (the `matrix.py` half only,
> as planned), #5175, and now #4390, #3911, #4938, #5121, #5244, #4918 and
> #4634. Rows are marked **adopted** below; `ADOPTED.md` records what was
> taken and what was left behind on each. The one remaining candidate is
> #5157 (PDBList mirror support).

| | upstream | defect | effort |
|---|---|---|---|
| 1 | [#3897](https://github.com/biopython/biopython/pull/3897) **adopted** | `DisorderedAtom.copy()` leaves `selected_child` pointing into the original, so the copy reads and transforms the original's coordinates | S |
| 2 | [#4866](https://github.com/biopython/biopython/pull/4866) **adopted** | `MeltingTemp.Tm_*` given a `SeqRecord` strips its repr down to ACGT and returns a plausible Tm for the garbage | S |
| 3 | [#4450](https://github.com/biopython/biopython/pull/4450) **adopted** | one non-standard month in a PDB `HEADER` aborts the whole parse with `ValueError: list.index(x)`, naming no file or field | S |
| 4 | [#3812](https://github.com/biopython/biopython/pull/3812) **adopted** | one residue with missing atoms destroys an entire HSExposure calculation, because `_get_cb` returns `(None, 0.0)` where the caller checks for `None` | S |
| 5 | [#5127](https://github.com/biopython/biopython/pull/5127) **adopted** | `polar_angle` is read from uninitialised memory wherever the radius is zero | S |
| 6 | [#4390](https://github.com/biopython/biopython/pull/4390) **adopted** | `auth_residues=False` staples auth insertion codes onto label numbering, giving residue ids that exist in neither scheme | S |
| 7 | [#5181](https://github.com/biopython/biopython/pull/5181) **adopted** | `calculate_consensus` raises `UnboundLocalError` on an all-zero column — take the `matrix.py` half only, see below | S |
| 8 | [#3911](https://github.com/biopython/biopython/pull/3911) **adopted** | the GenBank writer emits multi-line qualifiers that its own parser then rejects | M |
| 9 | [#4938](https://github.com/biopython/biopython/pull/4938) **adopted** | `PDBList` still fetches MMTF from a host RCSB decommissioned; BinaryCIF, which this fork already parses, cannot be fetched at all | S |
| 10 | [#5175](https://github.com/biopython/biopython/pull/5175) **adopted** | nine `assert`s validating input in production code, which vanish under `python -O` — the same defect as §0.8 | S |
| 11 | [#5121](https://github.com/biopython/biopython/pull/5121) **adopted** | mmCIF parsing is ~50% slower than it needs to be for want of a six-line fast path | S |
| 12 | [#5244](https://github.com/biopython/biopython/pull/5244) **adopted** | every SAM file written has `TLEN` 0 on every record | M |
| 13 | [#4918](https://github.com/biopython/biopython/pull/4918) **adopted** | GenBank `LOCUS` lines with molecule type `NA` are rejected outright | S |
| 14 | [#5157](https://github.com/biopython/biopython/pull/5157) | `PDBList` hardcodes wwPDB paths, so EBI and other mirrors cannot be used | S |
| 15 | [#4634](https://github.com/biopython/biopython/pull/4634) **adopted** | no way to resolve polytomies in `Bio.Phylo`, which many downstream tools require | S |

Two defects were found while surveying and belong to no pull request:

- `Bio/motifs/matrix.py:209` indexes `counts[3]` unconditionally, so any motif
  over an alphabet of fewer than four letters raises `IndexError`.
  **Fixed** in fork PR #47, which pads the ranked counts with zeros.
- `Bio/SeqIO/PdbIO.py:444-449` keys `cif-seqres` off label chain ids while
  `cif-atom` goes through `MMCIFParser(auth_chains=True)`, so the two parsers
  report different chain ids for the same file.

### Do not adopt

- **[#5016](https://github.com/biopython/biopython/pull/5016)** — 99,000 lines,
  almost all of it an accidental repo-wide `black` run; the author has said he
  is abandoning the branch. The algorithms are real; the pull request is not.
- **[#5085](https://github.com/biopython/biopython/pull/5085)** — 1,795 lines of
  C for melting temperature. The thread establishes that most of the speedup is
  reachable in ~40 lines of Python and numpy. Take the idea, leave the C.
- **[#4994](https://github.com/biopython/biopython/pull/4994)** — changes the
  module-level `translate()` gap default. §0.2 deliberately took the narrower
  fix; reversing that is a policy decision, not a bug fix, and applied alone it
  would not even work.
- **[#4170](https://github.com/biopython/biopython/pull/4170)** — supersedes the
  §0.1 fix with the paper's neighbour-correlated tables. More correct, but it
  would change every flexibility value a second time in one release, and no
  upstream reviewer could validate the science. If ever, then in the same
  release as §0.1 or not at all.
- **[#5181](https://github.com/biopython/biopython/pull/5181)'s second half** —
  the alphabet check it grew under review rejects `Motif("ACGT", ...)` on any
  ambiguity code. Routine input. Take the two-line `matrix.py` fix only.
- **[#4073](https://github.com/biopython/biopython/pull/4073)** — the GenBank
  writer fix is wanted, but it is bundled with a `Scanner.py` rewrite that
  changes line-joining for every GenBank file read. Take #3911's shape instead.

## Suggested sequencing

1. **§0.11–0.14 before anything else.** These are memory-safety defects in C
   code, three of them reproduced as segfaults, and §0.11 is reachable from a
   malformed input file. All four are S-effort. §0.11 in particular should be
   treated as a security fix, not a bug fix. *(Done — all four are
   fully fixed, §0.13's hardening included.)*
2. **The rest of Tier 0.** Mostly S-effort, and it is where users are actively
   getting wrong answers. §0.1–0.4 are reproduced and unambiguous. Fix the
   `Tests/test_ProtParam.py` fixture that currently pins §0.1 in place.
3. **Make the suite trustworthy before relying on it** — §4.1 (five modules fail
   under `PYTHONWARNINGS`), §4.4 (the runner breaks in a clean venv) and the
   committed NCBI API key in §4.5. All S-effort, and everything below is
   verified by running the tests, so this comes first among the non-urgent work.
   *(Done — all three specifics are fixed, and §4.4's `test` extra and
   §4.5's skip accounting have since landed too.)*
4. **§3.1 and §3.6** — supply-chain pinning and making mypy actually see numpy
   are both S-effort and make every later change safer to land. *(Done in the
   parts that motivated this step; each entry lists its remainder.)*
5. **§1.6 then §1.3** — turn on `check_untyped_defs` with a ratchet baseline
   before annotating, so the annotation work is checked as it lands.
   *(§1.6 is done. §1.3 is about half done: the harness and the sequence
   spine are in, and the I/O and alignment modules are open.)*
6. **§1.7 (import cost)** and **§1.11 (FASTQ hot path)** are the highest
   user-visible-value M-effort items, and §1.11(a) is nearly free.
   *(Both are done, §1.11(b) included. `Restriction` laziness stays
   deferred, by decision.)*
7. **Tier 2 cruft** can proceed in parallel with anything; start with §2.1,
   which breaks nothing because it is already broken. *(Done — §2.1, §2.3
   and §2.4 are removed; §2.2 is deprecated, and its removal waits until
   after the 1.88 release; §2.5 is decided and done; §2.6's `colour`
   aliases now warn, and go in the release after 1.88.)*
8. **§4.2 and §4.3** — removing the cwd dependency and moving to a real runner
   unlock parallelism and per-test reporting, which makes every later item
   cheaper to verify. *(Done — both are fixed.)*
9. **§1.1, §1.2 and §1.4** are the L-effort structural items. Do §1.4 first —
   it is the smallest, and it fixes §0.7 properly. *(All three are done.
   §1.1 was re-scoped: the shim was declined, and the five `Bio.Align`
   reader PRs it was replaced by have merged. §1.2 has its registry,
   `register_format` and entry-point plugins.)*
10. **§1.9 then §1.10** — releasing the GIL is what makes the free-threading
    work worthwhile, but §1.10's static-state cleanup is a prerequisite for
    doing §1.9 safely under a free-threaded build. *(The aligner and kdtrees
    kernels release the GIL, and free threading is supported; the
    multi-phase-init migration was rejected, see §1.10. PR #196 releases
    the GIL in the `cluster.c` kernels too; letting Ctrl-C interrupt them
    is the open remainder.)*

---

## Addendum — the 2026-08-16 post-queue sweep

After the 2026-08-14 queue merged (every numbered item above reflects that),
six parallel reviews swept the tree once more: security, correctness, code
quality, consistency, canonical idioms including CI/CD, and documentation.
Everything below is new — each reviewer read this document and `TODO.md`
first and reported only what they do not already cover. Two findings were
fixed the same day rather than listed: PR #91's deprecation of
`Bio.Blast.NCBIWWW`/`NCBIXML` had been merged into its stacked base branch
rather than `main` and was invisible to the tracking (re-landed in
[#109](https://github.com/dbolser/BioPAIthon/pull/109)), and the `__all__`
sweep had declared 45 merely-imported stdlib and NumPy names as public API
across 14 packages (fixed in
[#110](https://github.com/dbolser/BioPAIthon/pull/110)).

The security review otherwise came back clean: TwoBit/alignmentcounts/pwm
bounds, Entrez XML (XXE and path traversal), BLAST XML entity handling,
SnapGene/Xdna length fields, BioSQL parameterization and network defaults
were each checked and found guarded or already tracked. The one exception is
A.10 below.

### A. New correctness defects, all reproduced

> **Status: all twelve are fixed**, in PRs #115 to #127.

Ranked by user impact. A.1 and A.4 also reproduce on stock Biopython 1.85,
so they are upstream's too — queued in `UPSTREAM.md` pending a
duplicate-check against the tracker.

1. **TwoBit-backed sequences mis-slice three ways.** *(Fixed in PR #120.)*
   `Bio/SeqIO/_twoBitIO.c:401` computes `size = (end - start) / step` —
   truncating, not ceiling — so any extended slice with
   `(end - start) % step != 0` silently drops the final base:
   `record.seq[0::3]` loses the last codon position unless the length is a
   multiple of 3. `Bio/SeqIO/TwoBitIO.py:132-133` ignores negative steps
   when computing the byte range, so `record.seq[::-1]` raises
   `RuntimeError`. And `:124-131` has no upper bound on integer indices, so
   `record.seq[len(record)]` returns a base decoded from the *next* record's
   packed data instead of `IndexError`. **Effort S–M · Impact high.**
2. **`SeqRecord.upper()`/`lower()` share `letter_annotations` lists with the
   parent.** *(Fixed in PR #121.)* `Bio/SeqRecord.py:1181,1232` shallow-copy the dict, so
   mutating the derived record's `phred_quality` corrupts the original —
   the same disease §0.3 fixed for `annotations` at these sites; this
   dimension was missed. Only these two methods; the other derivations
   rebuild or drop the lists. **Effort S · Impact medium-high.**
3. **`SeqFeature._shift`/`_flip` share qualifier value lists** *(fixed in
   PR #125)*
   (`Bio/SeqFeature.py:280,297`) — so features on `record[2:8]` or
   `record.reverse_complement()` alias their qualifier lists with the
   parent's. §0.3's status note above claims these methods "do not share";
   that is true at the dict level and false at the value level. Fold into
   the planned features-copy work. **Effort S · Impact medium.**
4. **`BgzfWriter` crashes on poorly-compressible data.** *(Fixed in PR #126,
   which caps block data at `0xff00` bytes.)*
   `Bio/bgzf.py:834-837`: ~64 KB of high-entropy input deflates to more
   than 65,536 bytes and raises `RuntimeError: TODO - Didn't compress
   enough`; htslib caps per-block input instead. The guard also admits
   compressed sizes 65,511–65,536, which then fail in `struct.pack("<H",
   ...)` with a misleading error. **Effort S–M · Impact medium.**
5. **`Medline.read()` returns the first record of many** *(fixed in PR #115;
   PR #170 also made `parse()` treat whitespace-only lines as blank)* and raises bare
   `StopIteration` on an empty file (`Bio/Medline/__init__.py:218-219`),
   breaking the exactly-one contract every sibling `read()` enforces and
   `AGENTS.md` documents. **Effort S · Impact medium.**
6. **`Bio/motifs/__init__.py:208` is `raise Exception(ValueError, "...")`**
   *(fixed in PR #118)*
   — a typo; the intended `ValueError` is never raised and `except
   ValueError` misses it. **Effort S · Impact low-medium.**
7. **The GFA2 parser discards the mandatory length field** *(fixed in
   PR #119)* for `*`
   sequences (`Bio/SeqIO/GfaIO.py:199-207`): `S\ts1\t100\t*` parses to
   `len(record) == 0`, while the GFA1 branch preserves its optional
   `LN:i:` equivalent. **Effort S · Impact low-medium.**
8. **`GC123("")` raises `ZeroDivisionError`** *(fixed in PR #127)*
   (`Bio/SeqUtils/__init__.py:207`) — the per-position divisions are
   guarded, the total is not. **Effort S · Impact low.**
9. **`MafIndex.close()` closes the SQLite connection but not the MAF
   handle** *(fixed in PR #116; §1.1 later dropped the port)* (`Bio/AlignIO/MafIO.py:308-317`) — the Windows file-deletion
   case its own docstring exists for still fails. Should travel with the
   §1.1 MafIndex port. **Effort S · Impact low.**
10. **`Bio/Nexus/cnexus.c:34-35` sets `PyErr_NoMemory()` on allocation
    failure and falls through** *(fixed in PR #123)* to write through the NULL pointer — the
    §0.12/§0.13 error-path class, in the one extension that sweep missed.
    Reachable only under memory exhaustion. **Effort S · Impact low.**
11. **Two independent exception classes are both named `NexusError`**
    *(fixed in PR #122)*
    (`Bio/Nexus/Nexus.py:61`, `Bio/Nexus/StandardData.py:13`), and
    `Nexus.py:1057` calls into StandardData, so `except NexusError` around
    a parse can miss the parse error. **Effort S · Impact low-medium.**
12. **`Bio/SCOP/__init__.py:266-270` prints a missing sunid to stdout and
    then crashes on the same lookup** *(fixed for the HIE file in PR #117;
    the same lookup for the CLA file is a follow-up below)* with a bare `KeyError`. **Effort S ·
    Impact low.**

### B. Consistency debts (decide once, then mechanical)

> **Status: all decided, and done where the decision called for a change.**

- `Bio.Align.read`'s docstring documents its first parameter as `source`;
  the signature says `handle` (`Bio/Align/__init__.py:4937-4941`), so the
  documented keyword is a `TypeError`. Its sibling `parse` says `source`.
  One-word fix, whichever way §1.2's registry work settles the convention.
  *Fixed in PR #136: the docstring now says `handle`.*
- `Bio/PDB/StructureAlignment.py:64-71` deprecates with stdlib
  `DeprecationWarning` — invisible by default — instead of
  `BiopythonDeprecationWarning`, and has no `DEPRECATED.rst` entry.
  *Fixed in PR #134, with the entry.*
- `Bio/SeqIO/GfaIO.py:35-95` issues 13 parse-time malformed-input warnings
  as plain `BiopythonWarning` where ~55 sites elsewhere use
  `BiopythonParserWarning`. *Fixed in PR #129.*
- Parser exception bases are split 11 `ValueError`-based / 13 plain
  `Exception`-based, and ~30 sites `raise Exception(...)` inside files
  whose own convention is `ValueError` (worst: `Bio/Align/bigbed.py:795`,
  `a2m.py:101-109`). Widening bases to `ValueError` only increases
  catchability. *Fixed: PR #145 made the malformed-input exceptions
  subclass `ValueError`, and PR #135 replaced the bare
  `raise Exception(...)` sites.*
- First-parameter naming is `handle` (majority) vs `source` vs `file`, and
  the format argument is `format` (five packages) vs `fmt` (three);
  `Bio.motifs.write()` alone among the `write()`s takes no file and returns
  a string. Docs-first migration; keyword renames need deprecation shims.
  *Decided: keyword parameters are never renamed, neither `handle`/`source`
  nor `format`/`fmt`. The fork keeps upstream's names, so code moves
  between the two unchanged. Only `Bio.Align.read`'s docstring was wrong,
  and PR #136 fixed it.*
- The fork's own #88 indexing protocol (`parse_id_from_header`,
  `Bio/SeqIO/Interfaces.py:59` and five format modules) is public-named but
  docstring-marked "(PRIVATE)" — decide its status once, before downstream
  format authors copy it. *Decided: it is public. PR #133 documented it
  as a hook for format authors, and §1.2's `SeqIO.index` support for
  registered formats relies on it.*
- Docstring parameter style is three-way (reST `:param` / `Arguments:`
  lists / numpydoc), with `Bio/PDB/Atom.py` using two styles in one file.
  Fold into the §1.3 annotation pass rather than sweeping separately.
  *Decided that way: `AGENTS.md` now says to convert only in a module
  being annotated anyway.*

### C. Structural quality (M-effort, high leverage)

> **Status: all done.**

- **The dN/dS and McDonald–Kreitman suite exists twice** —
  `Bio/Align/analysis.py` and `Bio/codonalign/codonseq.py` — and the copy
  has already cost double maintenance: PR #25 fixed three bugs in one, PR
  #35 re-fixed the identical three in the other. Whatever §2.5 decides,
  make `cal_dn_ds`/`mktest` delegate to `Bio.Align.analysis` and delete
  ~1,000 duplicated lines. **Effort M.** *Fixed in PR #154 (about 1,230
  lines), which also fixed two bugs the copies did not share.*
- **Newick trees have two full parser stacks** (`Bio/Nexus/Trees.py` +
  `Nodes.py` vs `Bio/Phylo/NewickIO.py`), and `Phylo.read(f, "nexus")`
  parses tree text with the *old* one, then converts node-by-node through a
  recursive helper that dies on deep trees. Scoped fix: feed the tree
  strings Nexus extracts directly to `NewickIO.Parser`; the full merge is
  §1.1-scale and separate. **Effort M (scoped).** *Fixed, scoped, in
  PR #143: `Phylo.read(f, "nexus")` now parses trees with `NewickIO`, so
  deep trees no longer raise `RecursionError`.*
- `Bio.PopGen.GenePop.LargeFileParser` has zero importers and zero tests,
  and parses missing alleles differently from the live parser in the same
  package (literal `0` vs `None`); `FileParser.FileRecord` never closes its
  handle. Deprecate the dead parser per house style; give `FileRecord` a
  `close()`/context manager. **Effort S.** *Fixed: PR #139 deprecated
  `LargeFileParser`, and PR #140 gave `FileRecord` `close()` and context
  manager support.*
- `read_PIC` (`Bio/PDB/PICIO.py:41-815`) is a 774-line function of eleven
  nested closures — decompose into a private parser class; and
  `_get_atom_radius` (`Bio/PDB/ResidueDepth.py:108-491`) is a 383-line
  elif transcription of MSMS's radius table that should be a data table.
  `Bio/PDB/internal_coords.py:4940`'s `MissingAtomError` is advertised and
  never raised. **Effort M each.** *Fixed: PR #158 split `read_PIC` into a
  private parser class, PR #151 made the radius chain a rule table, and
  PR #144 deprecated `MissingAtomError`.*

### D. CI/CD and packaging canon (all S)

> **Status: done, apart from the `Scripts/` bullet and one checkout.**
> The `style` job's checkout (`ci.yml:46`) still lacks
> `persist-credentials: false`; it is listed under "Follow-ups from the
> October 2026 sprint".

- Replace `tj-actions/changed-files` + the 12-line heredoc in
  `ci.yml:47-90` with pre-commit's native `--from-ref/--to-ref` — removes
  the third-party dependency (the CVE-2025-30066 one, currently SHA-pinned)
  outright. *Fixed in PR #150.*
- No `timeout-minutes` outside the test jobs and no `concurrency:` group on
  `release.yml`; a hung cibuildwheel burns six hours, and two tags pushed
  close together interleave publishes. *Fixed in PR #161.*
- `persist-credentials: false` on the ten checkouts (nothing pushes);
  decide `fail_ci_if_error` for the three codecov uploads; delete the no-op
  `zip-safe` line; `packages.find` could replace the 68-entry list but
  needs a wheel-contents diff first (three data-only directories are
  load-bearing). *Mostly fixed: PR #153 dropped `zip-safe` and the
  persisted credentials on every checkout but the `style` job's, and
  PR #148 moved to `packages.find`. Decided:
  `fail_ci_if_error` stays `false` on the codecov uploads, so a codecov
  problem never fails CI.*
- `Scripts/` still uses `getopt`/`optparse` in four files — low value,
  legacy demos. *Not done.*

### E. Documentation corrections

> **Status: all done.** Some of the PRs found more stale text; that is
> listed under "Follow-ups from the October 2026 sprint" below.

- `README.rst:36-38` claims upstream's documentation "still describes this
  fork accurately" — no longer true post-sweep (pairwise2, GenBank length
  errors, mmtf, extras). `:120` claims Python 3.15-rc testing no CI job
  performs. The extras exist only in NEWS; the documented dev install
  (`--group dev`) under-installs relative to `.[test]` and silently skips
  the optional-stack tests. *Fixed in PR #137.*
- `Doc/Tutorial/chapter_pdb.rst:87-117,317-334` still teaches
  `Bio.PDB.mmtf` — including URL fetches against a host that has been
  DNS-dead since 2024 — with no deprecation note, directly above its own
  BinaryCIF section. `chapter_testing.rst:323-359` quotes a `run_tests.py`
  excerpt that no longer exists. *Fixed in PRs #131 and #138.*
- `Doc/Tutorial/chapter_contributing.rst:17-18` sends this fork's bug
  reports to upstream's issue tracker, contradicting `README.rst:30`;
  `CONTRIBUTING.rst:62-64,150-151` and `chapter_introduction.rst` similarly
  present upstream links as this project's. *Fixed in PR #130.*
- `Bio/SeqIO/__init__.py`'s format table is missing five readable formats
  (`embl-cds`, `genbank-cds`, `twobit`, `fasta-blast`, `fasta-pearson`) and
  names "clustalw" where the format is "clustal". *Fixed in PR #132.*
- `NEWS.rst`'s in-progress section announces pairwise2's removal and then
  describes pairwise2 crash fixes as current behaviour further down — fold
  the fixes into the removal entry. *Fixed in PR #128, by dropping the
  crash-fix entries rather than folding them in: they fixed a module that
  is no longer shipped.*

## Follow-ups from the October 2026 sprint

The sprint merged 91 pull requests, every number from #112 to #202
(#124 is the upstream sync). Their reviewers logged 139 follow-up
defects along the way. Each was checked on
2026-10-13 against `main` at `4cdff0207`, by reading the code and, where
the report gave a one-line reproduction, by running it (Python 3.12).

- **Still real: 112.** Listed below, grouped by area, with related reports
  merged into one entry.
- **Unconfirmed: 1.** The Infernal `*[ 0]*` hypothesis from PR #180 needs
  a real Infernal file to settle. It is listed below and marked so.
- **Already fixed: 16**, in the PR that found them or a later one. For
  example the Infernal 12-column crash, the `FormatRegistry` race, and
  the README's `Bio.DocSQL` example.
- **Duplicates: 6.** The same defect logged from two PRs.
- **Not defects in this tree: 4.** Two process notes, a CPython bug in
  `tempfile.SpooledTemporaryFile.rollover()` (for CPython's tracker, not
  this one), and a note for this bookkeeping, done here (addendum E).

"Upstream too" means the same code is on upstream `master` (`fc0c15678`,
checked with `git grep`, not run there). Those entries are also queued in
`UPSTREAM.md`.

### Decisions recorded at the close

- **Keyword parameters are never renamed.** Not `handle`/`source`, not
  `format`/`fmt`. The fork keeps upstream's names. Only `Bio.Align.read`'s
  docstring was wrong, and PR #136 fixed it (addendum B).
- **`parse_id_from_header` is public**, a hook for format authors (PR #133,
  addendum B).
- **`Bio.codonalign`'s experimental warning is dropped** and its API
  committed (PR #141, §2.5).
- **Codecov uploads keep `fail_ci_if_error: false`** (addendum D).
- **§1.5(c), a `blast-xml` adapter over `Bio.Blast`, was not done.** It is
  recorded only.
- **§2.2's second step, removing `Bio.PDB.mmtf`, waits until after the
  BioPAIthon 1.88 release.**
- **§1.7's lazy `Bio.Restriction` stays deferred.**
- **`cif-seqres` chain ids**, as decided in PR #169: an `auth_chains`
  argument names chains by author id (`True`) or label id (`False`). The
  default stays label ids, and warns for files where the two differ,
  because a later release will switch it to author ids. Cross-references
  are matched by author id in every mode.
- **Format plugin entry-point groups are `biopaithon.*`** (§1.2).
- **§1.1's `AlignIO` shim was dropped** in favour of closing `Bio.Align`'s
  gaps (PR #142).

### Planned steps not written

These steps were planned for the sprint and never started:

- §1.12: `pairwise-linear-gotoh-global`, `pairwise-linear-sw` and
  `pairwise-linear-gotoh-local`.
- §1.10 (and §1.9's Ctrl-C half): `cluster-interrupt`.
- §1.3: `annotate-seqio-alignio`, `annotate-align-msa-io`,
  `annotate-alignment` and `annotate-aligners`.

Every other planned step merged, including the seven last checked:
`ft-release` (PR #201), `drop-session-chdir` (PR #200),
`pairwise-linear-nw-on` (PR #199), `cluster-release-gil` (PR #196),
`annotate-seqrecord` (PR #187), `align-phylip-relaxed` (PR #195) and
`register-format-plugins` (PR #198).

### Crashes in C extensions

All are reached through a caller's arguments or data, not by parsing a
file. The first two and the two in `_pairwisealigner.c` were reproduced
on 2026-10-13 as a segfault or an abort; the other two were read.

- `Bio/Cluster/cluster.c` (`clusterdistance`, method `"m"`): an index list
  longer than the matrix has rows overflows the heap.
  `clusterdistance(data, index1=[0] * 5000, index2=[1, 2], method="m")` on a
  3-row matrix segfaults. Upstream too. (PR #196)
- `Bio/Cluster/cluster.c` (`svd`, through `pca`): data holding `inf`, or
  with a column mean that overflows (`1e308`), makes it read `w[-1]`; both
  segfault. Upstream too. (PR #196)
- `Bio/Cluster/clustermodule.c` (`py_pca`): the private `_cluster.pca`
  never checks the eigenvalues array's length, so a short one is written
  past its end. Its `pc` and `coordinates` errors print `mean`'s shape,
  and read past the end of a 1-D shape. Upstream too. (PR #196)
- `Bio/Align/_pairwisealigner.c`: `del aligner.mode` and
  `del aligner.wildcard` segfault, because no setter checks for `NULL`.
  One check in `Aligner_set_locked` covers them all. Upstream too.
  (PR #190)
- `Bio/Align/_pairwisealigner.c` (`Aligner_str`): a substitution matrix
  whose class name is over 1 KB overflows a 1,024-byte stack buffer, so
  `str(aligner)` aborts with "buffer overflow detected". The name is then
  also used as a format string. Upstream too. (PR #173)
- `Bio/Cluster/clustermodule.c`: `PyTree_length` returns `int` but is
  installed as a `lenfunc`. Undefined behaviour in principle, harmless on
  x86-64 and arm64. Upstream too. (PR #189)

### Sequences and SeqIO

- `Bio/SeqIO/FastaIO.py:390`: `FastaBlastIterator` has the typo
  `self_line = line`, so every record gets the first record's title:
  `>a`, `>b`, `>c` read as `a`, `a`, `a`. Upstream too. (PR #132)
- `Bio/Seq.py` (`join`): the type check uses up a generator, so
  `Seq("N").join(x for x in ["A", "C"])` is `Seq('')`. Upstream too.
  (PR #176)
- The `parse_id_from_header` hooks (`FastaIO`, `QualityIO`, `PirIO`,
  `PhdIO`, `AceIO`) split raw bytes where parsing splits the decoded text.
  An id with a control or Unicode space in it (`>a\x1cb desc`) gets an
  index key that differs from `record.id`, so the lookup raises "Key did
  not match". Upstream too. (PR #133)
- The same hooks decode as UTF-8 while `SeqIO.parse` opens files in the
  locale's encoding, so non-ASCII ids could differ between the two on a
  non-UTF-8 locale. Inferred from the code, not run. (PR #133)
- `Bio/SeqIO/_index.py` (`SequentialSeqFileRandomAccess.__init__`):
  replacing a marker-proxied built-in iterator with a plain function makes
  `SeqIO.index` raise `AttributeError` (no `record_start_marker`) instead
  of a readable `ValueError`. (PR #185)
- `Bio/SeqIO/_index.py:464,544,557`: the UniProt XML (`:464`) and
  IntelliGenetics (`:544`, `:557`) indexers compare summed line lengths
  with `handle.tell()`, which is a virtual offset in multi-block BGZF, so
  valid input fails an `assert` (see §0.8). Upstream too. (PR #180)
- `Bio/SeqIO/SffIO.py`: the header reads of the flow characters, key
  sequence and padding are not length-checked, so a truncated header gives
  a misleading warning and then `struct.error`. A negative name padding is
  not checked either, so `read(-n)` takes the rest of the file. Upstream
  too. (PR #159)
- `Bio/SeqIO/SffIO.py` (`SffWriter.write_records`): nothing stops a second
  call, which writes a second header mid-file. Upstream too. (PR #159)
- `Bio/SeqIO/InsdcIO.py`: a locus name of 16 characters or fewer with a
  sequence of 10^11 bp or more raises `ValueError` instead of using the
  long LOCUS line. Upstream too. (PR #159)
- `Bio/SeqIO/GfaIO.py:70`: `parts[1] not in "AifZJHB"` is a substring
  test, so an empty or multi-letter tag type (`AB::x`, `AB:if:x`) passes
  without a warning. Upstream too. (PR #129)
- `Bio/SeqIO/QualityIO.py` (`_get_phred_quality`): annotated
  `list[float] | list[int]`, but returns whatever is stored, which can now
  be an `array.array`. (PR #167)
- `Bio/SeqIO/PdbIO.py` (`CifSeqresIterator`): the docstring says the
  sequences come from `_entity_poly_seq`; the code reads
  `_pdbx_poly_seq_scheme`. Upstream too. (PR #169)
- `Bio/SeqIO/PdbIO.py`: `cif-seqres` labels DNA and RNA chains `protein`
  and turns nucleotides into `X` (1LCD). Known upstream as #3125, and a
  `TODO` in the code. (PR #169)
- `Bio/SeqIO/PdbIO.py`: `cif-atom` cannot name chains by label id
  (upstream's draft #4928 adds `cif-atom-label`). Decide whether SeqIO
  needs a format name for that before the `cif-seqres` default switches.
  (PR #169)
- `Bio/SeqFeature.py`, four small ones: `ExactPosition(5) + 1.5` is
  `ExactPosition(6)` where the annotation says `float`; the
  double-complement `LocationParserError` lacks its `f` prefix and prints
  `{text}`; `Reference() == None` raises `AttributeError`; and
  `CompoundLocation(generator)` raises `TypeError`. Upstream too.
  (PR #188)
- `Bio/SwissProt/__init__.py:783`: features get `id=None` where every
  other parser uses `"<unknown id>"`. Decide when `Bio.SwissProt` leaves
  the mypy baseline. (PR #188)

### Alignments

- `Bio/Align/emboss.py:33,44,89,108`:
  `raise ValueError("Unexpected line: %s") % line` raises `TypeError`.
  Upstream too. (PRs #155 and #182)
- `Bio/Align/exonerate.py:209`: an error message never gets its `%`
  arguments. Upstream too. (PR #155)
- `Bio/Align/sam.py:97`: the writer writes topology as `PP:` instead of
  `TP:`, so it does not round-trip. Upstream too. (PR #155)
- `Bio/Align/analysis.py` (`calculate_dn_ds`): NG86, LWL85 and ML do not
  check that the two rows have the same number of codons; only YN00 does.
  Read, not run. Upstream too. (PR #155)
- `Bio/Align/analysis.py` (`_yn00`): identical sequences raise
  `RuntimeError: dictionary changed size during iteration`
  (`cal_dn_ds(s, s, method="YN00")`). Upstream too. (PR #154)
- `Bio/Align/analysis.py` (`_count_diff_YN00`): the two-difference branch
  counts the step from `codon1` to `codon` twice and never the step to
  `codon2`; the three-difference branch passes `k[1]` where `k[2]` looks
  intended. Check against PAML's `yn00.c` before changing. Upstream too.
  (PR #154)
- Two upstream defects that PR #154 fixed here: `codonalign.mktest`'s
  `"-" in all_codon` never matches, so a partly gapped codon raises
  `KeyError: '---'`; and `_get_pi`'s F1x4 counts with a `Counter`, so a
  base absent from both sequences raises `KeyError` (upstream, masked by
  `UPSTREAM.md` row 4). Queued in `UPSTREAM.md` only. (PR #154)
- `Bio/Align/bigbed.py:493`: the writer's zoom record count disagrees for
  a zero-length item at the end of a chromosome (now a `RuntimeError`);
  needs a real writer fix. Upstream too (an `assert` there). (PR #155)
- `Bio/Align/stockholm.py`: the writer raises `AttributeError` for an
  `Alignment` without `column_annotations`. Upstream too. (PRs #155 and
  #197)
- `Bio/Align/stockholm.py`: in the reader, a `#=GR` line with no value
  raises `IndexError`, and the width error ends with a stray `)`.
  Upstream too. (PR #197)
- `Bio/Align/stockholm.py`: the reader silently drops a record that lacks
  its `//` when the next record's header follows; `AlignIO` reads both.
  (PR #197)
- `Bio/Align/phylip.py`: the strict writer does not reject names that
  truncation to 10 characters makes equal, so it can write an ambiguous
  file. `AlignIO` raises. Upstream too. (PR #195)
- `Bio/Align/phylip.py`: a later header that understates the length
  reuses the first header's message, so an interleaved file reports a
  sequence line as a bad header: `"2 4\ns1        AAAA\ns2        CCCC\n\nAAAA\nCCCC\n"`
  raises "Expected two integers in the first line, received 'AAAA'".
  (PR #184)
- `Bio/Align/phylip.py`: in the sequential layout, a blank line after the
  header counts as one of the N records, so `"2 4\n\ns1        AAAA\ns2        CCCC\n"`
  raises "Expected all sequences to have length 4; found 0" instead of
  saying the record count is wrong. Only the message is wrong; nothing
  is misread. (PR #184)
- `Bio/AlignIO/PhylipIO.py`: `RelaxedPhylipWriter` rejects ASCII
  whitespace only, but the reader splits on any whitespace, so it can
  write a file it cannot read back. `Bio.Align`'s new writer checks
  `str.isspace()`. Upstream too. (PR #195)
- `Bio/Align/__init__.py`: `PairwiseAligner(scoring="blastp").score(list(
  "HEAGAWGHE"), list("PAWHEAE"))` raises "sequence contains letters not
  in the alphabet", because list input is mapped twice. Upstream too.
  (PR #181)
- `Bio/Align/_pairwisealigner_linear.h`: on x87 builds
  (`FLT_EVAL_METHOD != 0`), excess precision could make the linear and
  full paths differ. Consider not taking the linear path there. (PR #186)
- `Bio/Align/tabular.py:328`: BLASTN or BLASTP output with a `query seq`
  column, common in `-outfmt 7`, raises "Unknown program". Upstream too.
  (PR #135)
- `Bio/AlignIO/ClustalIO.py`: a file holding only the header line hangs
  forever. Upstream too. (PR #160)
- `Bio/AlignIO/MauveIO.py`: line 259 passes `line` as a second argument to
  `ValueError`, so the message is never formatted; and `_ids` is a class
  attribute, so ids leak between iterators. Upstream too. (PR #160)
- `Bio/AlignIO/FastaIO.py`: a bare `>>` line raises `IndexError`.
  Upstream too. (PR #160)

### BLAST and SearchIO

- `Bio/Blast/__init__.py:1277`: `qblast(format_type="HTML")` can never
  succeed, because the polling loop's `<!DOCTYPE html` branch also
  matches the result page. Upstream too. (PR #157)
- `Bio/Blast/_parser.py` (`_start_blastxml2`): a root element with no
  namespace fails with "not enough values to unpack", and a missing
  `schemaLocation` with `KeyError`. Upstream too. (PR #157)
- `Bio/Blast/_parser.py` (`XMLHandler._externalEntityRefHandler`): it
  clears itself only on the call that first loads the DTD, so how an
  external general entity is handled depends on whether a v1 file was
  parsed earlier in the process. (PR #157)
- `Bio/SearchIO/BlastIO/blast_xml.py:678` (`BlastXmlIndexer.get_raw`):
  loops forever if the file was truncated after its `index_db` was built.
  Upstream too. (PR #156)
- `Bio/SearchIO/InfernalIO/infernal_text.py:489`: `assert gap_len > 0`
  probably rejects valid `*[ 0]*` output. Unconfirmed: it needs a real
  Infernal file (see §0.8). (PR #180)
- `Bio/SearchIO/__init__.py` (`index_db`): a lone `pathlib.Path` is not
  wrapped in a list, as a `str` is. Upstream too. (PR #164)

### Structures

- `Bio/PDB/binary_cif.py` (`_string_array_decoder`): ignores the column
  mask. On current RCSB BinaryCIF it raises `IndexError`, or returns wrong
  values without error (`pdbx_PDB_ins_code` on 1GBT). The MMTF
  deprecation sends users to this parser. Upstream too. (PR #131)
- `Bio/PDB/binary_cif.py` (`source.endswith`) and `FragmentMapper`
  (`fdir + "/"`) do not accept `os.PathLike`, unlike the other parsers.
  Upstream too. (PR #163)
- `FastMMCIFParser`, `PDBMLParser` and the BinaryCIF parser do not
  upper-case `type_symbol`, so an element such as `Fe` raises in `Atom`.
  `MMCIFParser` does. Upstream too. (PR #165)
- `Bio/PDB/StructureAlignment.py:150`: a custom `aligner=` raises
  `UnboundLocalError`, because `aligner` is set only in the `else`
  branch. The Tutorial shows this call. Upstream too. (PRs #134 and #165)
- `Bio/PDB/StructureAlignment.py`: there is no warning-free way to pass
  your own alignment, though the Tutorial offers `fasta_align=` as the
  alternative. Decide what is deprecated, then name the replacement in
  the message. (PR #134)
- `Bio/PDB/PICIO.py` (`_default_dihedron`): two unreachable branches
  would crash if reached (`list.add`). Upstream too. (PR #158)
- `Bio/PDB/PICIO.py`: `_process_dihedron` is annotated `-> set` but
  returns a tuple or `None`, and `read_PIC` is annotated `-> Structure`
  but can return `None`. `_ake_recurse` also repeats its recursion once
  per element. (PR #158)
- `Bio/PDB/internal_coords.py:4617`: the `AtomKey` docstring lists
  `(52, None, 'G', 'CA', ...)` as accepted; it raises. Upstream too.
  (PR #135)
- Upstream still has `MissingAtomError`, which nothing raises and its
  docstring advertises; deprecated here in PR #144. Queued in
  `UPSTREAM.md` only. (PR #144)
- `Bio/PDB/kdtrees.c`: `repr()` of `Point` and `Neighbor` raises
  `SystemError` (`PyUnicode_FromFormat` has no `%g`), and
  `KDTree(coords, bucket_size=10)` ignores the keyword. Upstream too.
  (PR #191)

### Phylogenetics

- `Bio/Phylo/NewickIO.py`: only one space is allowed after `:`, so
  `((a:  0.1,b:\t0.2),c);` gives leaves named `0.1` and `0.2` with no
  branch lengths, losing `a` and `b`. NEXUS input is unaffected. Upstream
  too. (PR #143)
- `Bio/Phylo/CDAOIO.py`: `Phylo.read(StringIO(...), "cdao")` raises
  `AttributeError` on `handle.name`, and the "RDFlib v3.2.1 or later"
  error is swallowed by the module's own `except ImportError`. Upstream
  too. (PR #174)

### Clustering

- `Bio/Cluster/__init__.py` (`somcluster`): a plain list raises
  `AttributeError`, because it reads `data.shape` before checking the
  data. And `_cluster.Tree` and `_cluster.pca` ignore keyword arguments:
  `Tree(nodes=[...])` gives an empty tree. Upstream too. (PR #192)

### Other modules

- `Bio/Medline/__init__.py` (`parse`): records separated only by a line
  of six or more spaces still merge, and their PMIDs join as `"1 2"`.
  This is the trade-off of upstream #4557/#4558; a fix needs lookahead.
  Upstream too. (PR #170)
- `Bio/SCOP/__init__.py:287`: a CLA sunid missing from the DES file gives
  a bare `KeyError` (A.12 fixed the HIE case). Upstream too. (PR #165)
- `Bio/PopGen/GenePop/FileParser.py`: calling `FileRecord.start_read()`
  again leaks the old handle and appends the loci twice, and
  `remove_loci_by_position()` sorts and reverses the caller's list in
  place. Upstream too. (PR #140)
- `Bio/PopGen/GenePop/LargeFileParser.py`: the module docstring names a
  class `LargeRecord`; it is `Record`. Moot once the module goes.
  (PR #139)
- `Bio/Graphics/GenomeDiagram/_GraphSet.py`: `new_graph` silently accepts
  `centre=`, though `DEPRECATED.rst` says it went in 1.62. Upstream too.
  (PR #168)
- `Bio/Graphics/GenomeDiagram`: `Diagram.new_track` and `Track.new_set`
  `setattr` any keyword, so a misspelt one does nothing. And
  `FeatureSet.add_feature(colour=..., color=...)` takes whichever comes
  last; that goes with the alias. Upstream too. (PR #168)
- `Bio/Graphics/GenomeDiagram/_Diagram.py` (`write`): a `pathlib.Path`
  filename fails for PDF output (ReportLab): `write(Path("x.pdf"), "PDF")`
  raises `TypeError: Cannot use PosixPath(...) as a filename or file`.
  SVG works. Upstream too. (PR #164)
- `Bio/SVDSuperimposer/__init__.py` (`set`): 1-D arrays raise
  `IndexError` instead of its own mismatch error. Upstream too. (PR #135)
- The module-level deprecation warnings (`Bio.PDB.mmtf`,
  `Bio.Blast.NCBIXML` and `NCBIWWW`, `LargeFileParser`) have no
  `stacklevel=2`, so they point at `Bio`, not the importer. Change them
  together. (PR #139)

### Typing

- `Bio/SeqRecord.py`: `reverse_complement` and `translate` annotate `id`,
  `name`, `description`, `features`, `annotations`, `letter_annotations`
  and `dbxrefs` as `bool`, so the documented
  `record.reverse_complement(id="rc")` is an error downstream. (PR #187)
- `Bio/SeqRecord.py`: the `seq` setter is typed `Seq | MutableSeq`, but
  the runtime accepts `None`. (PR #187)
- `Bio/SeqRecord.py` (`_RestrictedDict.update`): takes only a mapping, so
  `letter_annotations.update(key=...)` type-checks but raises `TypeError`.
  (PR #187)
- `AGENTS.md`'s Typing section does not describe the C extension stubs,
  `stubtest` and its allowlist, or the rule that every extension needs a
  stub, and should say the bound-`TypeVar` rule is for runtime code, not
  `.pyi` files. (PR #191)

### Tests

- The Windows test step is not a gate (`continue-on-error`): 22 to 33
  tests fail there on `main`, mostly SQLite index files still open when
  their directory is removed and BLAST XML checked out with CRLF line
  endings. (PRs #171 and #194)
- `test_PDB_binary_cif.test_byte_swapped_buffers_rejected` fails on
  PyPy 3.10 with msgpack installed; CI's PyPy job has no msgpack. Not
  re-run. (PR #171)
- `Tests/test_PDB_PSEA.py`: with `psea` installed, two tests fail, because
  `run_psea` returns an absolute path; fix it with the `outdir=` and
  `endswith` bugs listed in PR #163. (PR #163)
- `Tests/test_Tutorial.py`: `python test_Tutorial.py --offline` fails,
  because `--offline` is never removed before `unittest.main()`. (PR #166)
- `test_SeqIO_index` and `test_SearchIO_write` write scratch files into
  `Tests/<dir>/`, so two runs at once collide; `test_Affy` does the same
  in `Tests/Affy`, and fails on a read-only checkout. (PRs #166 and #162)
- `test_Align_bigbed.TestAlign_big` reads its generated `.bb` files from
  the current directory. (PR #162)
- `Tests/test_UniProt.py` and `Tests/test_ExPASy_offline.py` have no
  `__main__` block, so running them directly runs nothing. (PRs #166 and
  #164)
- `Tests/support.py`: `DATA` is not `resolve()`d, so a test comparing it
  with a resolved path breaks in a symlinked checkout. (PR #164)
- `Tests/test_NCBI_qblast.py` mocks `NCBIWWW.urlopen` for the whole
  process at import, and its responses are used up in test order.
  (PR #194)
- `Tests/test_align_substitution_matrices.py`: tests 4 to 6 depend on
  state left by earlier tests; `test6_scores` alone fails. (PR #194)
- The `BiopythonExperimentalWarning` filters in `test_phenotype`,
  `test_phenotype_fit` and `test_Tutorial` filter nothing now. (PR #141)
- `Tests/expected_skips.txt`'s header mentions only
  `run_tests.py --check-skips`, not `python -m pytest --check-skips`.
  (PR #193)

### CI, build and tooling

- Codecov rejects every upload ("Token required because branch is
  protected"); there is no `CODECOV_TOKEN`. CI does not fail on it, by
  decision. Adding the secret, allowing tokenless uploads, or removing
  the upload steps is Dan's call. As reported, not re-checked.
  (PRs #153 and #194)
- `.github/workflows/ci.yml`: the style job's checkout keeps its
  persisted credentials; the other ten set
  `persist-credentials: false`. (PRs #153 and #150)
- `.github/workflows/release.yml`: the publish step has no
  `skip-existing: true`, so a partial upload cannot be finished by
  re-running the job. (PR #161)
- `pyproject.toml`: only `_pairwisealigner` declares `depends`, so
  editing `Bio/Cluster/cluster.h`, `_arraycore.h` or
  `Bio/_freethreading.h` does not rebuild the extensions that include
  them. (PR #173)
- The ruff S101 hook's file pattern leaves out `Bio/GenBank/`, which needs
  `noqa` on its 15 kept asserts first. And if RUF100 is ever turned on in
  the main ruff hook, which runs with `--fix`, it would strip all 132
  `noqa: S101` comments; merge the two hooks' rules then. (PR #180)
- `Bio/__init__.py:128-140`: every editable install warns "importing from
  inside the source tree", because it looks for a `biopython.egg-link`
  that PEP 660 installs do not create. (PR #193)

### Documentation

- `Doc/Tutorial/chapter_contributing.rst:14,24,103,119` presents
  upstream's mailing list as this project's. (PR #130)
- `Doc/Tutorial/chapter_introduction.rst:106` says
  `pip install biopython`, and links to the README, which says
  `pip install --pre biopaithon`. (PR #130)
- `Bio/SeqIO/__init__.py`: the format table has typos in the gfa1,
  fastq-solexa and tab descriptions. (PR #132)
- `Doc/Tutorial/chapter_pdb.rst:2140-2155`: says `retrieve_pdb_file`
  saves `pdb1fat.ent` in `.Z` format by default; the default is now
  mmCIF. (PR #131)
- `CONTRIBUTING.rst:131`: the dev install is still
  `pip install -e . --group dev`, where the README now says
  `pip install -e ".[test]"`, with a PyPy caveat. (PR #137)
- `DEPRECATED.rst:176`: the `Bio.SeqIO.TabIO` heading's underline is
  longer than the title. Cosmetic. (PR #134)

## How this document was produced

Seven independent reviews of the tree at `5d6fe8d22`, each assigned one
dimension: build/packaging/CI; the core sequence stack; the I/O and parser
architecture; static typing and API surface; accumulated deprecations and dead
code; performance and the C extensions; and the test suite. All seven are
represented below.

Findings were required to cite file and line, and to be verified by reading or
running code rather than inferred from general knowledge. The four bugs marked
**[reproduced]** in Tier 0 were additionally re-confirmed independently before
being written up here.

Improvements to this document are welcome from any contributor — see
[AGENTS.md](AGENTS.md).
