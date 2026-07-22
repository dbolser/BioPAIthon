# BioPAIthon Improvements Plan

A prioritised, evidence-backed plan for improving this codebase, produced by a
comprehensive review of the tree at the point BioPAIthon forked from Biopython
(`5d6fe8d22`, version `1.88.dev0`).

This is a living document. Each item cites the file and line it was found at so
it can be checked rather than believed. Where a claim was reproduced by running
code, it is marked **[reproduced]**; where it rests on reading the source, it is
marked **[from source]**. Nothing here has been fixed yet.

**Effort** is S (hours), M (days), L (weeks). **Impact** is this fork's judgement
of user-visible value.

---

## Tier 0 — Correctness bugs

These produce wrong answers or destroy diagnostic information. They are cheap to
fix and should go first. Several are silent, which is what makes them serious:
users get a plausible number rather than an error.

### 0.1 `ProteinAnalysis.flexibility()` ignores the residue at the centre of its window **[reproduced]**

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

### 0.2 `Bio.Seq.translate()` silently discards `gap=` for `Seq` inputs **[reproduced]**

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

### 0.3 Derived `SeqRecord`s share mutable annotations with their parent **[reproduced]**

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

### 0.4 `Seq.search()` never matches at the final position **[reproduced]**

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

### 0.5 `as_handle()` swallows `TypeError` raised by its caller's block **[from source, with repro]**

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

### 0.6 `SeqRecord.reverse_complement()` crashes on features with `UnknownPosition` **[from source]**

`Bio/SeqRecord.py:1406-1413`. The sort key catches `TypeError` "expected for
UnknownPosition" and returns `None` — but a key list containing `None` cannot be
ordered, so `list.sort` raises `TypeError: '<' not supported between instances
of 'NoneType' and 'NoneType'`. The `except` clause defeats its own purpose.

`UnknownPosition` is produced by the GenBank/EMBL parsers for `?` locations, so
this fires on real files with no obvious connection to the user's action.

**Fix:** return a total-order sentinel — `(0, int(start))` for known and `(1, 0)`
for unknown — so unknown positions sort stably to the end.
**Effort S · Impact medium**

### 0.7 `SeqIO.index()` and `SeqIO.parse()` disagree on the same file **[from source]**

`Bio/SeqIO/_index.py:219` re-derives the record ID independently of the parser
(`line[marker_offset:].strip().split(None, 1)[0]`), while
`Bio/SeqIO/FastaIO.py:245-249` explicitly handles a bare `>` line. On a FASTA
record with an empty description, `parse` yields `id=''` and `index` raises a
bare `IndexError: list index out of range` with no filename or line number.

`Bio/File.py:243-244` (`raise ValueError(f"Key did not match ({key} vs
{key2})")`) exists solely to police this duplication after the fact. See §1.4
for the structural fix.
**Effort S (this case) · Impact medium**

### 0.8 343 `assert` statements validate untrusted file content **[from source]**

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

### 0.9 Malformed GenBank is a warning where malformed EMBL is an error **[from source]**

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

### 0.10 Two more error handlers report the wrong cause **[from source]**

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

### 0.11 Heap buffer overflow in `_bcif_helper` driven by untrusted file metadata **[reproduced by review]**

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

### 0.12 Out-of-bounds read and a dead validation check in `_aligncore` **[reproduced by review]**

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

### 0.13 `ccealign` leaks ~116 KB per call and segfaults on tuple coordinates **[reproduced by review]**

`Bio/PDB/ccealignmodule.c:610-611,617-618,636-637` contain three `Py_INCREF`s on
references that are immediately stolen (by `Py_BuildValue("[NN]", ...)` and
`PyStructSequence_SetItem`). Measured on 200-residue chains over 100 calls:
**116.5 KB and 261 tracked objects leaked per `run_cealign()` call**. The
reviewer patched out exactly those three `Py_INCREF`s, rebuilt, and re-measured:
**0 objects, ~0 KB per call**, nothing else changed. An all-against-all
`CEAligner` run over 10,000 pairs leaks roughly 1.1 GB.

Also: `:660-661` builds a fresh heap type inside the loop (20 distinct
`CEAlignment` types per call, so `type(r[0]) is type(r[1])` is `False` and the
results will not pickle); `:381-386` calls `PyList_GetItem` with no NULL check,
so passing tuples instead of lists **segfaults**; `:685` discards the
`PyArg_ParseTuple` return; every `PyMem_RawMalloc` at `:107,110,193,196,451` is
unchecked; `:751` is missing `PyMODINIT_FUNC`.

**Fix:** delete the three `Py_INCREF`s; hoist the struct-sequence type to a
module-level singleton created in `PyInit_ccealign`; NULL-check everything.
**Effort S (leak) / M (hardening) · Impact high**

### 0.14 `cpairwise2.rint` writes 8 bytes into a 4-byte stack slot **[reproduced by review]**

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

### 1.1 `Bio.AlignIO` and `Bio.Align` are two complete parser stacks for the same formats

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

### 1.2 Five incompatible format-registration mechanisms, none extensible

`Bio/SeqIO/__init__.py:414-475` and `Bio/AlignIO/__init__.py:153-176` use dicts
of eagerly imported classes; `Bio/Align/__init__.py:4839-4849` derives a module
path from the format string via `importlib`; `Bio/SearchIO/_utils.py:34-63` uses
lazy `(module, class)` string tuples.

Consequences: `Bio.Align` format names must be valid Python module names, which
is *why* it cannot offer `phylip-relaxed` (§1.1) and why it says `tabular` where
SeqIO says `fasta-m10`. Case handling differs — `Align.read(f, "FASTA")` works,
`SeqIO.parse(f, "FASTA")` raises. There is no registration hook at all, so a
downstream package must mutate private dicts. And the eager imports mean
`import Bio.SeqIO` pulls in `Bio.AlignIO` → `Bio.Align` → NumPy: **you cannot
parse a FASTA file without NumPy installed.**

**Plan:** one `Bio/_io_registry.py` with lazy `format_name → "module:Class"`
entries, keyed on an explicit name so subtype names stay free-form. All four
packages resolve through it. Add `register_format()` plus an
`importlib.metadata` entry-point group so plugins work.
**Effort L · Impact high**

### 1.3 The typing story is worse than having no types at all

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

### 1.4 The indexing layer re-implements each parser's record-boundary logic

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

### 1.5 Three parallel BLAST stacks, and the tutorial teaches the superseded one

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

### 1.6 `check_untyped_defs` is off, hiding 2,172 real errors

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

### 1.7 Import cost: ~190 ms warm for `import Bio.SeqIO`, and 77 ms of it is one module body

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

### 1.8 No `__all__` on re-exporting packages, and no stubs for 13 C extensions

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

### 1.9 No C extension ever releases the GIL, so threads give zero speed-up

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

### 1.10 Free-threaded CPython is unsupported and blocked by mutable file-scope state

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

### 1.11 FASTQ parsing costs 4.9× the time and 3.9× the memory of the raw iterator

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

### 1.12 `PairwiseAligner.align()` allocates the full O(n·m) traceback matrix

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

### 2.1 Already-broken or already-empty — no user impact

- `Scripts/xbbtools/xbb_blastbg.py:20-24` and `xbb_blast.py:23,200` import
  `Bio.Blast.Applications`, **which no longer exists**. The xbbtools GUI raises
  `ImportError` on first use, and `MANIFEST.in:14` ships it in every sdist.
- `Bio/HMM/__init__.py` is a 5-line docstring-only file — all four submodules
  were removed in 1.86 — yet `pyproject.toml` still ships `Bio.HMM` as a package.
- `Tests/run_tests.py:69-73` still excludes doctests for five modules removed in
  1.86, and `:96` for `Bio.PDB.Vector`, removed in **1.74**. The whole 45-entry
  list is guarded by `if np is None:` (`:61`) — unreachable, since numpy is a
  hard dependency.
- `Bio/Align/__init__.py:4449,4507`: `# FIXME remove this after 1.87 is out`.
  1.87 is out.
- `pyproject.toml:58` still ships `Bio.Alphabet` six years after removal.

**Effort S · Impact medium** (nothing breaks — it is already broken)

### 2.2 `Bio.PDB.mmtf` — the format's server no longer resolves in DNS

RCSB decommissioned MMTF in July 2024; `mmtf.rcsb.org` does not resolve. The
dependency `mmtf-python` last released 2022-07-06. `Tests/test_mmtf_online.py:28`
calls `get_structure_from_url("4ZHL")` and can never pass again. 563 lines plus a
CI dependency (`ci-dependencies.txt:18`).

**Plan:** one release with a deprecation warning, then remove the subpackage,
both test files, the `pyproject.toml` entry and the CI dependency; repoint
`Tests/test_PDB_internal_coords.py:63-66` at an existing mmCIF fixture.
`Bio.PDB.binary_cif` is the successor.
**Effort S · Impact medium-high**

### 2.3 `Bio.pairwise2` — deprecated eight releases ago, still built and shipped

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

### 2.4 The 1.86 deprecation cohort is ripe for one batched removal

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

### 2.5 `Bio.codonalign` has warned "experimental" for twelve years

`Bio/codonalign/__init__.py:21-25` warns on every import; introduced 2014.
2,644 lines whose last substantive commit was a formatting pass. An import-time
warning that has fired unchanged for twelve years trains users to filter *all*
Biopython warnings, which then masks the real deprecations above.

**Plan:** decide and record. Either drop the warning and commit to the API, or
migrate alignment-building onto `Bio.Align.CodonAligner` and keep only
`cal_dn_ds` and `mktest`, which are what people actually use.
**Effort M · Impact medium**

### 2.6 Python-2-era shims, one of which silently breaks a documented promise

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

### 2.7 Unhashable value objects

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

### 3.1 Supply-chain exposure and a shell-injection sink

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

### 3.2 CI builds 15 wheels per run; none is installable, tested, or published

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

### 3.3 Caching defeats the "test against latest dependencies" intent

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

### 3.4 Two of five supported Pythons are tested; MySQL is started but unusable

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

### 3.5 Three CI systems that have already drifted

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

### 3.6 Tooling configuration is split and already inconsistent

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

### 3.7 No optional-dependency extras; unbounded `numpy`

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

### 3.8 The sdist ships ~108 MB of test data

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

### 4.1 Five test modules fail if `PYTHONWARNINGS` is set **[reproduced by review]**

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

### 4.2 The suite only runs from inside `Tests/`

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

### 4.3 The runner reports 501 "tests" that are really 3,050 test cases, and discards their output

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

### 4.4 The previously documented test command did not work

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

### 4.5 Skips are silent and unbounded, and a live API key is committed

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

### 4.7 Dead and unsound runner machinery

- `run_tests.py:62-104` — a 43-line `if np is None:` block excluding 40 modules
  from doctests. numpy is a hard dependency, so this is unreachable (also §2.1).
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

### 4.8 Test fixtures dominate the repository

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

## Suggested sequencing

1. **§0.11–0.14 before anything else.** These are memory-safety defects in C
   code, three of them reproduced as segfaults, and §0.11 is reachable from a
   malformed input file. All four are S-effort. §0.11 in particular should be
   treated as a security fix, not a bug fix.
2. **The rest of Tier 0.** Mostly S-effort, and it is where users are actively
   getting wrong answers. §0.1–0.4 are reproduced and unambiguous. Fix the
   `Tests/test_ProtParam.py` fixture that currently pins §0.1 in place.
3. **Make the suite trustworthy before relying on it** — §4.1 (five modules fail
   under `PYTHONWARNINGS`), §4.4 (the runner breaks in a clean venv) and the
   committed NCBI API key in §4.5. All S-effort, and everything below is
   verified by running the tests, so this comes first among the non-urgent work.
4. **§3.1 and §3.6** — supply-chain pinning and making mypy actually see numpy
   are both S-effort and make every later change safer to land.
5. **§1.6 then §1.3** — turn on `check_untyped_defs` with a ratchet baseline
   before annotating, so the annotation work is checked as it lands.
6. **§1.7 (import cost)** and **§1.11 (FASTQ hot path)** are the highest
   user-visible-value M-effort items, and §1.11(a) is nearly free.
7. **Tier 2 cruft** can proceed in parallel with anything; start with §2.1,
   which breaks nothing because it is already broken.
8. **§4.2 and §4.3** — removing the cwd dependency and moving to a real runner
   unlock parallelism and per-test reporting, which makes every later item
   cheaper to verify.
9. **§1.1, §1.2 and §1.4** are the L-effort structural items. Do §1.4 first —
   it is the smallest, and it fixes §0.7 properly.
10. **§1.9 then §1.10** — releasing the GIL is what makes the free-threading
    work worthwhile, but §1.10's static-state cleanup is a prerequisite for
    doing §1.9 safely under a free-threaded build.

---

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
