# Code-smell review

Reviewed on 2026-10-17 at commit `051cd495c`.

This review focuses on maintainability and test design. Vulnerability analysis
is out of scope. It is a targeted structural review of the APIs, parsers,
registration machinery, and test infrastructure, rather than an exhaustive
review of every module. No source code was changed during the review.

The seven findings below are refactoring opportunities, not seven confirmed
runtime defects. Surrounding code and relevant history were checked before
making these recommendations. Changes should preserve the public API and be
made in focused patches with appropriate tests.

## 1. Duplicated alignment input preparation

**Status: fixed** in PR #227 (2026-10-20). `PairwiseAligner` and
`CodonAligner` each now prepare their input in a `_prepare_sequences` helper
shared by `score()` and `align()`.

[`PairwiseAligner.align()` and `.score()`](Bio/Align/__init__.py) separately
implement sequence conversion, alphabet construction, and strand handling
(around lines 4527 and 4584 at the reviewed commit). `CodonAligner` repeats
37 identical preprocessing lines between its two methods, confirmed by
comparing their syntax trees.

This creates multiple places to update when supported input representations or
conversion rules change. Scoring and producing an alignment should share the
same preparation rules.

**Recommendation:** Extract a private preparation helper per aligner, with
tests ensuring scoring and alignment accept the same input representations.
Keep the original sequence objects available for the returned alignment.

## 2. Distribution tests depend on installed packages

**Status: fixed** in PR #225 (2026-10-19). The tests now query fabricated
metadata records with `sys.path` restricted to them and the standard library,
and cover Biopython both present and absent.

[`Tests/test_distribution_name.py`](Tests/test_distribution_name.py), around
lines 45–63, assumes upstream Biopython is absent and exactly one BioPAIthon
metadata record is visible. The reviewed checkout exposes old
`biopython.egg-info`, current `biopaithon.egg-info`, and installed distribution
metadata. Those assumptions caused all five reported failures in the offline
test run.

The implementation deliberately lets a real Biopython distribution answer
first. A test that assumes it is absent therefore depends on the developer's
environment rather than solely on the behavior under test.

**Recommendation:** Use controlled metadata fixtures for unit tests and an
isolated installation for integration tests. Cover both the presence and
absence of upstream Biopython explicitly.

## 3. The format registry's abstraction leaks into its callers

[`FormatRegistry`](Bio/_io_registry.py), beginning around line 132, combines
dictionary storage, lazy imports, discovery, provenance, and synchronization.
[`SeqIO.register_format()`](Bio/SeqIO/__init__.py), around lines 1388–1435,
must bypass its methods with `dict.get`, inspect `plugins`, and coordinate
another lock. Tests similarly bypass normal copying to avoid resolving
handlers; see [`stored()` and `restore()`](Tests/test_register_format.py).

Callers need to understand both the raw storage and the registry's lifecycle.
That makes changes to discovery or registration harder to contain within the
registry itself.

**Recommendation:** Add explicit internal operations for inspecting entries,
taking snapshots, and registering related roles while retaining the existing
mapping interface. Preserve the deliberate rule that imports happen outside
locks.

## 4. PAML execution depends on process-wide working-directory state

[`Paml.run()`](Bio/Phylo/PAML/_paml.py), around lines 85–132, changes the
caller's working directory while constructing control files and launching a
program. This couples path interpretation to execution order and complicates
embedding the library in a larger application.

PAML's relative-path requirements explain the existing design, but the child
process can have its own working directory without changing the parent's.

**Recommendation:** Resolve paths explicitly and pass `cwd=` to the
subprocess. Construct control-file paths explicitly as well. Keep the existing
relative-path behavior covered by compatibility tests.

## 5. Alignment slicing invokes a parser method on an unrelated object

[`AlignmentIterator.__getitem__()`](Bio/Align/interfaces.py), around lines
220–226, retrieves `_read_header` from the parser class, then calls it with an
`Alignments` list as `self`.

This creates an implicit constraint: header methods must work without their
actual parser instance. An otherwise ordinary refactor that delegates header
work to another instance method would need to account for this alternate
receiver.

**Recommendation:** Introduce an explicit header-data representation or
copying hook so parser extensions can use ordinary instance methods safely.
Preserve the metadata and other header attributes exposed on materialized
alignment collections.

## 6. SAM parsing combines responsibilities and repeats CIGAR processing

[`AlignmentIterator._read_next_alignment()`](Bio/Align/sam.py), beginning
around line 566, spans 317 lines. It handles fields, optional tags, clipping,
coordinates, sequence reconstruction, and object construction. Its branches
with and without an MD tag separately interpret CIGAR operations.

A change to a shared CIGAR rule requires checking both paths, while the size of
the method makes their common behavior harder to identify and test separately.

**Recommendation:** Extract shared tokenization first, then isolate coordinate
construction and MD-based reconstruction. Preserve the distinct information
available in the two paths and test both against the existing SAM fixtures.

## 7. PAML control-file parsing has copied rules that have drifted

[`CODEML`](Bio/Phylo/PAML/codeml.py),
[`BASEML`](Bio/Phylo/PAML/baseml.py), and
[`YN00`](Bio/Phylo/PAML/yn00.py) repeat comment removal, assignment splitting,
numeric conversion, and option resetting in their `read_ctl_file()` methods.
The implementations already differ in assignment splitting and
scientific-notation detection.

These shared lexical rules are maintained independently even though the
program-specific options are the main reason for separate readers.

**Recommendation:** Share lexical parsing and scalar conversion, retaining
program-specific handling for options such as `NSsites` and `model`. Check
existing behavior before changing accepted syntax or exception types.

## Validation and limitations

The offline suite was run from the repository root under Python 3.14.6:

```bash
.venv/bin/python -m pytest --offline -q
```

Pytest reported:

```text
5 failed, 4329 passed, 119 skipped, 20 warnings, 2257 subtests passed in 212.31s (0:03:32)
```

The five reported failures were confined to two distribution-test methods,
including three subtest failures:

- `TestDistributionName.test_biopython_resolves_to_biopaithon`: three
  case-variant subtests reported `AssertionError: '1.88.dev0' != '1.89.dev0'`;
  the method also reported `AssertionError: 'biopython' != 'biopaithon'`.
- `TestDistributionName.test_not_listed_twice` reported
  `AssertionError: 2 != 1`.

**Status:** those five failures no longer exist. At `4ca8354ba` (2026-10-22)
`Tests/test_distribution_name.py` passes, after PR #225.

The complete session log is at
`/tmp/biopaithon-audit-baseline.log`. This is a temporary local artifact and is
not included in the repository.

Several optional dependencies and external programs were unavailable, and
network tests were skipped by the offline setting. The run therefore does not
establish coverage of all optional functionality. Other Python versions and
platforms were not tested in this review.
