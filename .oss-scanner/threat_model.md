# Threat model

BioPAIthon is a fork of Biopython. The importable package is still `Bio`
(plus `BioSQL`). `AGENTS.md` has the project's own rules for changes.

## What this project does and where untrusted input enters

`Bio` is a Python library for computational molecular biology. Its main job
is reading biological file formats. Most of those files are downloaded from
public databases (NCBI, UniProt, the PDB and others) or written by other
programs. Treat every file, and every server response, as untrusted. The
calling program is trusted; the values it passes on from files are not.

Untrusted input enters through:

- File parsers: `Bio.SeqIO` (FASTA, FASTQ, GenBank, EMBL, SwissProt, UniProt
  XML, SFF, ABI, 2bit, SnapGene, ...), `Bio.AlignIO` and `Bio.Align`
  (Clustal, Stockholm, MAF, bigMaf, bigBed, PSL, SAM, ...), `Bio.SearchIO`
  and `Bio.Blast` (BLAST XML, HMMER, ...), `Bio.PDB` (PDB, mmCIF, BinaryCIF,
  PDBML), `Bio.Phylo` and `Bio.Nexus` (Newick, NEXUS, PhyloXML, NeXML),
  `Bio.GenBank`, `Bio.SwissProt`, `Bio.Entrez` (NCBI XML), `Bio.motifs`,
  `Bio.Cluster.read`, `Bio.bgzf`, and the indexes (`SeqIO.index`,
  `SeqIO.index_db`).
- DTDs and schemas named in a document: `Bio.Entrez.Parser` and the
  `Bio.Blast` XML parser resolve them from the copies in `Bio/Entrez/`, and
  `Bio.Entrez.Parser` can also download them from NCBI into a local cache.
- Network clients, whose responses are parsed like files: `Bio.Entrez`,
  `Bio.ExPASy`, `Bio.Blast.qblast`, `Bio.PDB.PDBList`,
  `Bio.PDB.alphafold_db`, `Bio.TogoWS`, `Bio.UniProt`, `Bio.KEGG.REST`.
- The C extensions, listed in `pyproject.toml`. Four are fed file content:
  - `Bio.Align._aligncore` (`Bio/Align/_aligncore.c`): the printed
    alignment rows of most `Bio.Align` parsers.
  - `Bio.PDB._bcif_helper` (`Bio/PDB/bcifhelpermodule.c`): BinaryCIF
    column decoding.
  - `Bio.SeqIO._twoBitIO` (`Bio/SeqIO/_twoBitIO.c`): 2bit sequence data.
  - `Bio.Nexus.cnexus` (`Bio/Nexus/cnexus.c`): NEXUS comments and quoting.

  The other eight take Python arguments, which often hold data read from a
  file: `Bio.Align._pairwisealigner`, `Bio.Align._codonaligner`,
  `Bio.Align._alignmentcounts`, `Bio.Align.substitution_matrices._arraycore`,
  `Bio.motifs._pwm`, `Bio.Cluster._cluster`, `Bio.PDB.kdtrees`,
  `Bio.PDB.ccealign`.

## Components that matter most / least

- Most: the C extensions, above all the four fed file content, and the
  parsers of downloaded files.
- Less: the deprecated modules in `DEPRECATED.rst` (`Bio.PDB.mmtf`,
  `Bio.Blast.NCBIWWW`, `Bio.Blast.NCBIXML`,
  `Bio.PopGen.GenePop.LargeFileParser`), and `Scripts/` (examples, not
  installed).
- Out of scope: `Tests/` and `Doc/`, and bugs in dependencies (numpy, scipy
  and the rest). There is no vendored third-party code.

## How to exercise it

- `/src` is an editable install, so edits to Python files under `/src` take
  effect at once. After changing C code, rebuild in place with
  `cd /src && pip install --no-deps --no-build-isolation -e .`. The image's
  `CC` adds `-g -fno-omit-frame-pointer`.
- `python-asan` is `python` with a second copy of `Bio` from `/opt/asan`,
  whose extensions are built with ASan and UBSan. Use it for anything that
  touches C. It runs `python -P`, so `/src/Bio` cannot shadow `/opt/asan`.
  That also leaves a script's own directory off `sys.path`, so run tests
  under it with pytest or `run_tests.py`, not as `python-asan test_X.py`.
  It passes `--capture=sys` to pytest, so a sanitizer report is printed
  even when it aborts the run. After changing `/src`, refresh that copy with
  `rm -rf /src/build && CC="$CC -fsanitize=address -fsanitize=undefined" pip install --no-deps --no-build-isolation --upgrade --target /opt/asan /src`.
- Tests: `cd /src/Tests && python -m pytest --offline test_SeqIO.py`, one
  module at a time. The whole offline suite is
  `python -m pytest --offline -n 2`: about 3 minutes, and about 10 under
  `python-asan`. Run it from `Tests/`, not `/src`.
- One test, `test_distribution_name.py::TestDistributionName::test_not_listed_twice`,
  fails under `python-asan` and when run from `/src`. Both see a second
  `biopaithon` install: `/opt/asan`, or the build's `biopaithon.egg-info`.
- `Tests/<Format>/` (`Tests/GenBank`, `Tests/PDB`, `Tests/Nexus`,
  `Tests/TwoBit`, ...) hold real sample files, good seeds for malformed
  ones.
- Known open bugs are listed in `IMPROVEMENTS.md`. Please do not report
  those again.

## How you rate severity

- Critical: code execution from parsing a file or a server response. The
  public precedent is `eval()` in the `Bio.Nexus` NEXUS parser (fixed in
  Biopython 1.88). Also memory corruption from parsing a file where code
  execution is plausible, such as a write whose place or content the file
  controls.
- High: other memory corruption from parsing a downloaded file (out-of-bounds
  write or read, use after free, double free). Also a file that makes `Bio`
  read, write or fetch something it was not designed to. The public
  precedent is CVE-2025-68463 in `Bio.Entrez.Parser` (fixed in 1.87): an XML
  document could name a DTD or schema at any URL, and the parser fetched it.
- Medium: memory corruption reachable only through unusual API arguments.
  Resource exhaustion from a small file: a hang, or an allocation far larger
  than the file.
- Low: resource use that grows with the size of the input. Memory that leaks
  on every call into a C extension.
- Not a vulnerability: a Python exception (`ValueError` and the like) on
  malformed input. An interpreter crash (segfault, abort) or a hang is still
  a bug worth reporting, whatever its severity.
- Silently wrong scientific results, from malformed or valid input, are not
  security issues, but are welcome as bug reports.

## Reports

- A minimal reproducer: a short Python snippet and the smallest input file
  that triggers it.
- The file and line at fault, and the `python-asan` output for memory bugs.
- A patch with a regression test in `Tests/`, written as a
  `unittest.TestCase` (pytest collects nothing else), that fails before the
  fix and passes after. Build data paths from `support.DATA`.

## Anything to leave alone

- Network tests (those that call `requires_internet.check()`), and live
  behaviour of the network clients. The scan has no network.
- Small bugs in the deprecated modules listed above.
- Unpickling: loading a pickle runs code by design.
- The `BiopythonWarning` about importing from inside the source tree. This
  editable install sets it off.
- LeakSanitizer, which `python-asan` turns off: its reports at exit are
  dominated by CPython itself (see `.github/workflows/ci.yml`).
  `Tests/test_C_extension_memory.py` measures leaks instead.
