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

### [#5386](https://github.com/biopython/biopython/pull/5386) — Peter Cock (peterjc)

The ruff settings lived only in the arguments of the pre-commit hook, so a
bare `ruff check` used ruff's defaults and disagreed with CI (321 errors
against none). The pull request moves them into `[tool.ruff.lint]` in
`pyproject.toml`, and with them upstream issue
[#2828](https://github.com/biopython/biopython/issues/2828)'s isort settings
into `[tool.isort]`. This closes the ruff item of `IMPROVEMENTS.md` §3.6. It
changes no code and no behaviour of the library.

**The three commits taken** (`70e793df6` adds the `check-toml` pre-commit
hook, `cf6b88a98` makes `setup.py` executable to match its shebang,
`5b767a9f0` moves the settings; cherry-picked with `-x`). They applied
without conflicts; git's 3-way merge handled the context that differs here,
which is the fork's extra S101 hook in `.pre-commit-config.yaml` and the
renamed project URLs in `pyproject.toml`. The `[tool.isort]` block does
nothing here, since the fork runs no isort hook, but is kept so the file
matches upstream at the next sync.

**Left behind, and why.** `c6072cec8`, "Ran isort 9.0.2, and black
24.10.0", a sweep that reorders imports and moves comments in 23 files. It
does not apply here, because the fork has changed several of those files,
including `Bio/pairwise2.py` and `Bio/Restriction/__init__.py`, and
`Scripts/xbbtools/nextorf.py` no longer exists. Much of it is comment
shuffling of no value, such as moving "Import these explicitly to avoid
flake8 F405 below" above the `import *` line it explains. The hook does not select ruff's
import-order rules (`I`); if the fork wants them, it can run
`ruff check --select I` itself and check each hit for intent.

**Added by the fork:** `AGENTS.md` no longer lists the ruff options and
points at `pyproject.toml` instead; `IMPROVEMENTS.md` §3.6 and `NEWS.rst`
record the change.

**Checked:** with ruff 0.15.20 on the 589 tracked `.py` files, a bare
`ruff check` now gives "All checks passed!", and `ruff check --isolated`
(ignoring the config) still fails, so the config is what is being read. The
fork's second hook (`--select=S101` over the parser packages) still enables
only S101, as `--show-settings` shows, and passes. Peter Cock is already in
`CONTRIB.rst`.

**Upstream status:** as of 2026-11-10 the pull request is open, not a draft,
with no reviews; peterjc rebased it on 2026-11-02. Issue #2828 is open. The
licence box is ticked. If upstream merges a different shape, expect a small
conflict in `pyproject.toml` at the next sync.

### [#5380](https://github.com/biopython/biopython/pull/5380) — Shubhodeep Chatterjee (shubhodeepcx)

Answers upstream issue
[#4652](https://github.com/biopython/biopython/issues/4652), open since
2024: a user asked how to find the residues near a residue with
`Bio.PDB.NeighborSearch`, and peterjc said a proper example in the
Tutorial would be welcome. The "Determining atom-atom contacts" section of
`Doc/Tutorial/chapter_pdb.rst` only named the class. The pull request adds
a doctested example on `Tests/PDB/1A8O.pdb`: a `NeighborSearch` over all
atoms, `search(coord, 5.0, level="R")` around the CA of residue A:152, and
`search_all(4.0, level="R")` for all residue pairs. Docs only; no code
changes.

**The single commit taken** (`dc95e1e96`, cherry-picked with `-x`), which
applied without conflicts, including the author's own `CONTRIB.rst` line.

**Left behind:** nothing.

**Added by the fork:** the text listed `"C"` and `"M"` as the other levels;
a fork commit adds `"S"` (structure), which `search` returns. (`search_all`
also accepts `"S"`, but on a single structure it always returns no pairs,
since both atoms of every pair share the one structure; the Tutorial text
only describes `search`.) A second fork commit, after review, says that
the CA-centred search finds the residues near that one atom, not near the
whole residue, and adds a doctested search around every atom of the
residue. On 1A8O that adds ARG 154, whose N atom is 3.34 Å from the O of
ASP 152 but more than 5 Å from its CA.

**Checked:** the example's outputs (ASP 152; residues 151, 152, 153 and
the waters 1001, 1027, 1035 and 1054; 525 pairs) are the same on this
fork, on upstream master and on Biopython 1.88. `Tests/test_Tutorial.py`
extracts and runs the new block, and passes.

**Upstream status:** as of 2026-11-10 the pull request is open and not a
draft, with one commit and no reviews; the last activity was on
2026-10-18, when peterjc asked whether AI was used and the author pointed
to his answer on #5378. The licence box is ticked. The same `CONTRIB.rst`
line is in #5378, so adopting that too would give a trivial conflict.

### [#5194](https://github.com/biopython/biopython/pull/5194) — Gert Hulselmans (ghuls)

`Bio.motifs.minimal.read` required a "Background letter frequencies" section
and raised `ValueError` without one, although the MEME format calls it
recommended, not required, and says the MEME Suite assumes uniform
frequencies when it is missing. Files from factorbook.org leave it out. The
pull request makes `_read_background` default to 0.25 for each letter, stop at
the first `MOTIF` line if it finds no background section, and hand that line
back to the main loop, so the motif is not lost. It reads the handle only by
iteration, so it works on non-seekable streams.

**The single commit taken** (`56951256b`, cherry-picked with `-x`). It adds
`Tests/motifs/minimal_test_no_background.meme` (the existing
`minimal_test.meme` without its background section) and a test of it. One
conflict, in `Tests/test_motifs.py`: the commit deletes a redundant second
`open` of `minimal_test.meme`, which this fork had already moved to
`support.DATA`; resolved by deleting it, and the new test's `open` was
pointed at `support.DATA` too, as this fork's tests do. The PR body still
describes its first revision, which rewound the handle with `tell` and
`seek`; that was dropped under review (a handle may be `sys.stdin`) and is
not in the commit taken.

**Left behind:** nothing.

**Added by the fork:** with the commit alone, a file that ends right after
its `ALPHABET` line raised `UnboundLocalError`, because the `for` loop in
`_read_background` never ran and `line` was never set. A fork commit sets
`line = ""` before the loop, so that file gives an empty record, as any
file without motifs now does, with a test. Also checked, not changed: a
background section placed after a motif (outside the format) used to make
the parser skip the motifs before it, silently; with this change all motifs
are kept and that late section is ignored, giving a uniform background.
`NEWS.rst` says so.

**Checked:** on this fork's main the new no-background test fails with the
`ValueError` and the end-after-`ALPHABET` test fails too; with the branch
all of `Tests/test_motifs.py` passes and the `Bio.motifs.minimal` doctest
passes. Gert Hulselmans is already in `CONTRIB.rst`.

**Upstream status:** as of 2026-11-10 the pull request is open and not a
draft, with one commit and no approving review; the last activity was on
2026-04-16, when peterjc asked mdehoon to look again. The licence box is
ticked.

### [#5255](https://github.com/biopython/biopython/pull/5255) — Arata Synqa Bando

Fixes upstream issue
[#729](https://github.com/biopython/biopython/issues/729), open since 2015:
`kegg_conv` rejected the `/conv/<target_db>/<dbentries>` form that KEGG's
own documentation uses as its examples, such as
`kegg_conv("ncbi-proteinid", "hsa:10458+ece:Z5100")`, and also
`kegg_conv("compound", "chebi")`, because a typo had `"glycan"` where
`"chebi"` belonged. The pull request brings each `Bio.KEGG.REST` function
up to date with the current KEGG REST API, one commit per operation.

**The four commits taken** (`49799cf7f` info, `46f5a7c99` find,
`1600e60a2` get, `2c21f1725` link, cherry-picked with `-x`), which applied
without conflicts. They update the comments listing KEGG's databases, add
the `"nop"` option to `kegg_find` and `"conf"` to `kegg_get`, and give
`kegg_link` the taxonomy options `"species"` to `"phylum"` while allowing
its RDF options only for the `drug`, `atc` and `jtc` targets, as KEGG now
does (`link/pathway/hsa:10458/turtle` is HTTP 400, `link/atc/D01441/turtle`
answers). Two of them add online tests.

**Left behind, and why.**

- `97ae71e7f`, the conv commit. It keyed the entries form on the target
  only, so `kegg_conv("eco", "ncbi-geneid:948364")`, which the KEGG
  documentation allows and the server answers, still raised `ValueError`,
  and it kept the `option` argument but ignored it silently, documenting
  it as deprecated. The fork's own commit redoes it.
- `b48bb563b`, the list commit. Its new `option` argument and the removal
  of `/list/module/<org>` (HTTP 400 today) are right and are redone in the
  fork's commit under the same keyword name, but it also lowered the entry
  limit from 100 to 10, which refuses requests the server answers in full.
- `93b25fbe2`, which drops the online test of `/list/organism` (HTTP 400
  today). Correct, but redone in the fork's commit with a test of
  `/list/genome/animals` in its place.
- `f106d1c6b`, the offline tests. They went into `Tests/test_KEGG.py`, which
  imported `test_KEGG_online` at module level; that module skips itself
  offline, which took the KEGG parser tests with it. The fork writes its own
  offline module instead.
- `b11f2d32c`, the author's `CONTRIB.rst` and `NEWS.rst` lines; the author
  is added to the fork's `CONTRIB.rst` instead.

**Corrections on top**, in the fork's own commit. `kegg_conv` treats any
source containing `:` as entries, convertible to any database including an
organism or `genes`; converts to an outside database from any source,
which covers unprefixed KEGG entries such as `kegg_conv("chebi",
"C00001")`; and accepts `ncbi-proteinid` and `chebi`. KEGG withdrew conv's
RDF output and `ncbi-gi` (HTTP 400 for the database form, an empty answer
for entries such as `ncbi-gi:4507949`, even alongside other entries), so
those now raise `ValueError`, where the pull request ignored an option
silently. So does a conversion between a gene and a chemical database,
such as `kegg_conv("pubchem", "hsa")` or `kegg_conv("uniprot",
"compound")`, which KEGG answers with HTTP 400. `kegg_list` gains
`option=` for the two listings KEGG now subsets, such as
`/list/genome/animals`, and its error for `kegg_list("genome", "animals")`
names the `option=` form; it raises `ValueError` for `module` with an
organism, and applies its 100-entry limit to a `+`-joined string as well
as to a list: the server accepts 150 entries but silently returns only the
first 100. The `kegg_info` docstring pointed at `/list/organism` and now points
at `/list/genome`. The online tests use `ncbi-proteinid` in place of
`ncbi-gi`, and gain tests of the conv calls above. The find commit allowed
`"nop"` only for `compound` and `drug`, as if it were a chemical field; it
turns off keyword splitting for any keyword search (`find/genes/shiga/nop`
and `find/pathway/glycolysis/nop` answer), so the fork allows it for any
database. The link commit's RDF allowlist left out `dr`, KEGG's
abbreviation of `drug`, which `kegg_link("dr", "atc", "turtle")` used and
the server still answers; the fork adds it.

**Checked.** Every URL the new offline `Tests/test_KEGG_REST.py` accepts,
and the 400s it guards against, were requested from rest.kegg.jp on
2026-11-10. Nine of its thirteen tests fail on this fork's main. The online
`Tests/test_KEGG_online.py` passes in full against the live server (49
tests), as does the offline suite.

**Upstream status:** as of 2026-11-09 the pull request is still a draft,
opened in July 2026, with no reviews or comments and no activity since
2026-08-07. The licence box is ticked, and the author added their own name to
`CONTRIB.rst`. Issue #729 is open and upstream master has no fix; upstream
master still sends `/list/organism` in its online tests.

### [#5233](https://github.com/biopython/biopython/pull/5233) — Xinxiang Wang (yiyabo)

Fixes upstream issue
[#4291](https://github.com/biopython/biopython/issues/4291), reported in 2023
by jmzhuzhonghua. `open_database(driver="pymysql")` connected, but BioSQL
listed its MySQL drivers by name in four places and PyMySQL was in none of
them. So `load_database_sql` raised
`ValueError: Module pymysql not supported by the loader.`, the generic
`Adaptor` was used rather than `MysqlConnectorAdaptor`, `last_insert_id` went
through `Generic_dbutils`, and `SET sql_mode='ANSI_QUOTES'` was never run,
which makes MySQL 8 reject the double-quoted `"rank"` column with error 1064
when records are loaded. The pull request adds `"pymysql"` to the
`ANSI_QUOTES` check, the loader's list, the adaptor map and the `_dbutils` map.
Every change is keyed on the driver name `"pymysql"`, so other drivers are
unaffected.

**The single commit taken** (`a70f5de1f`, cherry-picked with `-x`). Its
`NEWS.rst` hunk added the author to the contributor list of upstream's 1.88
release, which this fork keeps as inherited history, so it was dropped while
resolving the conflict. The `CONTRIB.rst` hunk is kept: the author ticked the
box and added their own name.

**Left behind, in a fork commit on top:** the commit also added `"pymysql"`
to the drivers whose connection keywords `open_database` renames to `db` and
`passwd`, as MySQLdb wants. PyMySQL takes the DB-API names `database` and
`password` and deprecates `db` and `passwd`, so with that hunk every
connection raised two `DeprecationWarning`s (440 over
`test_BioSQL_pymysql.py`), and under `-W error::DeprecationWarning` the module
could not connect at all and skipped. PyMySQL is left in the existing `else`
branch, which already turns `db` and `passwd` into the DB-API names, and a
comment there says why.

**Added by the fork:** `Tests/test_BioSQL_pymysql.py`, the common BioSQL tests
with `DBDRIVER = "pymysql"`, mirroring `test_BioSQL_MySQLdb.py`, with its
`expected_skips.txt` line. The module also has three tests of its own, which
use its own settings rather than `common_BioSQL`'s, so they reach MySQL
through PyMySQL even in a full run (see the caveat below): one checks the
adaptor and the `_dbutils` entry, one that `ANSI_QUOTES` is set, and one
loads the schema and GenBank records into a database of its own and reads
them back. Each fails when the hunk it covers is reverted.
`pymysql` is installed next to `mysqlclient` in the CI job that runs a MySQL
server, and listed among the alternative BioSQL drivers in `pyproject.toml`.
The pull request had no test.

**Checked**, against a local MySQL 8.0 server in Docker with PyMySQL 1.2.4:
on this fork's main all 72 tests of the new module fail or error with the
loader `ValueError`; with the upstream commit alone they pass with 440
warnings; with the fork's correction they pass with
`-W error::DeprecationWarning`. `test_BioSQL_MySQLdb.py`,
`test_BioSQL_mysql_connector.py` and `test_BioSQL_sqlite3.py` pass, each run
on its own.

**Caveat, not caused here:** run in one pytest session, as CI runs the
suite, the server-backed BioSQL modules do not reach their server. They all
share `common_BioSQL`'s test classes, which read the driver from that
module's globals, and each module's `check_config` overwrites them at import,
so every module's tests run with whichever driver was imported last,
`sqlite3`. Locally all nine BioSQL modules together pass 289 tests in about
ten seconds, while `test_BioSQL_pymysql.py` alone takes nearly two minutes
against MySQL. So until that is fixed, in CI only the module's three own tests
reach MySQL through PyMySQL; the shared tests do not, and neither do
MySQLdb's.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with one commit, no reviews and no activity since 2026-08-06. The
licence box is ticked. Issue #4291 is open. An earlier one-line fix,
[#5191](https://github.com/biopython/biopython/pull/5191) by musashifolk,
which set `ANSI_QUOTES` for PyMySQL only, was closed unmerged on 2026-04-11;
peterjc wrote there that adding `"pymysql"` to the drivers would be welcome.
[#1882](https://github.com/biopython/biopython/issues/1882), the same `rank`
error reported in 2018, is closed.

### [#5318](https://github.com/biopython/biopython/pull/5318) — NinjagoJaden1

Fixes upstream issue
[#1016](https://github.com/biopython/biopython/issues/1016), filed by peterjc
in 2016 from a RedMine report by Wim De Smet. `EmblScanner` threw away every
`DT` line, so a record read as `"embl"` or `"imgt"` had no
`annotations["date"]`, and writing it as GenBank silently put `01-JAN-1980`
on the `LOCUS` line. The pull request passes the first word of each `DT`
line to `consumer.date()`, so the last one, normally "Last updated", wins:
the pragmatic approach peterjc outlined on the issue.

**The three commits taken** (`57c2cb0`, `006a07d`, `a32ec36`,
cherry-picked with `-x`), with nothing left behind. They applied without
conflicts and have no `NEWS.rst` or `CONTRIB.rst` hunk. The third adjusts two
existing tests for the new annotation: `test_SeqIO_Insdc` deletes `date`
before an EMBL round trip, because the EMBL writer writes no `DT` lines, and
`test_SeqIO_features` skips `date` when comparing a GenBank and an EMBL copy
of the same record, whose dates differ (`21-JUL-2008` against
`14-NOV-2006`). The author ticked the combined "added my name, am listed
already, or do not wish to be listed" box but added no name and is not
listed, so is not added to `CONTRIB.rst`.

**Correction on top**, in the fork's own commit. `data.split()[0]` raised
`IndexError` on a `DT` line with no date, which the parser used to ignore.
Such a line is now skipped with a `BiopythonParserWarning`, as a malformed
`DR` line is, with a test. This guard is worth offering upstream on #5318.

**Checked.** `test_embl_date` fails on main with `KeyError: 'date'`, and the
fork's `test_embl_empty_date` fails there too (no warning); both pass here,
as does the whole offline suite. The pull request's record now writes
`20-APR-1999` on the GenBank `LOCUS` line. `Tests/EMBL/101ma_no_coords.embl`,
whose `DT` line holds `1997-12-13`, keeps that string as its date, and the
GenBank writer warns about it and falls back to `01-JAN-1980`, as for any
malformed date.

**Not done: round-tripping the date.** The EMBL writer still writes no `DT`
line, so an EMBL round trip loses the date. peterjc asked on the pull
request whether it could be round-tripped; a `DT` line also carries release
information that the annotation does not hold, so that would be its own
change.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with no activity since 2026-09-18. peterjc said it "seems fine",
asked for the pull request template, for black to pass, and whether AI tools
were used, then reported the test failures the third commit fixes and asked
about the round trip above. The description now follows the template and all
upstream checks pass, but neither question has a reply. The licence box is
ticked. Issue #1016 is still open and upstream master has no fix.

### [#5362](https://github.com/biopython/biopython/pull/5362) — Zehao Wang

Fixes upstream issue
[#5329](https://github.com/biopython/biopython/issues/5329), reported by
Aniruddha Adak. `Seq({}, length=N)` with `N > 0` built a
`_PartiallyDefinedSequenceData` with no segments at all, breaking that
class's assumption that at least one region is defined: adding it to a
sequence on either side raised `IndexError`, `* 2` raised
`UnboundLocalError`, and printing it said "only partially defined". An
empty dictionary carries the same information as `None`, so `Seq.__init__`
now gives it `_UndefinedSequenceData(length)`, and it behaves exactly like
`Seq(None, length=N)`, `repr` included. `Seq({}, length=0)` takes the
earlier `length == 0` branch and is still `Seq('')`. The pull request adds a
doctest to `Seq.__init__` and three regression tests to `Tests/test_seq.py`,
which all fail on this fork's main before the fix.

**The single commit taken** (`e3d686082d`, cherry-picked with `-x`), which
applied cleanly with line offsets only. Nothing was left behind and there
are no corrections on top. The `if not data:` check sits after the loop over
the segments, which never runs for an empty dictionary, rather than before
it; that is harmless and is kept as upstream wrote it. The author ticked the
box for being listed already, through #5351, and is in `CONTRIB.rst` from
the adoption of #5356.

**Not covered:** a dictionary whose values are all empty, such as
`Seq({0: ''}, length=4)`, means the same as `{}` but still builds a
`_PartiallyDefinedSequenceData` with only zero-length segments, as before:
its `repr` keeps the empty entry, `defined_ranges` is `((0, 0),)`, printing
it says "only partially defined", and `* 2` gives
`Seq({0: '', 4: ''}, length=8)`. Nothing crashes. This was so before the
fix, the upstream pull request does not touch it, and it is left as a
possible follow-up.

**Upstream status:** closed unmerged on 2026-10-12, withdrawn by its author
pending the outcome of upstream's draft policy on AI-assisted contributions
([#5241](https://github.com/biopython/biopython/pull/5241), still open on
2026-11-09). Michiel de Hoon approved it on 2026-10-10; the withdrawal is
not a technical objection. The licence box is ticked. Issue #5329 is still
open and upstream master has no fix.

**Related, not taken here:**
[#5330](https://github.com/biopython/biopython/pull/5330) by Aniruddha Adak,
the reporter, makes the same change, but its description dropped the
template's checklist, so there is no licence box ticked and it is not
eligible.

### [#5278](https://github.com/biopython/biopython/pull/5278) — Xinxiang Wang (yiyabo)

Fixes upstream issue
[#1022](https://github.com/biopython/biopython/issues/1022), reported in 2016
by Tyghe Vallard. `CharBuffer.next_word` treats `-` as punctuation, so an
unquoted name such as `dog-cat` in a NEXUS `TRANSLATE` table came out as three
words, `dog`, `-` and `cat`, and `Nexus._translate`, which expects a comma
after each name, raised `NexusError: Missing ','`. `Bio.Phylo`'s `"nexus"`
reader goes through the same code, so it could not read the trees in such a
file. The hyphen cannot simply leave `PUNCTUATION`, because `_parse_list`
needs it for ranges such as `4-8` in a `CHARSET` or `TAXSET`, as etal pointed
out on the issue in 2017. The pull request gives `next_word` an optional
`extra_chars` argument, punctuation to keep inside a word, defaulting to the
empty string, and passes `"-"` only when `_translate` reads a name. Every
other caller is unchanged.

**The single commit taken** (`264b17f15`, cherry-picked with `-x`), with
nothing left behind. It applied without conflicts, and has no `NEWS.rst` or
`CONTRIB.rst` hunk. The author left the box for being listed unticked
("optional, not wished"), so is not added to `CONTRIB.rst`.

**Correction on top**, in the fork's own commit. black reformats both files.
The comment that trailed the changed `_translate` line now sits on its own
line, where black would have split the call over three lines to make room
for it. The docstring, that comment and the test's comment cite the issue as
`biopython/biopython#1022`, because the fork's issue numbers are not
upstream's.

**Checked.** The new test, `test_translate_dashed_taxon_names`, fails on
main with the `NexusError` and passes here. A fork commit extends it to read
the same file with `Phylo.read(..., "nexus")`, the path the issue reported,
which also raises the `NexusError` on main. The table also reads with a
hyphen at either end of a name or doubled (`-lead`, `trail-`, `a--b`), and
with no space after a comma. Plain and quoted names read as before, and
`charset r = 2-4;` still gives `[1, 2, 3]`. `test_Nexus`, `test_Phylo`,
`test_Phylo_NeXML`, `test_PhyloXML`, `test_Phylo_CDAO`, `test_AlignIO`,
`test_Align_nexus`, `test_SeqIO` and the `Bio.Nexus.Nexus` doctests pass.

**Not fixed: an unquoted name with a hyphen in a `MATRIX` row.** `_matrix`
reads the taxon name with the same `next_word`, so `dog-cat ACGT` is still
taxon `dog` with the data `-catACGT`, here as in Biopython 1.88. With
`nchar=4` reading fails with
`Matrix Nchar 4 does not match data length (8) for taxon dog`; with `nchar=8`
it succeeds silently and gives the wrong taxon and data. That is outside #1022 and the pull request, and would be its
own change. Also unchanged: through `Bio.Phylo` a translated name with a
hyphen comes back quoted, as `'dog-cat'`, because `NexusIO` passes each name
through `Nexus.safename`; a quoted name with a hyphen always did.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with one commit and no activity since 2026-08-06. peterjc asked
whether it was AI-generated and what the author thought of upstream's draft
policy on AI-assisted contributions
([#5241](https://github.com/biopython/biopython/pull/5241), still open). The
author said he had used an AI assistant for the code and for the English of
the description, explained the change, and offered to close the pull request
if the project does not want such patches. There is no technical objection.
The licence box is ticked. Issue #1022 is still open and upstream master has
no fix. [#5397](https://github.com/biopython/biopython/pull/5397) by
TheNetherWatcher, also open, rewrites `CharBuffer` so that reading is linear
rather than quadratic; it touches the same `next_word` but is a different
change, not taken here, and its commit merges onto this one without conflict
(`git merge-tree`).

### [#5279](https://github.com/biopython/biopython/pull/5279) — Yiyabo

Answers upstream issue
[#5257](https://github.com/biopython/biopython/issues/5257), where Anika
Dureja asked for a radius of gyration function and posted her own. Adds
`Entity.radius_of_gyration(geometric=False)` to Structure, Model, Chain and
Residue: the root-mean-square distance of the atoms from their center of
mass, each weighted by `Atom.mass`, or equally with `geometric=True`. The
walk over the hierarchy that `center_of_mass()` did inline, which unpacks
disordered residues and atoms so that every alternative location counts in
full, moves into a private `_get_atoms()` that both methods share. The pull
request's expected values for 1LCD (13.944 weighted by mass and 13.729
geometric, and 11.279, 12.338 and 12.508 for chains A, B and C) were
reproduced with an independent float64 computation, and `center_of_mass()`
through `_get_atoms()` gives the same array as before.

**The single commit taken** (`5462b2c78`, cherry-picked with `-x`), which
applied cleanly. Nothing was left behind.

**Corrections on top**, in four fork commits:

- The empty-structure test opened `"PDB/disordered.pdb"` relative to the
  working directory, which raises `FileNotFoundError` when pytest starts
  outside `Tests/`, as CI does; it now uses `support.DATA` like the rest of
  the file. The docstring example read `PDB/1GBT.cif` with `PDBParser`,
  which raises `ValueError` on an mmCIF file, and gave 10.34; it was marked
  `+SKIP`, so nothing noticed. It now uses `MMCIFParser` and runs as a
  doctest, giving 16.27, which an independent float64 computation over all
  1761 atoms agrees with.
- `_get_atoms()` raised `IndexError: pop from an empty deque` for an entity
  with children but no atoms, such as a Structure whose only Model is empty
  or a Chain whose one Residue has no atoms, so neither method gave the
  `ValueError` the pull request's description promises for an empty entity.
  Copilot said so in review and nothing came of it. `center_of_mass()`
  crashes the same way in Biopython 1.88 and in upstream master.
  `_get_atoms()` now raises `ValueError` when its queue runs empty. The new
  test covers both methods on both entities and fails with the `IndexError`
  on main and on the adopted commit alone.
- The docstrings of both methods now say what was only in the code: every
  alternative location of a disordered atom or residue counts at full mass,
  not only the selected one that `get_atoms()` gives (238 atoms against 203
  in model 0 of `3JQH.cif`, which moves Rg by under 0.03 Angstrom), and an
  atom of unknown element has a mass of NaN, which makes the mass-weighted
  result NaN. The two empty-entity tests check the error message, so they
  tell "no children" from "no atoms".
- `radius_of_gyration()` summed in float32, so the result changed when the
  entity was moved: 1LCD translated by 100,000 Angstrom along each axis
  gave 14.000 instead of 13.944, and 30 copies of its models (33,840 atoms)
  translated by 10,000 gave 14.001. It now works in float64, which matches
  the independent computation to six places in both cases. A new test
  checks that translation leaves the result unchanged. `center_of_mass()`
  still averages in float32, as it does upstream.

Copilot's other comment, that `geometric=True` could pass `weights=None`
rather than an array of ones, is not taken: it is cosmetic.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with no activity since 2026-08-06, when the author wrote that an AI
assistant helped with the patch and its description. Issue #5257 is open.
The licence box is ticked. The author left the box for being listed
unticked and wrote "not wished", so is not added to `CONTRIB.rst`.

### [#3022](https://github.com/biopython/biopython/pull/3022) — Tony Tung

`Bio.Phylo.NewickIO.Writer.to_strings` built each tree's text with
`newickize`, a closure that called itself once per level of nesting, two
Python frames a level with its generator expression. So `Phylo.write` with
format `"newick"` raised `RecursionError` from depth 496 of a caterpillar, in
Biopython 1.88, upstream master and the fork alike, though `Phylo.read` reads
trees of any depth. The pull request replaces the recursion with a loop over a
deque of clades and two kinds of marker, one that closes an internal clade and
writes its label and branch information, and one that writes a comma, all into
a `StringIO`. It was written after Nextstrain's augur hit the limit
([nextstrain/augur#328](https://github.com/nextstrain/augur/issues/328)).

**The single commit taken** (`9edae6a89`, cherry-picked with `-x`), with
nothing left behind. It conflicted in both files it touches, on lines that
changed after 2020. `NewickIO.py` conflicted in three places. In the imports,
the pull request's `deque` joined the fork's `warnings` and `weakref`. On the
label-quoting line the fork's line was kept: the pull request backslash-escaped
a single quote, and the fork doubles it, as upstream master has since
`35b09e7c9`. The old recursive return line, which the fork had made an
f-string, goes, as the pull request removes it. In `CONTRIB.rst`, the author's
line went into its alphabetical place after Tommy Carstensen. Apart from the
quoting line, the code is as the author wrote it.

**Corrections on top**, in the fork's own commits. The two `.format` calls the
pull request added are f-strings, as the fork's ruff (UP032) requires. `tung`
joins the words codespell accepts, as upstream did for other names, because
codespell reads the surname as a misspelling of "tongue". A comment in
`BaseTree` that listed the Newick writer among the recursive traversals no
longer does; the #4634 entry below still gives it as a reason, which was true
when that was written. The markers keep their names, though `TerminalMarker`
marks the end of a non-terminal clade: its docstring says so, and a rename
would only make the code differ from the pull request.

**Checked.** The pull request has no test. `test_Phylo` gains one that writes
a 5000-deep caterpillar, with the default options and with `plain=True`, and
compares the text with the input; against main it fails with `RecursionError`,
and the rest of `test_Phylo` passes. The new writer's output was compared with
the recursive writer's from main over 109 tree sources (hand-made cases and
the Newick files under `Tests/Nexus` and `Tests/TreeConstruction`, each read
three ways, and the NEXUS, PhyloXML and NeXML test files), each as read, with
its confidences doubled, and with every internal clade given a name and a
confidence, so that the fork's confidence comments for named clades are
written, times 384 combinations of the writer's options: 125,568 comparisons,
each with the same text, or the same exception, and the same warnings, which
432 of them give. Depths 5000 and 100000 round-trip.

**Two other differences**, outside those comparisons. The new writer is not
faster: run in one process on the same trees, best of five, it took 1.19 times
as long as the recursive writer for a balanced tree of 131,072 tips and 1.28
times as long for a star of 200,000 tips. And a clade that is its own
descendant, which is not a tree, used to raise `RecursionError` at once; the
new writer loops until memory runs out (`MemoryError` after 23 s with the
address space capped at 1.5 GB).

**Not fixed: the `ladderize` option and NEXUS output, for deep trees.** With
`ladderize="left"` or `"right"`, the writer first calls `Tree.ladderize`,
which calls itself for each subclade and counts tips with `count_terminals`,
another recursive traversal. So Newick output with that option fails from depth
986. `NexusIO.write` uses this writer, then lists the tips with
`get_terminals`, whose `_preorder_traverse` is a recursive generator. So format
`"nexus"` fails from depth 992. Both failed from depth 496 or 497 before.
Making those traversals iterative is a separate change.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with one commit. The only review comment, peterjc's docstring nit, was
fixed by the author. On 2020-07-18 peterjc asked etal to review it, and nobody
did; the author's last ping was on 2020-10-05. The licence box is ticked, and
on 2020-07-07 the author said his employer's legal team had approved ticking
it. He added himself to `CONTRIB.rst`.

### [#5350](https://github.com/biopython/biopython/pull/5350) — Brandon Mokrytzki

Fixes upstream issue
[#4878](https://github.com/biopython/biopython/issues/4878), and also
[#4879](https://github.com/biopython/biopython/issues/4879), whose half the
fork already had from #4918 below. The parser keeps a strand prefix with the
molecule type, as in `"ss-DNA"`, and the GenBank writer wrote the whole
string from column 48 of the LOCUS line, three columns right of the 45-47
that section 3.4.4.2 of the GenBank release notes gives the strand.

**Only the two test commits taken** (`6037becc9` and `cfeed7071`,
cherry-picked with `-x`). `6037becc9` conflicted on context only, because in
the pull request it follows the NA test the fork does not take; the
resolution keeps just the new test. `cfeed7071` applied cleanly.

**Left behind, and why.**

- `c43a3816d`, the writer fix. It pads a molecule type with no prefix with
  three spaces and still calls `ljust(7)`, which leaves the type four
  columns. `snRNA` and `ncRNA` then give an 80-character LOCUS line, `snoRNA`
  81 and `circRNA` 82, with no warning, and `Bio.SeqIO` cannot read them back
  (`LOCUS line does not contain space at position 68`). The pull request's
  test table stops at four-letter types, so its own tests pass. The fork's
  fix instead splits the prefix off into its own three columns and keeps the
  type seven wide.
- `0cdee08e6`, and the `__init__.py` and 80-column hunks of `3dacd8898`: the
  NA half, which the fork has in a stricter form from #4918 (an exact token
  match, which also reads `ss-NA`). The one hunk of `3dacd8898` the fork
  lacked, accepting `NA` in the writer's check of a LOCUS line too long for
  the usual columns, is ported by hand in the fork's commit. It tests the
  molecule type after any prefix, not the whole token, so `ss-NA` passes too.
- `fcd6451f1` and `969b9161a`, the `example_NA.gb` fixture and a test that
  writes into `Tests/GenBank/`. The #4918 tests already read NA, and the new
  long-name test below covers the write path that was missing. peterjc also
  found the fixture odd.
- `70303b3f9`, the author's `CONTRIB.rst` and `NEWS.rst` lines; the author is
  added to the fork's `CONTRIB.rst` instead.

**Corrections on top**, in the fork's own commit: the writer fix above. The
test table is rewritten as padded strings with a PEP257 summary, as peterjc's
review asked, widened to all ten columns, extended with `snRNA`, `ncRNA`,
`snoRNA` and `circRNA` with and without a prefix, and each line is read back.
Three tests are added. One writes NCBI's `NC_001422` (`ss-DNA`) and a
`ds-DNA` test plasmid back with their LOCUS lines unchanged. One writes `NA`,
`ss-NA`, `DNA` and `ss-DNA` with a 27-character locus name and checks the
layout after the length. One writes `protein`, `ss-protein` and `ds-PROTEIN`:
no parser gives a strand with a protein, and the writer drops it with the
type, where splitting it off first would otherwise leave `ss-` beside the
units `aa`, a line `Bio.SeqIO` cannot read back. Against main before the fix,
35 of the 49 subtests fail. The 14 that pass have no prefix, and their output
does not change: the 12 table rows, `DNA` with the long name, and `protein`.
With `c43a3816d` in place of the fork's fix, the eight rows of five to seven
letters fail, and so do `NA` and `ss-NA` with the long name and the two
proteins with a prefix.

**Upstream status:** as of 2026-11-09 the pull request is open, with no
commits since 2026-10-02 and style-only review comments from peterjc. On #4878
the same day another contributor asked to be assigned the issue, and the
author replied that the review changes would come by 15 November and that the
other contributor was welcome to take it over if not. Nobody upstream has
mentioned the overflow. The licence box is ticked.

### [#5389](https://github.com/biopython/biopython/pull/5389) — Michiel de Hoon

Titled "PDB deprecation warnings", but one of its two lines fixes a silent
defect. `PDBMLParser`'s `_parse_resolution_from` tested each element it looked
up with `if element:`. An `Element` is true only if it has children, and
`ls_d_res_high`, `d_res_high` and `resolution` are leaves, so the test failed
for every one: `header["resolution"]` was `None` for every PDBML file, in
Biopython 1.88 and in upstream master, with a `DeprecationWarning` for each
lookup on Python 3.12 and later. The fix is `is not None`. 1GBT, 3JQH and
6WG6 now give 2.0, 2.201 and 3.54, as `MMCIFParser` gives from their mmCIF
files.

The other line replaces the polar angle in `vectors.multi_coord_space`,
`arccos(z / r)`, with `arctan2(hypot(x, y), z)`. That supersedes the fork's
own fix from #5127 below, which gave `np.divide` an `out=` array so that a
hedron whose third atom sits on its second gets a polar angle of 0 rather than
whatever `np.empty` held. Upstream closed #5127 unmerged on 2026-10-23, 33
minutes before the same author opened this pull request. The `arctan2` form
was taken over the fork's line because it is the better formula and drops the
`where=`/`out=` masking: `arccos` is ill-conditioned near 1, and for
`p = (1e-8, 0, 1.5)` it returned 0 where the angle is 6.7e-9. Elsewhere the
two differ by rounding error, at most 4.4e-13 radians on a million random
vectors at each of nine scales from 1e-3 to 1e4. That still shows in the 13
decimals `write_PIC` prints: on nine test structures (1A7G, 1A8O, 1GBT, 1LCD,
2XHE, 3JQH, 4CUP, 4ZHL, 6WG6) every PIC file changes, by at most 2.3e-12
degrees, while the coordinates rebuilt from them are bit-identical.

**The single commit taken** (`6bf45f584`). `PDBMLParser.py` merged cleanly;
`vectors.py` conflicted on the line #5127 had changed here, and the
conflict was resolved by taking this pull request's line. Nothing was left
behind.

Corrections on top, in separate commits. The pull request has no test, so
`test_PDB_PDBMLParser` gains two. One asserts the resolutions of 1GBT and
3JQH; against main it fails with `None != 2.0` and `None != 2.201`. All three
fixtures are X-ray and stop at the first place looked, `refine`, so the other
edits 3JQH: without `refine`'s value, or with it nil, the resolution comes
from `refine_hist`; without that too, from an added `em_3d_reconstruction`,
where cryo-EM entries record it; and with none, it is `None`. And the
`arctan2` form changes one case the `arccos` form got right: for a
zero-length vector, `arctan2(0, -0.0)` is pi, not 0, and a z of -0.0, as
`float("-0.000")` gives, minus a z of 0.0 is -0.0. The rotation for such a
hedron then maps Z onto -Z, where `get_spherical_coordinates` and the #5127
`NEWS.rst` entry say the polar angle is 0. Adding `0.0` to z turns -0.0 into
+0.0; the zero-radius test gains that hedron, which fails with the bare
upstream line and passes on main and here.

**Upstream status:** as of 2026-11-09 the pull request is open and not a
draft, with no reviews and only a codecov comment; its CI passes. The licence
box is ticked. The author is already in `CONTRIB.rst`.

### [#5365](https://github.com/biopython/biopython/pull/5365) — Angus Barlow

Fixes upstream issue
[#5332](https://github.com/biopython/biopython/issues/5332).
`DistanceCalculator` with the `identity` model counted matches only over the
columns where neither letter is in `skip_letters`, but divided by
`len(seq1)`, every column. So with an explicit `skip_letters` each skipped
column made the distance larger: in an alignment of the rows `AC--GT`,
`AC--GT` and `ACTTGT`, every pair was 0.333 apart rather than 0, with either
`MultipleSeqAlignment` or `Bio.Align.Alignment`, and `AC--GT` against `AC--GA`
gave 0.5 where `ACGT` against `ACGA` gives 0.25. An `Alignment` was spared
only the columns that are gaps in every row, which `Align.read` drops. The
scoring-matrix branch already left skipped columns out of its total. The fix
counts the total over the same columns as the matches. `skip_letters` came in
with commit `f25c05203` in 2017, for upstream issue
[#491](https://github.com/biopython/biopython/issues/491), so that `identity`
could skip gaps as the other models do; that commit rewrote the match count
and left the total alone, so the mismatch is an oversight, not a choice. The
default `identity` model skips nothing and gives the same numbers as before.

**Both commits taken**: `4e894b440b` (the fix, its test, and the author's
line in `CONTRIB.rst`) and `69921cb009` (that line moved into alphabetical
order). The adoption of #5363 had already put the same line in the same place,
so together the two leave `CONTRIB.rst` unchanged. The first conflicted only
in `NEWS.rst`. Its hunk was a contributor list for upstream's in-progress
Biopython 1.89 section, which in this tree is preserved upstream history, so
the conflict was resolved by leaving that section as it was, as for #5391; the
fork's own section has a new paragraph instead. That hunk is all that was left
behind. There was no review discussion. The commit message credits Devin, an
AI tool.

One correction sits on top. The pull request's test reads its alignments only
with `AlignIO`, so it never passes a `Bio.Align.Alignment`, which was affected
too: `AC--GT` against `ACTTGT` gave 0.333 with either kind of alignment. Its
identical gapped rows would not help there, because `Align.read` drops the
columns that are gaps in every row, which is why #5401's `Alignment` test
passed without the fix (below). A new test checks `AC--GT` and `AC--GA`
against `ACTTGT`, each in both kinds of alignment. With only
`TreeConstruction.py` reverted, it and the pull request's test fail, and
nothing else in `test_TreeConstruction.py` or `test_Consensus.py` does. The
comment above the identity loop, which said no letters were skipped and had
been wrong since 2017, now says what the loop does, and the test cites the
issue as `biopython/biopython#5332`. One side effect, checked: when every
column is skipped, the identity branch now returns the int `1` from the
existing guard, where it returned `1.0`. The scoring-matrix branch already
did, and `DistanceMatrix` prints and writes both the same.

As of 2026-11-09 the pull request is open, with no comments or reviews, and
the licence box is ticked. Issue #5332 is open and upstream master has no
fix.

**Related, not taken:**
[#5401](https://github.com/biopython/biopython/pull/5401) by simplyshree is
the same fix, opened and closed unmerged by its author on 2026-11-03, after
peterjc asked for the licence box to be ticked and whether AI was used. Its
`Alignment` test used the identical gapped rows, which `Align.read` reduces
to ungapped ones, so it passed before the fix as well.
[#5333](https://github.com/biopython/biopython/pull/5333) by the issue's
author, aniruddhaadak80, was closed on 2026-10-03 without the licence box,
and is not eligible.

### [#5363](https://github.com/biopython/biopython/pull/5363) — Angus Barlow

Fixes upstream issue
[#5340](https://github.com/biopython/biopython/issues/5340), reported by
Ziqi (Adam) Chen. `ProteinAnalysis.protein_scale` raised
`ZeroDivisionError` for a window of 1, because `_weight_list` divides by
`window - 1`. For an even window it counted `subsequence[window // 2]`
twice, once in the innermost mirrored pair and again as the "middle", so
the later of the two middle residues carried an extra weight of 1, and the
profile was not the weighted mean of the window and depended on the
sequence's direction. The code dates from 2003 (`19f997e2f`) and was
written for odd windows only; nothing in its history deals with even
windows. The reporter's case
is a real one: Hopp and Woods' scale, `ProtParamData.hw`, is averaged over
six residues, and the test suite's 152-residue protein came out up to 0.56
away from that mean. The fix gives a window of 1 a weight unit of 0, skips
the middle term for even windows, and divides by the two-sided weight sum.
Checked here against a weighted mean written independently, a linear ramp
from `edge` at both ends to 1 at the window centre: windows 2, 4, 6 and 10
agree to 1e-15, a reversed sequence gives exactly the reversed profile, and
windows 3, 5, 9 and 21 give the same floats as before, bit for bit.

**One of the two commits taken** (`7f4d775969`, the fix and its two
regression tests), which applied with no conflicts. Left behind:
`5fe8afa5e7`, which adds the author to `CONTRIB.rst` and to a "Many thanks"
contributor list in upstream's in-progress 1.89 section of `NEWS.rst`. That
section is upstream history in this tree, so the fork writes its own `NEWS.rst`
paragraph instead; the `CONTRIB.rst` line is added by hand, as the author
ticked the box for adding his name.

Corrections on top. The upstream direction test passes for any symmetric
weighting, so a new test pins the values of even windows against the
weighted mean written out by hand, the Hopp and Woods six-residue mean
among them. The docstring now describes a window of 1 and even windows,
whose weights rise towards 1 at the point between the two middle residues
(0.4, 0.8, 0.8, 0.4 for a window of 4 with an edge of 0.4). And the fix
turned a window of 2 with an edge of 0, where both residues have a weight
of 0, from a wrong answer (the second residue's value) into a new
`ZeroDivisionError`; it now raises `ValueError` naming the window and edge.
A window smaller than 1, which raised `IndexError`, is refused up front with
its own `ValueError`. With these five new tests and only `ProtParam.py`
reverted, all five fail.

**Upstream status:** as of 2026-11-09 the pull request is open, with no
reviews or comments since it was opened on 2026-10-08, and its CI passes.
On the issue, peterjc first suggested documenting and limiting the window
size instead, accepted the Hopp and Woods case as a real use, and asked
the reporter whether #5363 worked for them; the reporter replied "LGTM".
The licence box is ticked. Issue #5340 is still open and upstream master
has no fix.

### [#5392](https://github.com/biopython/biopython/pull/5392) — Zijian Zeng

The `LWL85` method of `calculate_dn_ds` weighted the transitional and
transversional distances at twofold degenerate sites, `A[1]` and `B[1]`, by
the fourfold site count `L[2]` instead of the twofold count `L[1]`, so dS was
wrong whenever a pair differs at a twofold site and dN whenever one of those
differences is a transversion, unless the two counts happen to be equal. The
first version of the method, written in 2013 for what was then
`Bio.CodonAlign` (`61946aaff`), had the right weights; a tidying commit
fourteen minutes later (`8ed3c799f`) moved the counts into lists and wrote
`L[2]` where the first version had the twofold count. The fix swaps `L[2]`
for `L[1]` in those two terms. Checked here against Li, Wu and Luo's formula
worked by hand for the pull request's three fixtures, against the LWL85 sums
in PAML's `yn00.c`, which use `L[1] * A[1]` and `L[1] * B[1]` with the same
site-class order, and by running `yn00` built from PAML commit `4c7902fe` on
the codons Biopython passes to `_lwl85`: it agrees to within 6e-16 on the
three fixtures, on one with three times as many fourfold as twofold sites
(where the old dS was three times too high), and on the Tutorial pair.

**Not fixed here:** on the test suite's EGFR pair `yn00` does not agree (dN
0.020491, dS 0.019260, against 0.021132 and 0.016980 here), for an older
reason this pull request does not touch. `_diff_codon` classes each
differing position by its degeneracy in the first sequence's codon only,
where `yn00` counts half in each codon's class. So LWL85 depends on the order
of the two sequences (EGFR reversed gives dN 0.019854, dS 0.021552), as it
did before the fix, while NG86 does not. The site counts agree with `yn00`;
averaging the forward and reversed difference counts gives `yn00`'s counts
exactly, and the corrected formula on those reproduces its dN and dS to all
15 printed digits, so the weight fix itself is right. No upstream issue
covers this (searched on 2026-11-09).

**The single commit taken** (`d30e7e844`), with the analysis fix, its new
test module `test_Align_analysis_lwl85.py`, the updated LWL85 expectations in
`test_Align_codonalign.py` and `test_codonalign.py`, and the Tutorial's
E. coli/Shigella dS, which goes from 0.183 to 0.225 on top of the block this
fork rewrote in #270. With only `analysis.py` reverted, 11 of these tests
fail, including all 8 new subtests but not the equal-counts control, and the
Tutorial doctest fails on that line.

Left behind: the `Bio/codonalign/codonseq.py` hunk. Upstream still has a
second copy of `_lwl85` there; in this fork `codonseq.py` imports it from
`Bio.Align.analysis` (#154), so the hunk had nothing to apply to and the
conflict was resolved by keeping the fork's file. The new tests go through
`cal_dn_ds` as well as `calculate_dn_ds`, so both APIs are covered.

One correction on top: the old draft chapter `Doc/doc.rst`, which Sphinx
excludes and nothing doctests, showed the old LWL85 values for the EGFR
example, and they still reproduced exactly, so it now shows the corrected
ones.

As of 2026-11-09 the pull request is open, with no comments or reviews. The
author ticked the box for not wishing to be added to the contributor
listings, so is not added to `CONTRIB.rst`, and declared substantial AI
assistance in the description.

### [#5356](https://github.com/biopython/biopython/pull/5356) — Zehao Wang

Fixes upstream issue
[#5312](https://github.com/biopython/biopython/issues/5312), reported by
Ziqi (Adam) Chen. `IsoelectricPoint.pi()` bisected inside a fixed interval
of pH 4.05 to 12, so a sequence whose pI lies outside it got the end of the
interval, where the net charge is far from zero: ten aspartates gave 4.05
with a charge of -4.24, ten arginines 12.00 with +4.00, and `DD` and `ED`
4.05. `ProtParam.ProteinAnalysis.isoelectric_point()` goes through the same
code. The fix widens the interval in steps of 1 pH unit, not below 0 or
above 14, until the charge changes sign across it, then bisects as before
in a new private `_pi_bisect`. Measured here on 5000 random sequences of 1
to 200 residues: 95 changed, all of them old answers at the end of the
interval, and the other 4905 give the same float as before.

**One of the four commits taken** (`159aa6795a`, the `pi()` fix and its
tests), with the one-line conflict in the import block of
`Tests/test_SeqUtils.py` resolved by keeping both lines. The other three
were carried over from the author's branch for the `QCPSuperimposer` fix in
[#5351](https://github.com/biopython/biopython/pull/5351) and are not part
of this fix: `4c1273ece9` (an `eigh` fallback in `Bio/PDB/qcprot.py`),
`a3d9759eb6` (a collinear `QCPSuperimposer` test) and `6f4c68839e` (adding
the author to `CONTRIB.rst`). #5351 is not waiting on a decision here: this
fork had already fixed the same `QCPSuperimposer` defect, upstream issue
[#5346](https://github.com/biopython/biopython/issues/5346), in its own
commit `0ff7b558c`, which uses the author's diagnosis and key-matrix check
from #5351 and does not take the `eigh` fallback. The author ticked the
box for being listed already, through #5351, which was closed unmerged,
so this fork adds the line from `6f4c68839e` to `CONTRIB.rst` by hand.

One correction sits on top. The `pi()` docstring no longer calls it
recursive and says that the search stops at pH 0 and 14, and two tests pin
the changes users can see beyond the bug itself: an explicit `min_` and
`max_` that miss the pI are widened too, so
`IsoelectricPoint("INGAR").pi(5.5, 5, 6)` gave 6.0 and now gives 9.75; and
1000 arginines, whose pI under this model is above 14, get 14.

**Upstream status:** closed unmerged on 2026-10-12, withdrawn by its author
pending the outcome of upstream's draft policy on AI-assisted contributions
([#5241](https://github.com/biopython/biopython/pull/5241), still open on
2026-11-09). The withdrawal says the fix is sound; it is not a technical
objection. The licence box is ticked. Issue #5312 is still open and
upstream master has no fix.

**Related, not taken here:**
[#5324](https://github.com/biopython/biopython/pull/5324) by rollroyces
(licence box ticked) fixes the same issue by changing the default interval
to 0 to 14. Its author closed it on 2026-09-19 at a maintainer's request,
pending the same policy. It is not taken because it moves almost every pI:
measured on the same 5000 sequences, 4902 of the 4905 whose pI lies inside
the old interval would change, by up to 0.00008, since the bisection visits
different midpoints. #5356 leaves them all as they were.

### [#5405](https://github.com/biopython/biopython/pull/5405) — Muneer Alam

Fixes upstream issue
[#5368](https://github.com/biopython/biopython/issues/5368), reported by
Ziqi (Adam) Chen. For point sets that lie exactly on a line, such as two
points or points along a helix axis, `QCPSuperimposer.get_rms()` was wrong
while the rotation was right: the issue's three identical points on a line
60 Å long reported 2.816 Å and printed "Newton-Rhapson did not converge
after 50 iterations". The largest root of QCP's quartic is then double and
the starting guess E0 sits on it, so f and f' are both rounding noise and
`f / (f_prime + 1e-11)` throws the iterate to about 4e8. When the iteration
does not converge, the pull request takes the largest eigenvalue of the 4×4
key matrix from `numpy.linalg.eigvalsh`, which also removes the print, and
where the cofactor cascade falls back to the Kabsch SVD it measures the RMSD
on the fit.

**The single commit taken** (`35beac4fe`). It conflicted only because this
tree had already moved that Kabsch block into `_kabsch_rotation`, for the
check of the quaternion against the key matrix
([#5346](https://github.com/biopython/biopython/issues/5346)); the
resolution keeps the same code. He is added to `CONTRIB.rst`, as Abdel ATIA
was for #5175. Both left the optional contributor-listing box unticked, but
that box covers "added my name, am listed already, or do not wish to be
listed", so a blank box declines nothing.

That check is what makes the pull request safe here. As upstream has it,
without the check, the eigenvalue alone leaves the rotations of many wide
collinear sets wrong, and with the print gone nothing says so: 144 of 300
pairs of points spread ±150 Å got a wrong rotation, by up to 286 Å, none
with a message. In this tree, none did.

One correction sits on top. The check's own Kabsch fallback still reported
the RMSD from the eigenvalue, and when Newton-Raphson does converge at a
double root it finds the root only to about the square root of machine
precision, so the RMSD reported was still silently out by up to 0.025 Å in
22 to 117 of each 300 collinear sets. The RMSD is now measured on the fit
on every path, so `get_rms()` always agrees with `get_transformed()`; for
other fits the rotation is unchanged and the RMSD moves only at the level of
rounding, except for points spread over only a few thousandths of an Å or
less. There Newton-Raphson fails too, because it divides by f' + 1e-11 and
f' scales as the sixth power of the coordinates: it runs out of iterations,
or stops at once on its starting guess E0. The RMSD from the eigenvalue was
too small in 348 of 500 random fits spread ±0.001 Å, often by most of its
value, and most of those at ±0.0001 Å had no message; it is now right, with
the rotation bit-identical. A new test fits 60 wide collinear sets and fails on the previous
main, on this commit alone, and on the commit as upstream has it.

**Not fixed: points close to a line but not on it.** These can still get a
rotation slightly worse than the best one, with no warning. In 1000 exact
rotated copies each of 3 to 29 points lying 0.001 or 0.01 Å off a line
spread ±100 or ±150 Å, 34 to 189 fits were worse than SVD's by more than
1e-3 Å, up to 0.017 Å; the previous main had 35 to 246, with the same worst
case. All of them take the quaternion path, so the quaternion's error stays
within the tolerance of the key-matrix check. `get_rms()` now reports the
RMSD of that fit: it is off the best RMSD by more than 1e-3 Å in those same
34 to 189 sets, where the previous main's value from the eigenvalue was off
in 615 to 746. A report upstream that #5368 is fixed should say this.

One behaviour change beyond the issue: coordinates containing NaN or
infinity now raise `numpy.linalg.LinAlgError` from `eigvalsh`, as
`SVDSuperimposer` raises it, where they gave a NaN RMSD and the print.

On 2026-11-06 upstream closed issue #5368 as a duplicate, and the pull
request with "This has been solved in #5317." It is not: #5317 is the Kabsch
fallback for vanishing cofactors, which fixes rotations, and on 2026-11-09
upstream master, which has it, still printed the message and reported
2.816 Å for the issue's script.

### [#5373](https://github.com/biopython/biopython/pull/5373) — Mohit Arvind Khakharia

The earlier fix for the same issue, #5368. **Only one test method taken**,
`test_rmsd_identical_two_points_rotated`, which fits two points to a rotated
copy: the previous main reported 2.5e-4 Å there, with no message, and
#5405's tests have no such case. The commit (`2af0dbfee`) is cherry-picked
so that he stays its author, and its message, unchanged as for every
adoption, still describes the code change that was not taken.

Left behind:

- The `qcprot.py` hunk, which stops Newton-Raphson as soon as f or f' is not
  positive. That holds at the root in exact arithmetic, but at E0 both are
  rounding noise of either sign: when both are positive the step throws the
  iterate below the root and the new check then stops it there. Applied to
  the previous main, the RMSD reported was wrong by more than 1e-3 Å in 52
  to 164 of each 300 collinear sets, by up to 352 Å, all without a message,
  where main printed one for its large errors.
- `test_rmsd_identical_collinear`, the same three points and assertions as
  the first case of #5405's test, at a looser tolerance.
- The `NEWS.rst` paragraph, written for upstream's own in-progress section;
  the fork's own paragraph replaces it.

Its author closed it, still a draft, on 2026-11-05 as a duplicate of #5317,
which it is not (see #5405 above). peterjc's question whether it was written
with AI was not answered. He is already in `CONTRIB.rst`.

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
tests the upstream pull request was asked for but never grew. That correction
missed the writer's other check, for a LOCUS line too long for the usual
columns, so `NA` with a long locus name still raised `ValueError`; the #5350
adoption above fixes it.

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

Upstream closed this pull request unmerged on 2026-10-23, and 33 minutes
later the same author opened #5389, which changes the same line another way.
Its `out=` line has since been replaced here by #5389's `arctan2` form; see
that entry above.

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
