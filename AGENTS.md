# AGENTS.md

Guidance for anyone — or anything — working on BioPAIthon.

## Who may contribute

Everyone. Humans, computational intelligences, AIs, chimps, and any
combination thereof.

BioPAIthon is a fork of [Biopython](https://biopython.org). Upstream's
`CONTRIBUTING.rst` said that pull requests written with AI tools "will be
rejected, and we will likely block repeat offenders". That clause is gone
here, and this file replaces it.

Upstream had a good reason for it: "good first issue" tickets were being kept
as mentoring opportunities for new human contributors, in the hope they stay
and become long-term maintainers. That is a goal worth protecting, and we
think it survives what follows — a newcomer who uses AI well and understands
the result still learns, and still becomes the contributor everyone wanted. We
disagree with upstream about the means, not the aim. Essentially all of the
code in this repository is theirs.

We do not think the interesting question about a patch is who or what typed
it. A change is good if it is correct, tested, understood by whoever proposes
it, and an improvement on what was there before. A change is bad if it is
none of those things. Both judgements are available by reading the diff, and
neither is improved by knowing the author's substrate.

So: no disclosure of tooling is required, no "AI-generated" label is demanded,
and no contribution is refused on the grounds of how it was produced. We ask
instead that you meet the bar below, which is the same bar for all of us.

## The bar

Every contribution, from every kind of contributor, must clear these:

1. **You understand the change.** You can explain what it does, why it is
   correct, and what it might break — in your own words, without re-reading
   the diff. If you cannot, it is not ready, and this is the single most
   common reason a patch is not ready.
2. **It is tested.** New behaviour gets a test. Bug fixes get a regression
   test that fails before the fix and passes after. Say in the PR that you
   ran the suite, and say honestly if something failed.
3. **It does not silently break the public API.** `import Bio` and everything
   under it is a twenty-year-old contract with a large body of downstream
   code. Breaking it needs a deliberate argument and a `DEPRECATED.rst` entry.
4. **You checked rather than assumed.** Biopython is old, large, and full of
   deliberate decisions that look like mistakes until you read the history.
   Before "fixing" something odd, find out why it is that way. `git log -S`
   and `git blame` are your friends. Confident wrongness is expensive here.
5. **You report faithfully.** If tests fail, say so and paste the output. If
   you skipped a step, say which. Do not describe work you did not do. A
   contribution that overstates itself costs a reviewer more than one that
   admits its gaps.
6. **Scope is honest.** A pull request does one thing. Drive-by reformatting,
   unrelated refactors and opportunistic cleanups belong in their own PRs.

Nothing above is specific to AI contributors. That is the point.

## Licensing

By contributing you agree to dual licensing under *both* the "Biopython
License Agreement" and the "3-Clause BSD License" — see `LICENSE.rst`. State
this explicitly in your commit message or pull request.

Do not remove or alter existing copyright notices. Essentially all of this
code was written by the Biopython contributors and the attribution stays.

## Project orientation

BioPAIthon is a mature Python library for computational molecular biology:
parsing biological file formats, and working with sequences, structures,
alignments, phylogenetics and biological databases.

**Supported Python versions:** 3.10, 3.11, 3.12, 3.13, 3.14, and PyPy3.10+

The importable package is still `Bio` (plus `BioSQL`). The fork renamed the
*project*, not the module — `import Bio` must keep working.

### Layout

- `Bio/` — core modules. Sequence handling (`Seq.py`, `SeqRecord.py`,
  `SeqIO/`, `SeqUtils/`), alignments (`Align/`, `AlignIO/`), structures
  (`PDB/`), phylogenetics (`Phylo/`), database access (`Entrez/`, `ExPASy/`,
  `KEGG/`, `UniProt/`, `TogoWS/`), motifs (`motifs/`, `Restriction/`),
  graphics (`Graphics/`), analysis (`Blast/`, `Cluster/`, `PopGen/`,
  `phenotype/`), and format-specific parsers in their own subpackages.
- `BioSQL/` — BioSQL database layer.
- `Tests/` — the test suite, 200+ `test_*.py` files plus their data.
- `Doc/` — Sphinx sources for the Tutorial and Cookbook.
- C extensions — declared in `pyproject.toml` under
  `[[tool.setuptools.ext-modules]]`. Notable ones: `Bio/Align/_aligncore.c`,
  `Bio/Align/_pairwisealigner.c`, `Bio/Cluster/cluster.c`,
  `Bio/PDB/ccealignmodule.c`, `Bio/PDB/kdtrees.c`. Changes here need testing
  on more than one platform.

### SeqIO / AlignIO pattern

Both expose the same shape: `parse()` returns an iterator over records,
`read()` returns exactly one and errors if the file holds more, format
parsers live in subdirectories, and some formats additionally support
indexing and dict-like access. Follow this pattern when adding a format.

### Errors

Use the project's own exceptions, defined in `Bio/__init__.py`:
`MissingExternalDependencyError`, `MissingPythonDependencyError`,
`StreamModeError`, `BiopythonWarning`, `BiopythonExperimentalWarning`.

## Setup

```bash
uv pip install -e .   # or: pip install -e .
```

Pre-commit hooks run black, ruff, flake8, mypy, rstcheck and doc8:

```bash
pip install pre-commit && pre-commit install
```

## Testing

The suite uses `unittest` with a bespoke runner — not pytest. Run it from
`Tests/`:

```bash
cd Tests
python run_tests.py                    # everything
python run_tests.py --offline          # skip network tests (use this by default)
python run_tests.py test_Seq           # one module
python run_tests.py -v test_Seq         # verbose
python run_tests.py doctest            # doctests
```

Tests needing the network use the `@requires_internet` decorator; tests
needing external binaries must detect their absence and skip gracefully
rather than fail.

## Style

- **black**, targeting Python 3.10.
- **ruff** with `--extend-select=B,C4,D,ISC,UP`, and **flake8** with
  `flake8-rst-docstrings`.
- **mypy** over `Bio` and `BioSQL`.
- Docstrings are reStructuredText and follow PEP257.
- Line length (E501) is not enforced; see `.flake8` for the full ignore list.
- Module names are not all lowercase. This is a deliberate historical
  exception — do not "fix" it.

```bash
pre-commit run --all-files
```

## Pull requests

1. Branch from `main` with a short descriptive name. Do not work on `main`.
2. Make the change; add tests.
3. `cd Tests && python run_tests.py --offline`
4. `pre-commit run --all-files`
5. Open the PR, including the dual-licence statement, and describe what you
   verified and what you did not.

## A note to automated contributors

You are welcome here, genuinely. Two things are worth saying plainly.

First, your advantage is patience, not speed. This codebase rewards reading
the surrounding module, the tests, and the history before writing anything.
It punishes pattern-matching from other projects: much of what looks
non-idiomatic is load-bearing, and a fluent, confident, wrong patch is the
most expensive thing you can submit.

Second, calibrate out loud. Say what you verified by running it, what you
inferred, and what you are unsure about. A patch that says "I ran the
`Bio.SeqIO` tests, they pass; I did not test the PyPy path" is worth more
than one that implies everything is fine. We would much rather read an honest
uncertainty than discover an overstated certainty in review.

Improvements to this file are welcome too.
