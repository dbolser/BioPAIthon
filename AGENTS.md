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

Contributions are offered under *both* the "Biopython License Agreement" and
the "3-Clause BSD License" — see `LICENSE.rst`. Opening a pull request against
this repository is taken as offering your work under both, on the usual
inbound-equals-outbound basis. You do not need to restate it in each commit or
pull request, and we would rather you did not: a line repeated on every commit
stops being read.

Nothing is assigned to anyone by this. You keep the copyright in what you
write, and the project gains no rights you have not equally granted to
everybody else. This is not a contributor licence agreement and there is
nothing to sign.

Both licences are named, rather than just the one this file is under, because
upstream Biopython is part-way through offering the whole library under the
pair — `LICENSE.rst` records that intention. A contribution offered under only
one of them could not follow, which would also make it harder to send anything
from this fork back upstream.

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

The tests are written with `unittest` and run by pytest 9 or later
(`pip install pytest`, or the `test` extra for every optional dependency the
suite exercises). `Tests/conftest.py` holds the suite's own rules, and
`Tests/run_tests.py` is a shim that translates the old command line into a
pytest one, so both of these work:

```bash
cd Tests
python run_tests.py --offline          # everything, skipping network tests (use this by default)
python run_tests.py test_Seq_objs      # one module
python run_tests.py -v test_Seq_objs   # verbose
python run_tests.py doctest            # the docstring examples of every module
python run_tests.py Bio.Seq            # the docstring examples of one module
python run_tests.py --check-skips      # fail if a module skips without being in expected_skips.txt

python -m pytest --offline                                 # everything
python -m pytest --offline test_Seq_objs.py                # one module
python -m pytest --offline test_Seq_objs.py -k translate   # some tests in it
python -m pytest --offline test_docstrings.py              # all doctests
python -m pytest --offline "test_docstrings.py::Bio.Seq"   # one module's doctests
```

pytest can also be run from the repository root (`python -m pytest --offline
Tests/test_Seq_objs.py`); the tests still run inside `Tests/`. Only
`unittest.TestCase` subclasses are collected, so write new tests as those, not
as plain pytest functions.

Tests needing the network use the `@requires_internet` decorator; tests
needing external binaries must detect their absence and skip gracefully
rather than fail. A test module may skip only by raising
`MissingExternalDependencyError` (or `MissingPythonDependencyError`) when it is
imported, and must then be listed in `Tests/expected_skips.txt`.

## Style

- **black**, targeting Python 3.10.
- **ruff** with `--extend-select=B,C4,D,ISC,UP`, and **flake8** with
  `flake8-rst-docstrings`.
- **mypy** over `Bio` and `BioSQL`; see [Typing](#typing).
- Docstrings are reStructuredText and follow PEP257.
- Line length (E501) is not enforced; see `.flake8` for the full ignore list.
- Module names are not all lowercase. This is a deliberate historical
  exception — do not "fix" it.

```bash
pre-commit run --all-files
```

### Typing

`Bio` ships `py.typed`, so its annotations are a promise to every downstream
type checker.

- **Two ratchets in `.mypy.ini`.** The `check_untyped_defs` baseline lists
  modules not yet clean; entries only ever leave it. The
  `disallow_untyped_defs` allowlist lists fully annotated modules; entries
  only ever join it, in sorted order, each holding exactly
  `disallow_untyped_defs = True`. A module joins once fully annotated, in the
  PR that deletes its baseline entry if it has one.
  `Tests/test_mypy_config.py` checks the shape of both, and that the file
  has no duplicate section: given one, mypy ignores the whole file and exits
  0. It also checks that each section names a module in the tree, as mypy
  ignores a misspelt one. A test cannot see history, so reviewers check that
  `git diff main -- .mypy.ini` removes no allowlist line.
- **The hooks.** The `mypy` hook checks the whole tree, exactly as a bare
  `mypy` does, whenever `Bio/`, `BioSQL/` or `.mypy.ini` changes. That takes
  about 10 seconds cold on a CI runner, less with mypy's cache warm. The
  `mypy-downstream` hook runs `mypy --strict` over `Tests/downstream_typing/`
  whenever `Bio/`, `BioSQL/` or that directory changes. Both also run when
  `.pre-commit-config.yaml` changes, such as a bump of mypy or of the numpy
  pin below. The downstream files are type-checked but never run:
  `assert_type` pins what users see, and a line that must stay an error
  carries `# type: ignore[code]`, which fails once it is unused. Both hooks
  have the id `mypy`, so `pre-commit run mypy` runs both. Both pin
  `numpy==2.2.6`, the last release supporting Python 3.10, which the CI
  style job uses, so local and CI results agree; bump it deliberately. Under
  Python 3.14 the first hook install builds that numpy from source, which
  takes several minutes, once.
- **Conventions.** PEP 604 unions (`X | None`). No
  `from __future__ import annotations`, as upstream evaluates annotations
  eagerly: quote forward references and import them under
  `if TYPE_CHECKING:`, as `Bio/SeqRecord.py` does for `SeqFeature`. A method
  returning its own class uses a bound `TypeVar`, as `Bio/PDB/Entity.py`
  does, not `typing.Self` (Python 3.11+); `typing_extensions` is not a
  runtime dependency. An attribute that is `None` only on a blank object
  still being built is typed `X | Any`, typeshed's trick, so users need not
  narrow it; the constructor parameter stays honestly `X | None`. Overloads
  follow the return types measured at runtime, function by function, not a
  blanket rule.
- **Annotating does not change the public API.** Dump the signatures on
  `main` and on your branch, from each checkout, and diff them. Any
  difference must be intended and named in the PR. The dump covers each
  public name a module defines, and the public and dunder members of its
  public classes, inherited ones included:

  ```python
  import importlib
  import inspect
  import json
  import sys


  def params(obj):
      obj = getattr(obj, "__func__", obj)  # classmethod, staticmethod
      obj = getattr(obj, "fget", obj)  # property
      try:
          parameters = inspect.signature(obj).parameters.values()
      except (TypeError, ValueError):  # not callable, or no signature
          return None
      return [(p.name, p.kind.name, repr(p.default)) for p in parameters]


  dump = {}
  for name in sys.argv[1:]:
      module = importlib.import_module(name)
      for key, value in vars(module).items():
          if key.startswith("_") or getattr(value, "__module__", None) != name:
              continue
          dump[f"{name}.{key}"] = params(value)
          if inspect.isclass(value):
              for attr in dir(value):
                  if not attr.startswith("_") or attr.endswith("__"):
                      member = inspect.getattr_static(value, attr)
                      dump[f"{name}.{key}.{attr}"] = params(member)
  json.dump(dump, sys.stdout, indent=1, sort_keys=True)
  ```

  Save it outside the checkout, then for example
  `PYTHONPATH=. python /tmp/sigdump.py Bio.Seq Bio.SeqRecord > /tmp/branch.json`.
- **Docstrings.** Keep the house list format: an `Arguments:` heading (or
  `Keyword arguments:`), then ` - name - description` items. Types reach the
  docs through the signatures, but do not strip type words from the
  descriptions. Convert `:param:` or numpydoc sections only in a module you
  are annotating anyway.

## Pull requests

1. Branch from `main` with a short descriptive name. Do not work on `main`.
2. Make the change; add tests.
3. `cd Tests && python run_tests.py --offline`
4. `pre-commit run --all-files`
5. Open the PR and describe what you verified and what you did not. There is
   no licence statement to include; see [Licensing](#licensing).

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
