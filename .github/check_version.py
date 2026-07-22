#!/usr/bin/env python3
"""Check that the release tag matches ``Bio.__version__``.

PyPI versions are immutable: once a version is uploaded, that exact string can
never be reused, even if the release is deleted. The git tag does not set the
version - ``pyproject.toml`` reads it from ``Bio.__version__`` - so a mismatch
silently publishes something other than what the tag says. This script turns
that into a loud failure before anything is built.

Usage::

    python .github/check_version.py            # reads GITHUB_REF_NAME
    python .github/check_version.py v1.88.dev0 # or an explicit tag

Exits 0 when the ref is not a version tag, so the same job can run on
``workflow_dispatch`` without special-casing in the workflow.
"""

from __future__ import annotations

import os
import pathlib
import re
import sys

from packaging.version import InvalidVersion
from packaging.version import Version

INIT = pathlib.Path(__file__).resolve().parent.parent / "Bio" / "__init__.py"


def package_version() -> str:
    """Return __version__ from Bio/__init__.py without importing Bio.

    Importing would drag in NumPy and the compiled extensions, which are not
    necessarily built in the job that runs this check.
    """
    match = re.search(
        r'^__version__\s*=\s*["\']([^"\']+)["\']', INIT.read_text(), re.MULTILINE
    )
    if match is None:
        sys.exit(f"::error::Could not find __version__ in {INIT}")
    return match.group(1)


def main(argv: list[str]) -> int:
    """Compare the tag against the package version and report."""
    tag = argv[1] if len(argv) > 1 else os.environ.get("GITHUB_REF_NAME", "")

    if not tag.startswith("v"):
        print(f"Ref {tag!r} is not a version tag; skipping the version check.")
        return 0

    declared = package_version()
    try:
        tag_version = Version(tag[1:])
        pkg_version = Version(declared)
    except InvalidVersion as err:
        sys.exit(f"::error::Not a valid PEP 440 version: {err}")

    if tag_version != pkg_version:
        sys.exit(
            f"::error::Tag {tag} does not match Bio.__version__ = {declared}.\n"
            f"The tag does not set the version - pyproject.toml reads it from "
            f"Bio/__init__.py - so this would publish {declared}, not {tag[1:]}. "
            f"PyPI versions are immutable, so fix one of them before releasing."
        )

    print(f"OK: tag {tag} matches Bio.__version__ = {declared} (as {pkg_version}).")
    if pkg_version.is_prerelease:
        print(
            f"::notice::{pkg_version} is a pre-release. `pip install biopaithon` "
            f"will not select it; users need `pip install --pre biopaithon`."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
