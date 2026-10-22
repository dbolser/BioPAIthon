# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests that removed modules explain themselves instead of vanishing.

``Bio.Application`` and the command line wrappers built on it were removed in
release 1.86.  Code that imports them is still common -- a GitHub code search
finds thousands of public files -- and without a stub each one gets a bare
``ModuleNotFoundError`` naming a module and nothing else.  ``Bio.Alphabet``,
removed in 1.78, has always been kept as a stub for this reason; these follow
it.

This is a deliberate difference from upstream Biopython, which ships no stub
for any of them.
"""

import ast
import os
import unittest

import Bio
from conftest import REMOVED_MODULES

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
DOC_CONF = os.path.join(TESTS_DIR, os.pardir, "Doc", "conf.py")

# A phrase each stub's error message must contain, so that a message cannot
# be reworded into uselessness without failing here.  The modules themselves
# are listed once, in Tests/conftest.py, which the doctest collector also
# reads.
REPLACEMENTS = {
    "Bio.Align.AlignInfo": "SummaryInfo",
    "Bio.Alphabet": "molecule_type",
    "Bio.Application": "subprocess",
    "Bio.Align.Applications": "MuscleCommandline",
    "Bio.Blast.Applications": "NcbiblastnCommandline",
    "Bio.Emboss.Applications": "WaterCommandline",
    "Bio.HMM": "hmmlearn",
    "Bio.Phylo.Applications": "PhymlCommandline",
    "Bio.Sequencing.Applications": "Novoalign",
    "Bio.pairwise2": "PairwiseAligner",
}


def is_stub(path):
    """Return True if the module at path does nothing but raise ImportError."""
    with open(path, "rb") as handle:
        tree = ast.parse(handle.read(), path)
    body = tree.body
    if ast.get_docstring(tree, clean=False) is not None:
        body = body[1:]
    return (
        len(body) == 1
        and isinstance(body[0], ast.Raise)
        and isinstance(body[0].exc, ast.Call)
        and isinstance(body[0].exc.func, ast.Name)
        and body[0].exc.func.id == "ImportError"
    )


def find_stubs():
    """Return the names of the stub modules in the Bio package being tested."""
    package = os.path.dirname(os.path.abspath(Bio.__file__))
    stubs = set()
    for root, dirs, files in os.walk(package):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for filename in files:
            if not filename.endswith(".py"):
                continue
            path = os.path.join(root, filename)
            if not is_stub(path):
                continue
            relative = os.path.relpath(path, os.path.dirname(package))
            parts = relative[: -len(".py")].split(os.sep)
            if parts[-1] == "__init__":
                parts.pop()
            stubs.add(".".join(parts))
    return stubs


def doc_conf_excludes():
    """Return the Bio modules Doc/conf.py keeps out of the API documentation."""
    with open(DOC_CONF, "rb") as handle:
        tree = ast.parse(handle.read(), DOC_CONF)
    excluded = set()
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
            continue
        value = node.value
        if not value.startswith("../Bio/"):
            continue
        value = value[len("../") :]
        if value.endswith("/"):
            value = value[:-1]
        elif value.endswith(".py"):
            value = value[: -len(".py")]
        excluded.add(value.replace("/", "."))
    return excluded


class TestRemovedModulesList(unittest.TestCase):
    """The one list in conftest.py agrees with the tree and with Doc/conf.py."""

    def test_list_matches_the_stubs_in_the_tree(self):
        self.assertEqual(sorted(REMOVED_MODULES), sorted(find_stubs()))

    def test_each_stub_has_a_replacement_phrase(self):
        self.assertEqual(sorted(REMOVED_MODULES), sorted(REPLACEMENTS))

    def test_doc_conf_excludes_each_stub(self):
        # The docs build cannot import from Tests/, so Doc/conf.py keeps its
        # own list of what autodoc must not import.  Check it is not short.
        if not os.path.exists(DOC_CONF):
            self.skipTest("Doc/conf.py not present")
        missing = set(REMOVED_MODULES) - doc_conf_excludes()
        self.assertEqual(set(), missing)


class TestRemovedModules(unittest.TestCase):
    """Importing a removed module raises ImportError, not ModuleNotFoundError."""

    def test_raises_import_error(self):
        for name in REMOVED_MODULES:
            with self.subTest(module=name):
                with self.assertRaises(ImportError) as cm:
                    __import__(name)
                # A missing module raises ModuleNotFoundError, a subclass of
                # ImportError, so assertRaises alone would pass on no stub
                # at all.  The point of the stub is that it is *not* that.
                self.assertNotIsInstance(cm.exception, ModuleNotFoundError)

    def test_message_names_the_replacement(self):
        for name in REMOVED_MODULES:
            with self.subTest(module=name):
                with self.assertRaises(ImportError) as cm:
                    __import__(name)
                message = str(cm.exception)
                self.assertIn(name, message)
                self.assertIn(REPLACEMENTS[name], message)

    def test_still_importable_by_a_guarded_caller(self):
        """`except ImportError` around the old import must keep working."""
        for name in REMOVED_MODULES:
            with self.subTest(module=name):
                try:
                    __import__(name)
                except ImportError:
                    continue
                self.fail(f"{name} imported without raising")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
