# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Check that each C extension's .pyi stub declares everything it exports.

mypy's stubtest compares the members of each class a stub declares with the
runtime class, but it never notices a class missing from a stub. It skips a
runtime name whose __module__ is not the module's own name, taking it for an
import, and a C type's __module__ comes from its tp_name. No C type in Bio
has a fully dotted tp_name ("_pairwisealigner.PairwiseAligner" gives
"_pairwisealigner", "AlignmentCounts" gives "builtins"), so stubtest skips
them all. This test closes that gap: every extension module in
pyproject.toml must have a .pyi beside its C source, and every public name in
dir(module) must be declared at the top level of that stub.
"""

import ast
import importlib
import os
import unittest

try:
    import tomllib
except ImportError:  # Python 3.10
    tomllib = None

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PYPROJECT = os.path.join(ROOT, "pyproject.toml")

# Extensions not stubbed yet. Bio.Cluster._cluster gets its stub together with
# the annotations of Bio/Cluster/__init__.py. An entry that has gained a stub
# fails the test below, so delete it in the pull request that adds the stub.
UNSTUBBED = {"Bio.Cluster._cluster"}


def _extension_modules():
    """Return the dotted names of the C extensions in pyproject.toml."""
    with open(PYPROJECT, "rb") as handle:
        config = tomllib.load(handle)
    return [ext["name"] for ext in config["tool"]["setuptools"]["ext-modules"]]


def _stub_path(module_name):
    """Return the path of the .pyi stub for the dotted module name."""
    return os.path.join(ROOT, *module_name.split(".")) + ".pyi"


def _declared_names(path):
    """Return the names a stub file declares at its top level."""
    with open(path, encoding="utf-8") as handle:
        tree = ast.parse(handle.read(), filename=path)
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


@unittest.skipIf(tomllib is None, "needs tomllib, new in Python 3.11")
class StubCompletenessTests(unittest.TestCase):
    """Check each C extension has a stub declaring every public name it has."""

    def test_every_extension_has_a_stub(self):
        """Every extension in pyproject.toml, bar UNSTUBBED, has a .pyi stub."""
        extensions = _extension_modules()
        # Guard against passing vacuously, if pyproject.toml were restructured.
        self.assertTrue(extensions, "found no extension module in pyproject.toml")
        unstubbed = sorted(
            name for name in extensions if not os.path.isfile(_stub_path(name))
        )
        self.assertEqual(
            unstubbed,
            sorted(UNSTUBBED),
            "each C extension needs a .pyi stub named after it, beside its "
            "C source; UNSTUBBED must list exactly those still without one",
        )

    def test_stubbed_extensions_declare_public_names(self):
        """Every public name of a stubbed extension is in its stub."""
        stubbed = [
            name for name in _extension_modules() if os.path.isfile(_stub_path(name))
        ]
        # Guard against passing vacuously, if pyproject.toml were restructured.
        self.assertTrue(stubbed, "found no extension module with a .pyi stub")
        for name in stubbed:
            with self.subTest(module=name):
                module = importlib.import_module(name)
                public = {n for n in dir(module) if not n.startswith("_")}
                missing = sorted(public - _declared_names(_stub_path(name)))
                self.assertEqual(
                    missing, [], f"{_stub_path(name)} does not declare these"
                )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
