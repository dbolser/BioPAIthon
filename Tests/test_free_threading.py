# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests that the C extensions declare that they can run without the GIL.

On a free-threaded build of CPython, importing an extension that has not
declared itself safe without the GIL turns the GIL back on for the whole
process. Every extension that pyproject.toml builds must make that
declaration in its PyInit function, through Bio/_freethreading.h, and be
listed below in DECLARED. So a new extension must be safe without the GIL.
"""

import importlib.machinery
import os
import subprocess
import sys
import sysconfig
import unittest

try:
    import tomllib
except ImportError:  # Python 3.10
    tomllib = None

import support

import Bio

DECLARED = {
    "Bio.Align._aligncore",
    "Bio.Align._alignmentcounts",
    "Bio.Align._codonaligner",
    "Bio.Align._pairwisealigner",
    "Bio.Align.substitution_matrices._arraycore",
    "Bio.Cluster._cluster",
    "Bio.Nexus.cnexus",
    "Bio.PDB._bcif_helper",
    "Bio.PDB.ccealign",
    "Bio.PDB.kdtrees",
    "Bio.SeqIO._twoBitIO",
    "Bio.motifs._pwm",
}

# Loads one extension module from its file, bypassing the package __init__
# (which could import other extensions), and reports the GIL state after.
LOAD_ALONE = """\
import importlib.machinery
import sys

name, path = sys.argv[1:]
loader = importlib.machinery.ExtensionFileLoader(name, path)
spec = importlib.machinery.ModuleSpec(name, loader, origin=path)
loader.exec_module(loader.create_module(spec))
print(sys._is_gil_enabled())
"""

# Imports the subpackages that load the C extensions, as a user's program
# would, and reports the GIL state after.
IMPORT_ALL = """\
import sys

import Bio.Align
import Bio.AlignIO
import Bio.Cluster
import Bio.motifs
import Bio.Nexus.Nexus
import Bio.PDB
import Bio.SeqIO

print(sys._is_gil_enabled())
"""


def run_without_gil_override(*args):
    """Run Python with the given arguments, with no PYTHON_GIL set.

    PYTHON_GIL=0 would keep the GIL off whatever the modules imported
    declare.
    """
    env = dict(os.environ)
    env.pop("PYTHON_GIL", None)
    return subprocess.run(
        [sys.executable, *args], env=env, capture_output=True, text=True
    )


class ClassificationTests(unittest.TestCase):
    """Every extension is DECLARED."""

    @unittest.skipIf(tomllib is None, "tomllib is new in Python 3.11")
    def test_every_extension_is_declared(self):
        """Check DECLARED lists exactly the extensions built."""
        with open(support.DATA.parent / "pyproject.toml", "rb") as handle:
            config = tomllib.load(handle)
        built = {ext["name"] for ext in config["tool"]["setuptools"]["ext-modules"]}
        self.assertEqual(DECLARED, built)


@unittest.skipUnless(
    sysconfig.get_config_var("Py_GIL_DISABLED"), "requires a free-threaded build"
)
class LoadAloneTests(unittest.TestCase):
    """Each extension, loaded alone, leaves the GIL disabled.

    A test method is added below for each module in DECLARED, so that a
    verbose run names every module checked.
    """

    def gil_enabled_after_loading(self, name):
        package, _, module = name.rpartition(".")
        directory = os.path.join(os.path.dirname(Bio.__file__), *package.split(".")[1:])
        spec = importlib.machinery.PathFinder.find_spec(module, [directory])
        self.assertIsNotNone(spec, f"{name} is not built in {directory}")
        result = run_without_gil_override("-c", LOAD_ALONE, name, spec.origin)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(result.stdout.strip(), ("True", "False"), result.stderr)
        return result.stdout.strip() == "True"


def _add_load_alone_test(name):
    def test(self):
        self.assertIs(self.gil_enabled_after_loading(name), False)

    test.__doc__ = f"Check loading {name} alone leaves the GIL disabled."
    method = "test_" + name.replace(".", "_")
    setattr(LoadAloneTests, method, test)


for _name in sorted(DECLARED):
    _add_load_alone_test(_name)


@unittest.skipUnless(
    sysconfig.get_config_var("Py_GIL_DISABLED"), "requires a free-threaded build"
)
class ImportTests(unittest.TestCase):
    """Importing the subpackages, as a program would, leaves the GIL disabled."""

    def test_import_subpackages(self):
        """Check importing the subpackages that use C extensions keeps the GIL off.

        On failure, CPython's RuntimeWarning in stderr names the module that
        turned the GIL on.
        """
        result = run_without_gil_override("-c", IMPORT_ALL)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "False", result.stderr)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
