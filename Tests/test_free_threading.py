# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests of which C extensions declare that they can run without the GIL.

On a free-threaded build of CPython, importing an extension that has not
declared itself safe without the GIL turns the GIL back on for the whole
process. Every extension that pyproject.toml builds is listed below as
either DECLARED (its PyInit function makes that declaration, through
Bio/_freethreading.h) or PENDING (not yet). When an extension is made safe
and declared, move it from PENDING to DECLARED.
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
    "Bio.Nexus.cnexus",
    "Bio.PDB._bcif_helper",
    "Bio.PDB.ccealign",
    "Bio.PDB.kdtrees",
    "Bio.SeqIO._twoBitIO",
    "Bio.motifs._pwm",
}

PENDING = {
    "Bio.Align._aligncore",
    "Bio.Align._alignmentcounts",
    "Bio.Align._codonaligner",
    "Bio.Align._pairwisealigner",
    "Bio.Align.substitution_matrices._arraycore",
    "Bio.Cluster._cluster",
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


class ClassificationTests(unittest.TestCase):
    """Every extension is either DECLARED or PENDING, and not both."""

    @unittest.skipIf(tomllib is None, "tomllib is new in Python 3.11")
    def test_every_extension_is_classified(self):
        """Check DECLARED and PENDING together are the extensions built."""
        with open(support.DATA.parent / "pyproject.toml", "rb") as handle:
            config = tomllib.load(handle)
        built = {ext["name"] for ext in config["tool"]["setuptools"]["ext-modules"]}
        self.assertEqual(DECLARED | PENDING, built)

    def test_no_extension_is_both(self):
        """Check no extension is listed as both DECLARED and PENDING."""
        self.assertEqual(DECLARED & PENDING, set())


@unittest.skipUnless(
    sysconfig.get_config_var("Py_GIL_DISABLED"), "requires a free-threaded build"
)
class LoadAloneTests(unittest.TestCase):
    """Each extension, loaded alone, leaves the GIL as its list says.

    A test method is added below for each module in DECLARED and PENDING,
    so that a verbose run names every module checked. The PENDING ones
    show that loading an undeclared module does turn the GIL on here, so
    a pass for a DECLARED module is not an artefact of how it was loaded.
    """

    def gil_enabled_after_loading(self, name):
        package, _, module = name.rpartition(".")
        directory = os.path.join(os.path.dirname(Bio.__file__), *package.split(".")[1:])
        spec = importlib.machinery.PathFinder.find_spec(module, [directory])
        self.assertIsNotNone(spec, f"{name} is not built in {directory}")
        env = dict(os.environ)
        # PYTHON_GIL=0 would keep the GIL off whatever the module declares.
        env.pop("PYTHON_GIL", None)
        result = subprocess.run(
            [sys.executable, "-c", LOAD_ALONE, name, spec.origin],
            env=env,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(result.stdout.strip(), ("True", "False"), result.stderr)
        return result.stdout.strip() == "True"


def _add_load_alone_test(name, declared):
    def test(self):
        self.assertIs(self.gil_enabled_after_loading(name), not declared)

    if declared:
        test.__doc__ = f"Check loading {name} alone leaves the GIL disabled."
    else:
        test.__doc__ = f"Check loading {name} alone enables the GIL."
    method = "test_" + name.replace(".", "_")
    setattr(LoadAloneTests, method, test)


for _name in sorted(DECLARED):
    _add_load_alone_test(_name, declared=True)
for _name in sorted(PENDING):
    _add_load_alone_test(_name, declared=False)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
