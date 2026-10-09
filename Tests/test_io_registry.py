# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the lazy format tables of Bio.SeqIO, Bio.Align and Bio.Phylo."""

import importlib
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import types
import unittest
from io import StringIO

import support

from Bio import Align
from Bio import MissingPythonDependencyError
from Bio import Phylo
from Bio import SeqIO
from Bio._io_registry import _resolve
from Bio._io_registry import FormatRegistry
from Bio.Align import clustal

CLUSTAL = support.DATA / "Clustalw" / "opuntia.aln"
NEWICK = support.DATA / "Nexus" / "int_node_labels.nwk"

# Modules written to a temporary directory, so that each test knows whether
# they have been imported yet.
HANDLERS = "_io_registry_test_handlers"
KEYERROR = "_io_registry_test_keyerror"
MODULE_SOURCES = {
    HANDLERS: "class Handler:\n    class Inner:\n        pass\n",
    KEYERROR: 'raise KeyError("raised while importing")\n',
}


def setUpModule():
    global module_dir
    module_dir = tempfile.mkdtemp()
    for name, source in MODULE_SOURCES.items():
        with open(os.path.join(module_dir, name + ".py"), "w") as handle:
            handle.write(source)
    sys.path.insert(0, module_dir)
    importlib.invalidate_caches()


def tearDownModule():
    sys.path.remove(module_dir)
    for name in MODULE_SOURCES:
        sys.modules.pop(name, None)
    shutil.rmtree(module_dir)


def stored_values(registry):
    """Return the table's values as stored, without resolving them."""
    return {name: dict.__getitem__(registry, name) for name in registry}


class TestCaseForgettingHandlers(unittest.TestCase):
    """Start each test with the temporary handler module not yet imported."""

    def setUp(self):
        sys.modules.pop(HANDLERS, None)


class Resolve(TestCaseForgettingHandlers):
    """Spec strings resolve as importlib.metadata.EntryPoint.load does."""

    def test_module(self):
        module = _resolve(HANDLERS)
        self.assertIs(module, sys.modules[HANDLERS])

    def test_attribute_path(self):
        handler = _resolve(f"{HANDLERS}:Handler")
        self.assertIs(handler, sys.modules[HANDLERS].Handler)
        self.assertIs(_resolve(f"{HANDLERS}:Handler.Inner"), handler.Inner)

    def test_dotted_module(self):
        self.assertIs(_resolve("Bio.Align.clustal"), clustal)
        self.assertIs(
            _resolve("Bio.Align.clustal:AlignmentWriter"), clustal.AlignmentWriter
        )

    def test_missing(self):
        with self.assertRaises(ModuleNotFoundError):
            _resolve("_io_registry_test_no_such_module")
        with self.assertRaises(AttributeError):
            _resolve(f"{HANDLERS}:NoSuchAttribute")


class Registry(TestCaseForgettingHandlers):
    """FormatRegistry resolves values on first access and keeps the result."""

    def test_lazy_and_cached(self):
        registry = FormatRegistry({"fmt": f"{HANDLERS}:Handler"})
        # Membership, listing and length do not import anything.
        self.assertIn("fmt", registry)
        self.assertEqual(list(registry), ["fmt"])
        self.assertEqual(len(registry), 1)
        self.assertNotIn(HANDLERS, sys.modules)
        handler = registry["fmt"]
        self.assertIs(handler, sys.modules[HANDLERS].Handler)
        # The resolved value replaces the spec in the table.
        self.assertIs(dict.__getitem__(registry, "fmt"), handler)
        self.assertIs(registry["fmt"], handler)
        self.assertIs(registry.get("fmt"), handler)

    def test_package_prefix(self):
        marker = object()
        specs = {
            "short": "FastaIO.FastaIterator",
            "full": "Bio.Align.clustal:AlignmentWriter",
            "delegated": None,
            "object": marker,
        }
        registry = FormatRegistry(specs, package="Bio.SeqIO")
        self.assertEqual(
            stored_values(registry),
            {
                "short": "Bio.SeqIO.FastaIO:FastaIterator",
                "full": "Bio.Align.clustal:AlignmentWriter",
                "delegated": None,
                "object": marker,
            },
        )
        self.assertIs(registry["short"], SeqIO.FastaIO.FastaIterator)
        self.assertIs(registry["full"], clustal.AlignmentWriter)
        # Without a package, the short form is just a module name.
        registry = FormatRegistry(specs)
        self.assertEqual(stored_values(registry)["short"], "FastaIO.FastaIterator")

    def test_factory(self):
        calls = []

        def factory(name):
            calls.append(name)
            return types.SimpleNamespace(name=name)

        registry = FormatRegistry({"one": None, "two": None}, factory)
        first = registry["one"]
        self.assertEqual(first.name, "one")
        self.assertIs(registry["one"], first)
        self.assertEqual(calls, ["one"])
        self.assertEqual(registry["two"].name, "two")
        self.assertEqual(calls, ["one", "two"])

    def test_other_values_returned_as_is(self):
        marker = object()
        registry = FormatRegistry({"fmt": marker})
        self.assertIs(registry["fmt"], marker)
        registry["fmt"] = len
        self.assertIs(registry["fmt"], len)

    def test_values_and_items_resolve(self):
        registry = FormatRegistry(
            {"fmt": f"{HANDLERS}:Handler", "inner": f"{HANDLERS}:Handler.Inner"}
        )
        self.assertNotIn(HANDLERS, sys.modules)
        values = registry.values()
        handler = sys.modules[HANDLERS].Handler
        self.assertEqual(values, [handler, handler.Inner])
        self.assertEqual(registry.items(), [("fmt", handler), ("inner", handler.Inner)])

    def test_get_absent_name_returns_default(self):
        registry = FormatRegistry({"fmt": f"{HANDLERS}:Handler"})
        self.assertIsNone(registry.get("absent"))
        marker = object()
        self.assertIs(registry.get("absent", marker), marker)
        with self.assertRaises(KeyError):
            registry["absent"]
        self.assertNotIn(HANDLERS, sys.modules)

    def test_get_failing_name_propagates(self):
        """An error while resolving a present name is not 'absent'."""

        def factory(name):
            raise KeyError("raised by the factory")

        registry = FormatRegistry({"spec": KEYERROR, "factory": None}, factory)
        marker = object()
        for name, message in [
            ("spec", "raised while importing"),
            ("factory", "raised by the factory"),
        ]:
            with self.subTest(name=name):
                with self.assertRaises(KeyError) as cm:
                    registry.get(name, marker)
                self.assertEqual(cm.exception.args, (message,))
                with self.assertRaises(KeyError) as cm:
                    registry[name]
                self.assertEqual(cm.exception.args, (message,))

    def test_get_unhashable_name(self):
        registry = FormatRegistry({"fmt": None})
        with self.assertRaises(TypeError):
            registry.get(["fmt"])

    def test_concurrent_first_access_agrees_on_one_value(self):
        """Racing first accesses all get the single value that is stored."""
        nthreads = 8
        barrier = threading.Barrier(nthreads)

        def factory(name):
            time.sleep(0.01)  # widen the race window
            return type("Wrapper", (), {})

        registry = FormatRegistry({"fmt": None}, factory)
        results = []

        def worker():
            barrier.wait()
            results.append(registry["fmt"])

        threads = [threading.Thread(target=worker) for _ in range(nthreads)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(len(results), nthreads)
        stored = dict.__getitem__(registry, "fmt")
        for result in results:
            self.assertIs(result, stored)

    def test_assignment_during_resolution_wins(self):
        """A value assigned while a spec resolves is kept, not overwritten."""
        replacement = types.SimpleNamespace(name="replacement")

        def factory(name):
            registry[name] = replacement
            return types.SimpleNamespace(name="resolved")

        registry = FormatRegistry({"fmt": None}, factory)
        self.assertIs(registry["fmt"], replacement)
        self.assertIs(dict.__getitem__(registry, "fmt"), replacement)

    def test_spec_assigned_during_resolution_is_resolved(self):
        """A spec assigned while another resolves is itself resolved."""

        def factory(name):
            registry[name] = f"{HANDLERS}:Handler"
            return types.SimpleNamespace(name="resolved")

        registry = FormatRegistry({"fmt": None}, factory)
        handler = registry["fmt"]
        self.assertIs(handler, sys.modules[HANDLERS].Handler)
        self.assertIs(dict.__getitem__(registry, "fmt"), handler)


class BuiltinTables(unittest.TestCase):
    """Every built-in entry of Bio.SeqIO and Bio.Align resolves."""

    def test_seqio(self):
        for table in [SeqIO._FormatToIterator, SeqIO._FormatToWriter, SeqIO._converter]:
            self.assertIsInstance(table, FormatRegistry)
            for name in table:
                with self.subTest(name=name):
                    self.assertTrue(callable(table[name]))

    def test_align(self):
        self.assertIsInstance(Align._registry, FormatRegistry)
        self.assertEqual(list(Align._registry), list(Align.formats))
        for fmt in Align.formats:
            with self.subTest(fmt=fmt):
                module = Align._registry[fmt]
                if fmt == "phylip-relaxed":
                    # A format name that is not a module name.
                    expected = importlib.import_module("Bio.Align.phylip")._relaxed
                else:
                    expected = importlib.import_module(f"Bio.Align.{fmt}")
                self.assertIs(module, expected)
                self.assertTrue(hasattr(module, "AlignmentIterator"))

    def test_phylo(self):
        table = Phylo._io.supported_formats
        self.assertIsInstance(table, FormatRegistry)
        expected = {
            "newick": "NewickIO",
            "nexus": "NexusIO",
            "phyloxml": "PhyloXMLIO",
            "nexml": "NeXMLIO",
        }
        # cdao is offered wherever rdflib is installed, even an rdflib that
        # CDAOIO cannot use; using it then raises MissingPythonDependencyError.
        if importlib.util.find_spec("rdflib") is not None:
            expected["cdao"] = "CDAOIO"
        self.assertEqual(list(table), list(expected))
        for fmt, name in expected.items():
            with self.subTest(fmt=fmt):
                try:
                    module = table[fmt]
                except MissingPythonDependencyError as err:
                    self.skipTest(str(err))
                self.assertIs(module, importlib.import_module(f"Bio.Phylo.{name}"))
                self.assertTrue(callable(module.parse))
                self.assertTrue(callable(module.write))

    def test_import_seqio_stays_light(self):
        """Importing Bio.SeqIO needs neither Bio.Align nor importlib.metadata."""
        code = (
            "import sys, Bio.SeqIO\n"
            "heavy = ['numpy', 'Bio.Align', 'Bio.AlignIO', 'importlib.metadata']\n"
            "loaded = [name for name in heavy if name in sys.modules]\n"
            "assert not loaded, 'import Bio.SeqIO pulled in %s' % loaded\n"
        )
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code],
            capture_output=True,
            text=True,
            cwd=support.DATA,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


class PhyloImportsFormatsOnDemand(unittest.TestCase):
    """Bio.Phylo imports a tree format module only when it is used.

    Each check runs in a fresh interpreter, since this one has already
    imported whatever the rest of the test suite needed.
    """

    def run_python(self, code):
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code],
            capture_output=True,
            text=True,
            cwd=support.DATA,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_reading_newick_stays_light(self):
        """Reading a Newick tree needs none of NumPy, Bio.Align, Bio.Nexus or rdflib."""
        self.run_python(
            "import sys\n"
            "from Bio import Phylo\n"
            f"Phylo.read({str(NEWICK)!r}, 'newick')\n"
            "heavy = ['numpy', 'Bio.Align', 'Bio.Nexus', 'rdflib']\n"
            "loaded = [name for name in heavy if name in sys.modules]\n"
            "assert not loaded, 'reading Newick pulled in %s' % loaded\n"
        )

    def test_submodules_are_still_attributes(self):
        """Each submodule that importing Bio.Phylo used to bind is still there."""
        self.run_python(
            "import Bio.Phylo\n"
            "names = ['BaseTree', 'CDAO', 'NeXML', 'NeXMLIO', 'Newick', 'NewickIO',"
            " 'NexusIO', 'PhyloXML', 'PhyloXMLIO', '_cdao_owl', '_io', '_utils']\n"
            "for name in names:\n"
            "    module = getattr(Bio.Phylo, name)\n"
            "    assert module.__name__ == 'Bio.Phylo.' + name, module\n"
        )

    def test_without_rdflib(self):
        """Without rdflib, cdao is not offered and CDAOIO is not an attribute."""
        self.run_python(
            "import sys\n"
            "sys.modules['rdflib'] = None  # as if rdflib were not installed\n"
            "import Bio.Phylo\n"
            "from Bio import MissingPythonDependencyError\n"
            "assert 'cdao' not in Bio.Phylo._io.supported_formats\n"
            "assert not hasattr(Bio.Phylo, 'CDAOIO')\n"
            "try:\n"
            "    Bio.Phylo.CDAOIO\n"
            "except AttributeError as err:\n"
            "    assert isinstance(err.__cause__, MissingPythonDependencyError)\n"
            "try:\n"
            "    from Bio.Phylo import CDAOIO\n"
            "except MissingPythonDependencyError:\n"
            "    pass\n"
            "else:\n"
            "    raise AssertionError('imported CDAOIO without rdflib')\n"
        )

    def test_rdflib_imported_without_spec(self):
        """An rdflib already in sys.modules without a __spec__ counts as present."""
        self.run_python(
            "import sys, types\n"
            "sys.modules['rdflib'] = types.ModuleType('rdflib')  # no __spec__\n"
            "import Bio.Phylo\n"
            "assert 'cdao' in Bio.Phylo._io.supported_formats\n"
        )


class SeqIOErrorWhileResolving(unittest.TestCase):
    """Bio.SeqIO does not report a failing format as an unknown one."""

    def setUp(self):
        SeqIO._FormatToIterator["test-keyerror"] = f"{KEYERROR}:Iterator"
        SeqIO._FormatToWriter["test-keyerror"] = f"{KEYERROR}:Writer"

    def tearDown(self):
        del SeqIO._FormatToIterator["test-keyerror"]
        del SeqIO._FormatToWriter["test-keyerror"]

    def test_parse(self):
        with self.assertRaises(KeyError) as cm:
            SeqIO.parse(StringIO(""), "test-keyerror")
        self.assertEqual(cm.exception.args, ("raised while importing",))

    def test_write(self):
        with self.assertRaises(KeyError) as cm:
            SeqIO.write([], StringIO(), "test-keyerror")
        self.assertEqual(cm.exception.args, ("raised while importing",))

    def test_absent_format_is_still_unknown(self):
        with self.assertRaises(ValueError) as cm:
            SeqIO.parse(StringIO(""), "test-absent")
        self.assertEqual(cm.exception.args, ("Unknown format 'test-absent'",))


class AlignNameThatIsNotAModule(unittest.TestCase):
    """Bio.Align looks names up in its registry, not as module names."""

    def setUp(self):
        Align._registry["test-variant"] = types.SimpleNamespace(
            AlignmentIterator=clustal.AlignmentIterator,
            AlignmentWriter=clustal.AlignmentWriter,
        )
        self.alignment = Align.read(CLUSTAL, "clustal")

    def tearDown(self):
        del Align._registry["test-variant"]

    def test_not_a_builtin_format(self):
        self.assertNotIn("test-variant", Align.formats)

    def test_parse(self):
        alignments = list(Align.parse(CLUSTAL, "Test-Variant"))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0].shape, self.alignment.shape)

    def test_read(self):
        alignment = Align.read(CLUSTAL, "TEST-variant")
        self.assertEqual(alignment.shape, self.alignment.shape)

    def test_write(self):
        stream = StringIO()
        self.assertEqual(Align.write(self.alignment, stream, "test-VARIANT"), 1)
        expected = StringIO()
        Align.write(self.alignment, expected, "clustal")
        self.assertEqual(stream.getvalue(), expected.getvalue())

    def test_format(self):
        expected = self.alignment.format("clustal")
        self.assertEqual(self.alignment.format("Test-Variant"), expected)
        self.assertEqual(format(self.alignment, "test-variant"), expected)

    def test_import_error_is_not_an_unknown_format(self):
        Align._registry["test-broken"] = KEYERROR
        try:
            with self.assertRaises(KeyError) as cm:
                Align.parse(CLUSTAL, "Test-Broken")
            self.assertEqual(cm.exception.args, ("raised while importing",))
        finally:
            del Align._registry["test-broken"]


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
