# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for file formats which installed distributions declare as entry points.

Bio.SeqIO and Bio.Align look for entry points once per process, so each test
runs its checks in a fresh interpreter.  That interpreter finds three fake
distributions, fakeplugin, otherplugin and nameless, in a temporary directory
on sys.path, along with the modules their entry points name.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import support

FASTA = str(support.DATA / "Fasta" / "f002")
FASTA_ONE = str(support.DATA / "Fasta" / "f001")
CLUSTAL = str(support.DATA / "Clustalw" / "opuntia.aln")

# Distribution name to the text of its entry_points.txt file.  The metadata
# of the distribution "nameless" gives no Name.
DISTRIBUTIONS = {
    "fakeplugin": """\
[biopaithon.seqio.iterators]
starfasta = fakeplugin_formats:StarFastaIterator
fasta = fakeplugin_formats:StarFastaIterator
BadName = fakeplugin_formats:StarFastaIterator
dup = fakeplugin_formats:StarFastaIterator
same = fakeplugin_formats:StarFastaIterator
broken = fakeplugin_broken:BrokenIterator

[biopaithon.seqio.writers]
starfasta = fakeplugin_formats:StarFastaWriter [extra]

[biopaithon.align]
STAR-ALN = fakeplugin_aln
GoodAln = fakeplugin_aln
GOODALN = otherplugin_formats
""",
    "otherplugin": """\
[biopaithon.seqio.iterators]
dup = otherplugin_formats:OtherIterator
same = fakeplugin_formats:StarFastaIterator
weird = not a valid!!value
""",
    "nameless": """\
[biopaithon.seqio.iterators]
noname = fakeplugin_formats:StarFastaIterator
NoName = fakeplugin_formats:StarFastaIterator
""",
}

# The modules the entry points name.
MODULES = {
    "fakeplugin_formats": """\
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.FastaIO import FastaWriter


class StarFastaIterator(FastaIterator):
    def __next__(self):
        record = super().__next__()
        record.annotations["parsed by"] = "fakeplugin"
        return record


class StarFastaWriter(FastaWriter):
    pass
""",
    "fakeplugin_aln": (
        "from Bio.Align.clustal import AlignmentIterator\n"
        "from Bio.Align.clustal import AlignmentWriter\n"
    ),
    "fakeplugin_broken": 'raise ImportError("fakeplugin_broken needs a library")\n',
    "otherplugin_formats": (
        "from Bio.SeqIO.FastaIO import FastaIterator as OtherIterator\n"
    ),
}

# Start of the code run in each fresh interpreter.  Warnings given while Bio
# is imported are hidden.  After that, any warning other than a
# BiopythonWarning is an error, such as a DeprecationWarning from
# importlib.metadata, and every BiopythonWarning is printed, however often
# it is given.
PRELUDE = """\
import sys
import warnings
from io import StringIO

sys.path.insert(0, {plugins!r})

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import Align
    from Bio import BiopythonWarning
    from Bio import SeqIO

warnings.simplefilter("error")
warnings.simplefilter("always", BiopythonWarning)
"""

GROUP = "'biopaithon.seqio.iterators'"

# The warnings given when Bio.SeqIO looks for its iterator plugins, sorted as
# run_python returns them.  None is given for "same", which both distributions
# name the same, nor for "noname".
ITERATOR_WARNINGS = [
    f"Ignoring entry point 'BadName' of fakeplugin in group {GROUP}:"
    " Format string 'BadName' should be lower case",
    f"Ignoring entry point 'NoName' of a distribution with no name in group"
    f" {GROUP}: Format string 'NoName' should be lower case",
    f"Ignoring entry point 'fasta' of fakeplugin in group {GROUP}:"
    " 'fasta' is a built-in format",
    f"Ignoring entry point 'weird' of otherplugin in group {GROUP}:"
    " 'not a valid!!value' is not an object reference",
    f"Ignoring entry points 'dup' of fakeplugin and 'dup' of otherplugin in group"
    f" {GROUP}: they give format 'dup' different objects",
]

# The warning given when Bio.Align looks for its plugins.  Its names are
# case-insensitive, so the two entry points name the same format.
ALIGN_WARNINGS = [
    "Ignoring entry points 'GOODALN' of fakeplugin and 'GoodAln' of fakeplugin"
    " in group 'biopaithon.align': they give format 'goodaln' different objects"
]

ALL_WARNINGS = sorted(ITERATOR_WARNINGS + ALIGN_WARNINGS)


def setUpModule():
    global plugins
    plugins = tempfile.mkdtemp()
    for distribution, entry_points in DISTRIBUTIONS.items():
        directory = os.path.join(plugins, f"{distribution}-0.1.dist-info")
        os.mkdir(directory)
        name = "" if distribution == "nameless" else f"Name: {distribution}\n"
        with open(os.path.join(directory, "METADATA"), "w") as handle:
            handle.write(f"Metadata-Version: 2.1\n{name}Version: 0.1\n")
        with open(os.path.join(directory, "entry_points.txt"), "w") as handle:
            handle.write(entry_points)
    for name, source in MODULES.items():
        with open(os.path.join(plugins, name + ".py"), "w") as handle:
            handle.write(source)


def tearDownModule():
    shutil.rmtree(plugins)


class PluginTestCase(unittest.TestCase):
    """Run code in a fresh interpreter which can find the fake plugins."""

    def run_python(self, code):
        """Run the code, and return the BiopythonWarning messages it printed."""
        env = dict(os.environ)
        env.pop("PYTHONWARNINGS", None)
        result = subprocess.run(
            [sys.executable, "-c", PRELUDE.format(plugins=plugins) + code],
            capture_output=True,
            text=True,
            cwd=support.DATA,
            env=env,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return sorted(re.findall(r"BiopythonWarning: (.*)", result.stderr))


class SeqIOPlugin(PluginTestCase):
    """Every SeqIO function works with a format from an entry point."""

    def test_round_trip(self):
        warnings = self.run_python(
            f"""\
fasta = list(SeqIO.parse({FASTA!r}, "fasta"))
ids = [record.id for record in fasta]
records = list(SeqIO.parse({FASTA!r}, "starfasta"))
assert [record.id for record in records] == ids
for record in records:
    assert record.annotations["parsed by"] == "fakeplugin"
record = SeqIO.read({FASTA_ONE!r}, "starfasta")
assert record.annotations["parsed by"] == "fakeplugin"

assert SeqIO._FormatToWriter["starfasta"].__name__ == "StarFastaWriter"
expected = format(fasta[0], "fasta")
assert format(fasta[0], "starfasta") == expected
assert fasta[0].format("starfasta") == expected
handle = StringIO()
SeqIO.write(fasta, handle, "fasta")
expected = handle.getvalue()
handle = StringIO()
assert SeqIO.write(fasta, handle, "starfasta") == 3
assert handle.getvalue() == expected
handle = StringIO()
assert SeqIO.convert({FASTA!r}, "starfasta", handle, "starfasta") == 3
assert handle.getvalue() == expected

index = SeqIO.index({FASTA!r}, "starfasta")
assert list(index) == ids
assert index[ids[1]].annotations["parsed by"] == "fakeplugin"
index.close()
index = SeqIO.index_db(":memory:", {FASTA!r}, "starfasta")
assert sorted(index) == sorted(ids)
assert index[ids[2]].annotations["parsed by"] == "fakeplugin"
index.close()
"""
        )
        self.assertEqual(warnings, ITERATOR_WARNINGS)

    def test_index_db_reload(self):
        """A fresh process reopens an index_db database of a plugin format."""
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "starfasta.idx")
            self.run_python(
                f"SeqIO.index_db({filename!r}, {FASTA!r}, 'starfasta').close()\n"
            )
            self.run_python(
                f"""\
index = SeqIO.index_db({filename!r})
assert len(index) == 3
for record in index.values():
    assert record.annotations["parsed by"] == "fakeplugin"
index.close()
"""
            )


class AlignPlugin(PluginTestCase):
    """Bio.Align uses a format from an entry point, case-insensitively."""

    def test_star_aln(self):
        warnings = self.run_python(
            f"""\
expected = Align.read({CLUSTAL!r}, "clustal").format("clustal")
for name in ["star-aln", "Star-Aln"]:
    alignment = Align.read({CLUSTAL!r}, name)
    assert alignment.format(name) == expected, name
    assert format(alignment, name) == expected, name
    assert len(list(Align.parse({CLUSTAL!r}, name))) == 1, name
    handle = StringIO()
    assert Align.write(alignment, handle, name) == 1, name

import fakeplugin_aln

assert Align._registry["star-aln"] is fakeplugin_aln
assert "STAR-ALN" not in dict.keys(Align._registry)
assert Align._registry.plugins == {{"star-aln"}}
assert "star-aln" not in Align.formats
assert "goodaln" not in Align._registry
"""
        )
        self.assertEqual(warnings, ALIGN_WARNINGS)


class Precedence(PluginTestCase):
    """register_format wins over an entry point, whichever comes first."""

    def test_register_then_miss(self):
        warnings = self.run_python(
            """\
from Bio.Align import clustal
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.FastaIO import FastaWriter

SeqIO.register_format("starfasta", FastaIterator, FastaWriter)
Align.register_format("STAR-ALN", clustal)
for table in [SeqIO._FormatToIterator, SeqIO._FormatToWriter, Align._registry]:
    assert "absent" not in table
    assert "starfasta" not in table.plugins
    assert "star-aln" not in table.plugins
assert SeqIO._FormatToIterator["starfasta"] is FastaIterator
assert SeqIO._FormatToWriter["starfasta"] is FastaWriter
assert Align._registry["star-aln"] is clustal
"""
        )
        self.assertEqual(warnings, ALL_WARNINGS)

    def test_miss_then_register(self):
        warnings = self.run_python(
            """\
from Bio.Align import clustal
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.FastaIO import FastaTwoLineIterator
from Bio.SeqIO.FastaIO import FastaWriter

for table in [SeqIO._FormatToIterator, SeqIO._FormatToWriter, Align._registry]:
    assert "absent" not in table
assert "starfasta" in SeqIO._FormatToIterator.plugins
assert "starfasta" in SeqIO._FormatToWriter.plugins
assert "star-aln" in Align._registry.plugins

SeqIO.register_format("starfasta", FastaIterator, FastaWriter)
Align.register_format("STAR-ALN", clustal)
# Replacing the plugins did not import them:
assert "fakeplugin_formats" not in sys.modules
assert "fakeplugin_aln" not in sys.modules
assert SeqIO._FormatToIterator["starfasta"] is FastaIterator
assert SeqIO._FormatToWriter["starfasta"] is FastaWriter
assert Align._registry["star-aln"] is clustal
for table in [SeqIO._FormatToIterator, SeqIO._FormatToWriter, Align._registry]:
    assert not table.plugins & {"starfasta", "star-aln"}

# They are explicit registrations now:
try:
    SeqIO.register_format("starfasta", FastaTwoLineIterator)
except ValueError:
    pass
else:
    raise AssertionError("replaced a registration without replace=True")
"""
        )
        self.assertEqual(warnings, ALL_WARNINGS)


class Conflicts(PluginTestCase):
    """Entry points that cannot be used are skipped, with one warning each."""

    def test_skipped_entry_points(self):
        warnings = self.run_python(
            """\
from Bio.SeqIO.FastaIO import FastaIterator

table = SeqIO._FormatToIterator
assert "absent" not in table
assert table["fasta"] is FastaIterator
for name in ["BadName", "badname", "dup", "NoName", "weird"]:
    assert name not in table, name
    assert name not in table.plugins, name
# A distribution with no name, or with one bad entry point, still counts:
assert table.plugins >= {"starfasta", "same", "broken", "noname"}
assert dict.get(table, "same") == "fakeplugin_formats:StarFastaIterator"
"""
        )
        self.assertEqual(warnings, ITERATOR_WARNINGS)

    def test_plugin_failing_to_import(self):
        self.run_python(
            f"""\
for function in [SeqIO.parse, SeqIO.read]:
    try:
        function({FASTA!r}, "broken")
    except ImportError as error:
        assert str(error) == "fakeplugin_broken needs a library", error
    else:
        raise AssertionError("no ImportError")
assert len(list(SeqIO.parse({FASTA!r}, "starfasta"))) == 3
assert len(list(SeqIO.parse({FASTA!r}, "fasta"))) == 3
"""
        )

    def test_scan_failing(self):
        """A failing scan warns once, and means no plugins."""
        warnings = self.run_python(
            f"""\
import importlib.metadata


def entry_points(**params):
    raise RuntimeError("unreadable metadata")


importlib.metadata.entry_points = entry_points
for table in [SeqIO._FormatToIterator, SeqIO._FormatToWriter, Align._registry]:
    assert "starfasta" not in table
    assert "star-aln" not in table
    assert not table.plugins
assert len(list(SeqIO.parse({FASTA!r}, "fasta"))) == 3
"""
        )
        self.assertEqual(
            warnings,
            [
                "Could not look for file format plugins, so none are used:"
                " RuntimeError('unreadable metadata')"
            ],
        )


class Listing(PluginTestCase):
    """Listing a table's names includes the plugins, before any miss."""

    def test_first_listing(self):
        for listing in ["list(table)", "table.keys()", "len(table)"]:
            with self.subTest(listing=listing):
                warnings = self.run_python(
                    f"""\
table = SeqIO._FormatToIterator
listed = {listing}
assert "starfasta" in table.plugins
if isinstance(listed, int):
    assert listed == len(table.builtin) + len(table.plugins), listed
else:
    assert "starfasta" in listed
"""
                )
                self.assertEqual(warnings, ITERATOR_WARNINGS)

    def test_miss_inside_loop(self):
        """A miss inside a loop over the same table adds no names to it."""
        self.run_python(
            """\
table = SeqIO._FormatToWriter
names = []
for name in table:
    names.append(name)
    assert table.get(name + "-absent") is None
    assert name + "-absent" not in table
assert "starfasta" in names
"""
        )


class Threads(PluginTestCase):
    """Threads missing at once scan once, and warn once."""

    def test_eight_threads(self):
        warnings = self.run_python(
            """\
import importlib.metadata
import threading
import time

calls = []
entry_points = importlib.metadata.entry_points


def counting_entry_points(**params):
    calls.append(params)
    time.sleep(0.2)  # so that the other threads miss meanwhile
    return entry_points(**params)


importlib.metadata.entry_points = counting_entry_points
barrier = threading.Barrier(8)
found = []


def miss():
    barrier.wait()
    found.append(
        SeqIO._FormatToIterator.get("starfasta") is not None
        and SeqIO._FormatToWriter.get("starfasta") is not None
        and Align._registry.get("star-aln") is not None
    )


threads = [threading.Thread(target=miss) for i in range(8)]
for thread in threads:
    thread.start()
for thread in threads:
    thread.join()
# One scan served all three tables:
assert calls == [{}], calls
assert found == [True] * 8, found
"""
        )
        self.assertEqual(warnings, ALL_WARNINGS)

    def test_added_meanwhile(self):
        """A miss looks again if another thread added the plugins meanwhile."""
        warnings = self.run_python(
            """\
table = SeqIO._FormatToIterator
discover = table._discover


def discover_after_another_thread():
    # Another thread adds the plugins between this thread's first look and
    # its own call, which then has nothing left to do:
    discover()
    return discover()


table._discover = discover_after_another_thread
assert "starfasta" in table
"""
        )
        self.assertEqual(warnings, ITERATOR_WARNINGS)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
