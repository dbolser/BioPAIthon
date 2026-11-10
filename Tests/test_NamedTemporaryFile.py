# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests that parsers accept tempfile.NamedTemporaryFile handles.

On CPython the object returned by ``tempfile.NamedTemporaryFile`` is a wrapper
around the real file. Iterating over it works, but it has no ``__next__``
method, so calling ``next()`` on it raises ``TypeError: '_TemporaryFileWrapper'
object is not an iterator``. Parsers that called ``next(handle)`` failed on
such handles although the same file opened with ``open()`` parsed fine; see
Biopython issue #4926.

Each test here writes a file from the test suite into a NamedTemporaryFile,
parses it through the wrapper, and checks that the result is the same as
parsing the file opened with ``open()``.
"""

import unittest
from io import StringIO
from tempfile import NamedTemporaryFile

import support

from Bio import Compass
from Bio.Affy import CelFile
from Bio import motifs
from Bio import SeqIO
from Bio import SwissProt
from Bio.Emboss import Primer3
from Bio.PopGen import GenePop
from Bio.Sequencing import Ace
from Bio.Sequencing import Phd


class NamedTemporaryFileTests(unittest.TestCase):
    """Parse each file both from open() and from a NamedTemporaryFile."""

    def check(self, path, parse):
        """Return parse() of the file at path, checking both kinds of handle.

        The file is parsed once from ``open()`` and once from a
        NamedTemporaryFile holding a copy of it; both must give the same
        summary, which is returned so the caller can check its content.
        """
        path = support.DATA / path
        with open(path) as handle:
            expected = parse(handle)
            handle.seek(0)
            text = handle.read()
        with NamedTemporaryFile("w+t") as handle:
            handle.write(text)
            handle.seek(0)
            self.assertEqual(parse(handle), expected)
        return expected

    def test_SeqIO_ig(self):
        ids = self.check(
            "IntelliGenetics/TAT_mase_nuc.txt",
            lambda h: [record.id for record in SeqIO.parse(h, "ig")],
        )
        self.assertEqual(len(ids), 17)
        self.assertEqual(ids[0], "A_U455")

    def test_Compass_read(self):
        query = self.check("Compass/comtest1", lambda h: Compass.read(h).query)
        self.assertEqual(query, "60456.blo.gz.aln")

    def test_Compass_parse(self):
        hits = self.check(
            "Compass/comtest2", lambda h: [r.hit for r in Compass.parse(h)]
        )
        self.assertEqual(len(hits), 4)

    def test_CelFile_read_v3(self):
        shape = self.check(
            "Affy/affy_v3_example.CEL", lambda h: CelFile.read(h).intensities.shape
        )
        self.assertEqual(shape, (5, 5))

    def test_motifs_alignace(self):
        consensus = self.check(
            "motifs/alignace.out",
            lambda h: [str(m.consensus) for m in motifs.parse(h, "alignace")],
        )
        self.assertEqual(len(consensus), 16)

    def test_motifs_sites(self):
        consensus = self.check(
            "motifs/Arnt.sites", lambda h: str(motifs.read(h, "sites").consensus)
        )
        self.assertEqual(consensus, "CACGTG")

    def test_Ace_read(self):
        summary = self.check(
            "Ace/contig1.ace",
            lambda h: [(c.name, c.nreads) for c in Ace.read(h).contigs],
        )
        self.assertEqual(summary, [("Contig1", 2), ("Contig2", 14)])

    def test_Ace_parse(self):
        names = self.check(
            "Ace/consed_sample.ace", lambda h: [c.name for c in Ace.parse(h)]
        )
        self.assertEqual(names, ["Contig1"])

    def test_Phd_read(self):
        name = self.check("Phd/phd2", lambda h: Phd.read(h).file_name)
        self.assertEqual(name, "ML4924R")

    def test_Phd_parse(self):
        names = self.check("Phd/phd1", lambda h: [r.file_name for r in Phd.parse(h)])
        self.assertEqual(len(names), 3)

    def test_SwissProt_read(self):
        name = self.check(
            "SwissProt/F2CXE6.txt", lambda h: SwissProt.read(h).entry_name
        )
        self.assertEqual(name, "F2CXE6_HORVD")

    def test_SwissProt_parse(self):
        names = self.check(
            "SwissProt/multi_ex.txt",
            lambda h: [r.entry_name for r in SwissProt.parse(h)],
        )
        self.assertEqual(len(names), 8)

    def test_GenePop_read(self):
        loci = self.check("PopGen/c2line.gen", lambda h: GenePop.read(h).loci_list)
        self.assertEqual(loci, ["136255903", "136257048", "136257636"])

    def test_Primer3_read(self):
        count = self.check(
            "Emboss/bac_find.primer3", lambda h: len(Primer3.read(h).primers)
        )
        self.assertEqual(count, 5)

    def test_Primer3_parse(self):
        counts = self.check(
            "Emboss/NirK_full.primer3",
            lambda h: [len(r.primers) for r in Primer3.parse(h)],
        )
        self.assertEqual(counts, [5] * 16)


class EndOfFileTests(unittest.TestCase):
    """The end-of-file checks rewritten to use readline still work."""

    def test_Ace_read_empty(self):
        with self.assertRaisesRegex(ValueError, "Premature end of file"):
            Ace.read(StringIO(""))

    def test_SwissProt_read_more_than_one(self):
        with open(support.DATA / "SwissProt" / "multi_ex.txt") as handle:
            text = handle.read()
        with NamedTemporaryFile("w+t") as handle:
            handle.write(text)
            handle.seek(0)
            with self.assertRaisesRegex(ValueError, "More than one SwissProt"):
                SwissProt.read(handle)

    def test_SwissProt_read_empty(self):
        with self.assertRaisesRegex(ValueError, "No SwissProt record found"):
            SwissProt.read(StringIO(""))

    def test_Phd_read_more_than_one(self):
        with open(support.DATA / "Phd" / "phd1") as handle:
            with self.assertRaisesRegex(ValueError, "More than one PHD record"):
                Phd.read(handle)

    def test_CelFile_v3_truncated(self):
        # Used to leak a bare StopIteration rather than raise ValueError.
        for text in ("[CEL", "[CEL]\n"):
            with self.assertRaisesRegex(ValueError, "Version 3 CEL"):
                CelFile.read(StringIO(text))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
