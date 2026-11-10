# This code is part of the biop-ai distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Offline tests of the URLs Bio.KEGG.REST builds, and of what it rejects.

urlopen is replaced by a fake, so no request reaches rest.kegg.jp; the
online tests are in test_KEGG_online.py. The URLs accepted here were each
checked against the live server.
"""

import io
import time
import types
import unittest
from unittest import mock

from Bio.KEGG import REST
from Bio.KEGG.REST import kegg_conv
from Bio.KEGG.REST import kegg_find
from Bio.KEGG.REST import kegg_get
from Bio.KEGG.REST import kegg_info
from Bio.KEGG.REST import kegg_link
from Bio.KEGG.REST import kegg_list

BASE = "https://rest.kegg.jp/"


def fake_urlopen(url):
    """Return an empty response that remembers its URL, as urlopen does."""
    resp = io.BytesIO(b"")
    resp.url = url
    return resp


class KEGGRESTURLs(unittest.TestCase):
    """Check the URL each call requests, without the network."""

    def setUp(self):
        previous = REST._q.previous
        self.addCleanup(setattr, REST._q, "previous", previous)
        no_wait = types.SimpleNamespace(time=time.time, sleep=lambda seconds: None)
        for patcher in (
            mock.patch.object(REST, "urlopen", side_effect=fake_urlopen),
            mock.patch.object(REST, "time", no_wait),
        ):
            patcher.start()
            self.addCleanup(patcher.stop)

    def assertURL(self, handle, path):
        """Check the handle came from the given path under rest.kegg.jp."""
        with handle:
            self.assertEqual(handle.url, BASE + path)

    def test_info(self):
        self.assertURL(kegg_info("kegg"), "info/kegg")
        self.assertURL(kegg_info("hsa"), "info/hsa")

    def test_list_database(self):
        self.assertURL(kegg_list("pathway"), "list/pathway")
        self.assertURL(kegg_list("brite"), "list/brite")
        self.assertURL(kegg_list("genome"), "list/genome")
        self.assertURL(kegg_list("hsa"), "list/hsa")
        self.assertURL(kegg_list("T01001"), "list/T01001")

    def test_list_pathway_org(self):
        self.assertURL(kegg_list("pathway", "hsa"), "list/pathway/hsa")

    def test_list_option(self):
        self.assertURL(kegg_list("brite", option="jp"), "list/brite/jp")
        self.assertURL(kegg_list("genome", option="animals"), "list/genome/animals")
        self.assertURL(kegg_list("genome", option="562"), "list/genome/562")

    def test_list_dbentries(self):
        self.assertURL(kegg_list("hsa:10458+ece:Z5100"), "list/hsa:10458+ece:Z5100")
        self.assertURL(
            kegg_list(["hsa:10458", "ece:Z5100"]), "list/hsa:10458+ece:Z5100"
        )
        self.assertURL(kegg_list(["C01290", "G00092"]), "list/C01290+G00092")

    def test_list_rejects(self):
        # KEGG answers /list/module/<org> with HTTP 400
        with self.assertRaises(ValueError):
            kegg_list("module", "hsa")
        with self.assertRaises(ValueError):
            kegg_list("compound", "hsa")
        with self.assertRaises(ValueError):
            kegg_list("pathway", option="jp")
        with self.assertRaises(ValueError):
            kegg_list("brite", "hsa", option="jp")
        # The subset is option, not org
        with self.assertRaisesRegex(ValueError, "option='animals'"):
            kegg_list("genome", "animals")

    def test_list_limit(self):
        entries = [f"hsa:{i}" for i in range(1, 101)]
        self.assertURL(kegg_list(entries), "list/" + "+".join(entries))
        self.assertURL(kegg_list("+".join(entries)), "list/" + "+".join(entries))
        # KEGG silently returns only the first 100, so both forms refuse more
        entries.append("hsa:101")
        with self.assertRaises(ValueError):
            kegg_list(entries)
        with self.assertRaises(ValueError):
            kegg_list("+".join(entries))

    def test_find(self):
        self.assertURL(kegg_find("genes", ["shiga", "toxin"]), "find/genes/shiga+toxin")
        self.assertURL(
            kegg_find("compound", "C7H10O5", "formula"),
            "find/compound/C7H10O5/formula",
        )
        self.assertURL(kegg_find("compound", "(+)", "nop"), "find/compound/(+)/nop")
        # nop is not a chemical field: KEGG takes it for any keyword search
        self.assertURL(kegg_find("genes", "shiga", "nop"), "find/genes/shiga/nop")
        self.assertURL(
            kegg_find("pathway", "glycolysis", "nop"), "find/pathway/glycolysis/nop"
        )
        self.assertURL(
            kegg_find("genes", ["shiga", "toxin"], "nop"),
            "find/genes/shiga+toxin/nop",
        )
        # but the chemical fields only for compound and drug (HTTP 400)
        with self.assertRaises(ValueError):
            kegg_find("genes", "shiga", "formula")

    def test_get(self):
        self.assertURL(
            kegg_get(["hsa:10458", "ece:Z5100"], "aaseq"),
            "get/hsa:10458+ece:Z5100/aaseq",
        )
        self.assertURL(kegg_get("hsa00600", "conf"), "get/hsa00600/conf")
        with self.assertRaises(ValueError):
            kegg_get("hsa00600", "png")

    def test_conv_database(self):
        self.assertURL(kegg_conv("eco", "ncbi-geneid"), "conv/eco/ncbi-geneid")
        self.assertURL(kegg_conv("ncbi-geneid", "eco"), "conv/ncbi-geneid/eco")
        self.assertURL(kegg_conv("hsa", "ncbi-proteinid"), "conv/hsa/ncbi-proteinid")
        self.assertURL(kegg_conv("compound", "pubchem"), "conv/compound/pubchem")
        self.assertURL(kegg_conv("compound", "chebi"), "conv/compound/chebi")
        self.assertURL(kegg_conv("chebi", "drug"), "conv/chebi/drug")

    def test_conv_dbentries(self):
        self.assertURL(
            kegg_conv("ncbi-proteinid", "hsa:10458+ece:Z5100"),
            "conv/ncbi-proteinid/hsa:10458+ece:Z5100",
        )
        self.assertURL(
            kegg_conv("ncbi-proteinid", ["hsa:10458", "ece:Z5100"]),
            "conv/ncbi-proteinid/hsa:10458+ece:Z5100",
        )
        self.assertURL(
            kegg_conv("genes", "ncbi-geneid:948364"),
            "conv/genes/ncbi-geneid:948364",
        )
        # upstream issue 729: an organism as the target of the entries form
        self.assertURL(
            kegg_conv("eco", "ncbi-geneid:948364"), "conv/eco/ncbi-geneid:948364"
        )
        self.assertURL(kegg_conv("hsa", "up:Q9Y6N5"), "conv/hsa/up:Q9Y6N5")
        self.assertURL(
            kegg_conv("compound", "pubchem:3303"), "conv/compound/pubchem:3303"
        )
        self.assertURL(kegg_conv("chebi", "C00001"), "conv/chebi/C00001")
        self.assertURL(
            kegg_conv("pubchem", ["C00001", "D00001", "G00001"]),
            "conv/pubchem/C00001+D00001+G00001",
        )

    def test_conv_rejects(self):
        with self.assertRaises(ValueError):
            kegg_conv("compound", "glycan")
        with self.assertRaises(ValueError):
            kegg_conv("eco", "hsa")
        # KEGG answers these with HTTP 400, or with nothing at all
        with self.assertRaises(ValueError):
            kegg_conv("ncbi-geneid", "eco", "turtle")
        with self.assertRaises(ValueError):
            kegg_conv("ncbi-geneid", "hsa:10458", "n-triple")
        with self.assertRaises(ValueError):
            kegg_conv("eco", "ncbi-gi")
        with self.assertRaises(ValueError):
            kegg_conv("ncbi-gi", "hsa:10458")
        with self.assertRaises(ValueError):
            kegg_conv("hsa", "ncbi-gi:4507949")
        with self.assertRaises(ValueError):
            kegg_conv("hsa", ["hsa:10458", "ncbi-gi:4507949"])
        # Gene and chemical databases do not mix (HTTP 400 or an empty result)
        with self.assertRaises(ValueError):
            kegg_conv("pubchem", "hsa")
        with self.assertRaises(ValueError):
            kegg_conv("chebi", "T01001")
        with self.assertRaises(ValueError):
            kegg_conv("uniprot", "compound")

    def test_link(self):
        self.assertURL(kegg_link("pathway", "hsa"), "link/pathway/hsa")
        self.assertURL(
            kegg_link("pathway", ["hsa:10458", "ece:Z5100"]),
            "link/pathway/hsa:10458+ece:Z5100",
        )
        self.assertURL(
            kegg_link("genome", "taxid:562", "species"),
            "link/genome/taxid:562/species",
        )
        self.assertURL(kegg_link("atc", "D01441", "turtle"), "link/atc/D01441/turtle")
        # dr is KEGG's abbreviation of drug
        self.assertURL(kegg_link("dr", "atc", "turtle"), "link/dr/atc/turtle")
        # KEGG offers RDF output of link only for drug, atc and jtc
        with self.assertRaises(ValueError):
            kegg_link("pathway", "hsa:10458", "turtle")
        with self.assertRaisesRegex(ValueError, "kegg link request"):
            kegg_link("pathway", "hsa", "kingdom")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
