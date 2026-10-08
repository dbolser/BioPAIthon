# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests that Bio.SearchIO rejects malformed input with ValueError.

These checks used to be assert statements, which raised a bare
AssertionError and vanished under ``python -O``, leaving the malformed
input to be misread silently.
"""

import os
import subprocess
import sys
import tempfile
import unittest
from io import StringIO

import support

import Bio
from Bio import SearchIO
from Bio.SearchIO import HSP
from Bio.SearchIO import HSPFragment
from Bio.SearchIO import Hit
from Bio.SearchIO import QueryResult


def edited(filename, old, new):
    """Return a test file's text with one occurrence of old replaced by new."""
    with open(filename) as handle:
        text = handle.read()
    # Guard against the test data changing underneath the test, which
    # would quietly make the test pass without exercising anything.
    if old not in text:
        raise RuntimeError(f"{old!r} not found in {filename}")
    return text.replace(old, new, 1)


def parse_all(text, fmt, **kwargs):
    """Parse the given text completely."""
    return list(SearchIO.parse(StringIO(text), fmt, **kwargs))


class BlastXmlIndex(unittest.TestCase):
    """The BLAST XML indexer needs one <Iteration> block per run of lines."""

    def index_text(self, text):
        handle, path = tempfile.mkstemp(suffix=".xml")
        self.addCleanup(os.remove, path)
        with os.fdopen(handle, "w") as handle:
            handle.write(text)
        index = SearchIO.index(path, "blast-xml")
        self.addCleanup(index.close)
        return index

    def test_iteration_not_at_line_start(self):
        text = edited(
            support.DATA / "Blast" / "xml_2226_blastp_001.xml",
            "    <Iteration>\n",
            "    <x/><Iteration>\n",
        )
        with self.assertRaises(ValueError) as cm:
            self.index_text(text)
        self.assertIn(
            "Expected <Iteration> at the start of the line", str(cm.exception)
        )

    def test_iteration_opened_on_closing_line(self):
        text = edited(
            support.DATA / "Blast" / "xml_2226_blastp_001.xml",
            "</Iteration>\n    <Iteration>",
            "</Iteration><Iteration>",
        )
        with self.assertRaises(ValueError) as cm:
            self.index_text(text)
        self.assertIn("Found <Iteration> before the </Iteration>", str(cm.exception))

    def test_truncated_iteration(self):
        with open(support.DATA / "Blast" / "xml_2226_blastp_001.xml") as handle:
            text = handle.read()
        text = text[: text.index("</Iteration>")]
        with self.assertRaises(ValueError) as cm:
            self.index_text(text)
        self.assertIn("File ended before the </Iteration>", str(cm.exception))

    def test_iteration_without_query_def(self):
        text = edited(
            support.DATA / "Blast" / "xml_2226_blastp_001.xml",
            "      <Iteration_query-def>random_s00</Iteration_query-def>\n",
            "",
        )
        with self.assertRaises(ValueError) as cm:
            self.index_text(text)
        self.assertIn("Expected <Iteration_query-def>", str(cm.exception))


class BlatPsl(unittest.TestCase):
    """PSL columns must agree with the blocks they summarise."""

    def test_qstart_disagrees_with_blocks(self):
        text = edited(
            support.DATA / "Blat" / "psl_34_001.psl",
            "hg18_dna\t33\t11\t27\tchr4",
            "hg18_dna\t33\t12\t27\tchr4",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "blat-psl")
        self.assertIn("has qStart 12, but its blocks give 11", str(cm.exception))

    def test_too_many_hit_starts(self):
        text = edited(
            support.DATA / "Blat" / "psl_34_001.psl",
            "1\t16,\t11,\t61646095,\n",
            "1\t16,\t11,\t61646095,61646200,\n",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "blat-psl")
        self.assertIn(
            "1 block sizes, 1 query starts and 2 hit starts", str(cm.exception)
        )


class Exonerate(unittest.TestCase):
    """Exonerate alignments, vulgar lines and cigar lines must agree."""

    def test_vulgar_score_disagrees_with_header(self):
        text = edited(
            support.DATA / "Exonerate" / "exn_22_q_multiple.exn",
            "560974 + 4485 M 897 897",
            "560974 + 4486 M 897 897",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "exonerate-vulgar")
        self.assertIn("Expected score '4485'", str(cm.exception))
        self.assertIn("found '4486'", str(cm.exception))

    def test_cigar_hit_id_disagrees_with_header(self):
        text = edited(
            support.DATA / "Exonerate" / "exn_22_q_multiple.exn",
            "cigar: gi|296142823|ref|NM_001178508.1| 0 897 + gi|330443482|",
            "cigar: gi|296142823|ref|NM_001178508.1| 0 897 + gi|330443483|",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "exonerate-cigar")
        self.assertIn("Expected hit ID", str(cm.exception))

    def test_unknown_vulgar_label(self):
        text = edited(
            support.DATA / "Exonerate" / "exn_22_m_est2genome_vulgar.exn",
            "6150 M 1230 1230",
            "6150 X 1230 1230",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "exonerate-vulgar")
        self.assertIn("Unexpected vulgar label 'X'", str(cm.exception))

    def test_text_header_without_query_range(self):
        text = edited(
            support.DATA / "Exonerate" / "exn_22_m_affine_local.exn",
            "   Query range: 0 -> 1230\n",
            "",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "exonerate-text")
        self.assertIn("No query_start in the alignment header", str(cm.exception))

    def test_codon_rows_parse_under_optimize(self):
        """The cdna2genome codon rows do not depend on assert statements.

        The parser used to detect the extra codon rows by catching the
        AssertionError from an assert, so under ``python -O`` it never
        added them and failed with IndexError.
        """
        # An absolute path, so the subprocess need not run in Tests/.
        exn = os.fspath(support.DATA / "Exonerate" / "exn_22_m_cdna2genome.exn")
        code = f"""\
from Bio import SearchIO
for qresult in SearchIO.parse({exn!r}, "exonerate-text"):
    for hit in qresult:
        for hsp in hit:
            print(hsp.query_range_all, hsp.hit_range_all)
            print([str(frag.hit.seq) for frag in hsp.fragments])
"""
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        outputs = []
        for flags in ([], ["-O"]):
            result = subprocess.run(
                [sys.executable, *flags, "-W", "ignore", "-c", code],
                capture_output=True,
                text=True,
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(result.stdout)
        self.assertTrue(outputs[0])
        self.assertEqual(outputs[0], outputs[1])


class FastaM10(unittest.TestCase):
    """The query and hit of a FASTA -m 10 alignment must be the same type."""

    def test_query_and_hit_types_differ(self):
        text = edited(
            support.DATA / "Fasta" / "output002.m10", "; sq_type: p", "; sq_type: D"
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "fasta-m10")
        self.assertIn("Query sequence type 'D' does not match", str(cm.exception))


class Hmmer(unittest.TestCase):
    """HMMER text and tables must have the layout the parsers expect."""

    def test_hmmer3_domain_numbered_out_of_order(self):
        text = edited(
            support.DATA / "Hmmer" / "text_30_hmmscan_001.out",
            "  == domain 2    score: -1.8 bits",
            "  == domain 3    score: -1.8 bits",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "hmmer3-text")
        self.assertIn("Expected '  == domain 2'", str(cm.exception))

    def test_hmmer3_domtab_row_too_short(self):
        with open(support.DATA / "Hmmer" / "domtab_30_hmmscan_001.out") as handle:
            lines = handle.readlines()
        lines[3] = " ".join(lines[3].split()[:20]) + "\n"
        with self.assertRaises(ValueError) as cm:
            parse_all("".join(lines), "hmmscan3-domtab")
        self.assertIn(
            "Expected at least 22 space-separated columns, found 20", str(cm.exception)
        )

    def test_hmmer2_without_program_line(self):
        text = edited(
            support.DATA / "Hmmer" / "text_22_hmmpfam_001.out",
            "hmmpfam - search one or more sequences against HMM database\n",
            "",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "hmmer2-text")
        self.assertIn("Expected a program line", str(cm.exception))


class Infernal(unittest.TestCase):
    """Infernal output must be consistent and have the expected columns."""

    def test_text_hit_table_row_too_short(self):
        text = edited(
            support.DATA / "Infernal" / "cmsearch_114_U2_Yeast_full.txt",
            "681747 - .. 0.91    no 0.33\n",
            "681747 - .. 0.91    no\n",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "infernal-text")
        self.assertIn(
            "Expected 16 columns in the hit table row, found 15", str(cm.exception)
        )

    def test_text_scores_table_row_without_description(self):
        # A target with no description gives a 12-column row, which is valid.
        text = edited(
            support.DATA / "Infernal" / "cmsearch_114_U2_Yeast_noali.txt",
            "681747 -  cm    no 0.33  TPA_inf: Saccharomyces cerevisiae S288C"
            " chromosome II,\n",
            "681747 -  cm    no 0.33\n",
        )
        (qresult,) = parse_all(text, "infernal-text")
        self.assertEqual(qresult[0].id, "ENA|BK006936|BK006936.2")
        self.assertEqual(qresult[0].description, "")
        self.assertEqual(qresult[0][0].gc, 0.33)

    def test_text_scores_table_row_too_short(self):
        text = edited(
            support.DATA / "Infernal" / "cmsearch_114_U2_Yeast_noali.txt",
            "681747 -  cm    no 0.33  TPA_inf: Saccharomyces cerevisiae S288C"
            " chromosome II,\n",
            "681747 -  cm    no\n",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "infernal-text")
        self.assertIn(
            "Expected at least 12 columns in the hit scores table row, found 11",
            str(cm.exception),
        )

    def test_tab_hit_descriptions_differ(self):
        text = edited(
            support.DATA / "Infernal" / "cmsearch_114_5S_Yeast.tbl",
            "489469      +    no    1 0.52   0.0   88.8   1.6e-18 !   TPA_inf:"
            " Saccharomyces cerevisiae S288C chromosome XII,",
            "489469      +    no    1 0.52   0.0   88.8   1.6e-18 !   TPA_inf:"
            " Saccharomyces cerevisiae S288C chromosome XIII,",
        )
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "infernal-tab")
        self.assertIn(
            "Hit 'ENA|BK006945|BK006945.2' has description", str(cm.exception)
        )

    def test_tab_unknown_format_number(self):
        with open(support.DATA / "Infernal" / "cmsearch_114_5S_Yeast.tbl") as handle:
            text = handle.read()
        with self.assertRaises(ValueError) as cm:
            parse_all(text, "infernal-tab", _fmt=4)
        self.assertIn("_fmt must be 1, 2 or 3, not 4", str(cm.exception))


class ModelArguments(unittest.TestCase):
    """The object model checks the values it is given."""

    def test_start_after_end(self):
        frag = HSPFragment("hit", "query")
        frag.hit_end = 5
        with self.assertRaises(ValueError) as cm:
            frag.hit_start = 10
        self.assertIn("Coordinate 10 must be <= hit_end 5", str(cm.exception))

    def test_end_before_start(self):
        frag = HSPFragment("hit", "query")
        frag.query_start = 10
        with self.assertRaises(ValueError):
            frag.query_end = 5

    def test_coordinate_not_an_int(self):
        frag = HSPFragment("hit", "query")
        with self.assertRaises(TypeError):
            frag.hit_start = 1.5

    def test_slice_with_short_annotation(self):
        frag = HSPFragment("hit", "query", hit="ACGT", query="ACGT")
        frag.aln_annotation["similarity"] = "||"
        with self.assertRaises(ValueError) as cm:
            frag[0:3]
        self.assertIn("'similarity' to have length 3", str(cm.exception))

    def test_absorb_hit_for_another_query(self):
        qresult = QueryResult(id="query1")
        hit = Hit([HSP([HSPFragment("hit", "query2")])])
        with self.assertRaises(ValueError) as cm:
            qresult.absorb(hit)
        self.assertIn("Expected Hit with query ID 'query1'", str(cm.exception))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
