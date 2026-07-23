# Copyright 2019-2020 by Sergio Valqui. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.UniProt.GOA Module.

GOA files can be found here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
"""

import io
import os
import tempfile
import unittest

from Bio.UniProt import GOA


class GoaTests(unittest.TestCase):
    """Test for UniProt GOA Files."""

    def test_gaf_iterator(self):
        """Test GOA GAF file iterator."""
        # Test GAF 2.0
        recs = []
        with open("UniProt/goa_yeast.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)

        # Check number of records
        self.assertEqual(len(recs), 587)
        # Check keys are same as predefined fields
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GAF20FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A023PXA5")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "YAL019W-A")
        self.assertEqual(recs[0]["Qualifier"], [""])
        self.assertEqual(recs[0]["GO_ID"], "GO:0003674")
        self.assertEqual(recs[0]["DB:Reference"], ["GO_REF:0000015"])
        self.assertEqual(recs[0]["Evidence"], "ND")
        self.assertEqual(recs[0]["With"], [""])

        # Test GAF 2.1, it has the same fields as GAF 2.0
        recs = []
        with open("UniProt/gene_association.goa_yeast.1.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)

        # Check number of records
        self.assertEqual(len(recs), 300)
        # Check keys are same as predefined fields
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GAF20FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "P17536")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "TPM1")
        self.assertEqual(recs[0]["Qualifier"], [""])
        self.assertEqual(recs[0]["GO_ID"], "GO:0000001")
        self.assertEqual(recs[0]["DB:Reference"], ["PMID:10652251"])
        self.assertEqual(recs[0]["Evidence"], "TAS")
        self.assertEqual(recs[0]["With"], [""])

    def test_gpa_iterator(self):
        """Test GOA GPA file iterator."""
        recs = []
        with open("UniProt/goa_yeast.gpa.59.gpa") as handle:
            for rec in GOA.gpa_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 300)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPA11FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A023PXA5")
        self.assertEqual(recs[0]["Qualifier"], ["enables"])
        self.assertEqual(recs[0]["GO_ID"], "GO:0003674")
        self.assertEqual(recs[0]["DB:Reference"], ["GO_REF:0000015"])
        self.assertEqual(recs[0]["ECO_Evidence_code"], "ECO:0000307")
        self.assertEqual(recs[0]["With"], [""])
        self.assertEqual(recs[0]["Interacting_taxon_ID"], "")
        self.assertEqual(recs[0]["Date"], "20030730")
        self.assertEqual(recs[0]["Assigned_by"], "SGD")
        self.assertEqual(recs[0]["Annotation Extension"], [""])
        self.assertEqual(recs[0]["Annotation_Properties"], "go_evidence=ND")

    def test_gpi_iterator(self):
        """Test GOA GPI file iterator, gpi-version: 1.1."""
        recs = []
        with open("UniProt/gp_information.goa_yeast.28.gpi") as handle:
            for rec in GOA.gpi_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 300)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPI11FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB_Object_ID"], "A2P2R3")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "YMR084W")
        self.assertEqual(
            recs[0]["DB_Object_Name"],
            ["Putative glutamine--fructose-6-phosphate aminotransferase [isomerizing]"],
        )
        self.assertEqual(recs[0]["DB_Object_Synonym"], ["YM084_YEAST", "YMR084W"])
        self.assertEqual(recs[0]["DB_Object_Type"], "protein")
        self.assertEqual(recs[0]["Taxon"], "taxon:559292")
        self.assertEqual(recs[0]["Parent_Object_ID"], "")
        self.assertEqual(recs[0]["DB_Xref"], [""])
        self.assertEqual(recs[0]["Gene_Product_Properties"], ["db_subset=Swiss-Prot"])

    def test_gpi_iterator_one_two(self):
        """Test GOA GPI file iterator, gpi-version: 1.2."""
        recs = []
        with open("UniProt/goa_human_sample.gpi") as handle:
            for rec in GOA.gpi_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 9)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPI12FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A024R1R8")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "hCG_2014768")
        self.assertEqual(recs[0]["DB_Object_Name"], ["HCG2014768, isoform CRA_a"])
        self.assertEqual(recs[0]["DB_Object_Synonym"], ["hCG_2014768"])
        self.assertEqual(recs[0]["DB_Object_Type"], "protein")
        self.assertEqual(recs[0]["Taxon"], "taxon:9606")
        self.assertEqual(recs[0]["Parent_Object_ID"], "")
        self.assertEqual(recs[0]["DB_Xref"], [""])
        self.assertEqual(recs[0]["Gene_Product_Properties"], ["db_subset=TrEMBL"])

    def test_selection_writing(self):
        """Test record_has, and writerec.

        Adapted from Bio.UniProt.GOA.py by Iddo Friedberg idoerg@gmail.com.
        """
        recs = []
        filtered = []

        # Fields to filter
        evidence = {"Evidence": {"ND"}}
        synonym = {"Synonym": {"YA19A_YEAST", "YAL019W-A"}}
        taxon_id = {"Taxon_ID": {"taxon:559292"}}

        # Temporal file to test writerec
        f_number, f_filtered = tempfile.mkstemp()
        os.close(f_number)

        # Open a file and select records as per filter
        with open("UniProt/goa_yeast.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)
                # Filtering
                if (
                    GOA.record_has(rec, taxon_id)
                    and GOA.record_has(rec, evidence)
                    and GOA.record_has(rec, synonym)
                ):
                    filtered.append(rec)

        # Check number of filtered records
        self.assertEqual(len(filtered), 3)

        # Write the filtered records to a file using writerec
        with open(f_filtered, "w") as handle:
            # '!gaf-version: 2.1'
            handle.write("!gaf-version: 2.1 \n")  # Adding file header
            for rec in filtered:
                GOA.writerec(rec, handle)

        # Open and read the file containing the filtered records
        recs_ff = []  # Records from filtered file
        with open(f_filtered) as handle:
            for rec in GOA.gafiterator(handle):
                recs_ff.append(rec)

        # Delete test file
        os.remove(f_filtered)

        # Compare, recs saved by writerec and filtered recs
        self.assertEqual(filtered, recs_ff)


def gaf10_line(object_id, go_id, symbol="ABC1", synonym="ABC1|YAL001C"):
    """Return one GAF 1.0 data line (15 tab separated fields)."""
    return "\t".join(
        [
            "UniProtKB",  # DB
            object_id,  # DB_Object_ID
            symbol,  # DB_Object_Symbol
            "NOT|contributes_to",  # Qualifier
            go_id,  # GO_ID
            "PMID:12345|GO_REF:0000002",  # DB:Reference
            "IDA",  # Evidence
            "UniProtKB:P12345|InterPro:IPR000001",  # With
            "F",  # Aspect
            "Example protein",  # DB_Object_Name
            synonym,  # Synonym
            "protein",  # DB_Object_Type
            "taxon:559292|taxon:4932",  # Taxon_ID
            "20140101",  # Date
            "SGD",  # Assigned_By
        ]
    )


def gaf20_line(object_id, go_id, symbol="ABC1", synonym="ABC1|YAL001C"):
    """Return one GAF 2.0 data line (17 tab separated fields)."""
    return gaf10_line(object_id, go_id, symbol, synonym) + "\t\t"


def gpa10_line(object_id, go_id):
    """Return one GPA 1.0 data line (12 tab separated fields)."""
    return "\t".join(
        [
            "UniProtKB",  # DB
            object_id,  # DB_Object_ID
            "enables|NOT",  # Qualifier
            go_id,  # GO_ID
            "PMID:12345|GO_REF:0000002",  # DB:Reference
            "IDA",  # Evidence code
            "UniProtKB:P12345|InterPro:IPR000001",  # With
            "taxon:4932",  # Interacting_taxon_ID
            "20140101",  # Date
            "SGD",  # Assigned_by
            "part_of(GO:0005634)|occurs_in(CL:0000001)",  # Annotation_Extension
            "P12345-1",  # Spliceform_ID
        ]
    )


def gpi10_line(object_id):
    """Return one GPI 1.0 data line (11 tab separated fields)."""
    return "\t".join(
        [
            "UniProtKB",  # DB
            "Swiss-Prot",  # DB_subset
            object_id,  # DB_Object_ID
            "ABC1",  # DB_Object_Symbol
            "Example protein",  # DB_Object_Name
            "ABC1_YEAST|YAL001C",  # DB_Object_Synonym
            "protein",  # DB_Object_Type
            "taxon:559292",  # Taxon
            "PRO|Complex Portal",  # Annotation_Target_Set
            "20140101",  # Annotation_Completed
            "",  # Parent_Object_ID
        ]
    )


def handle(*lines):
    """Return a file like object holding the given lines."""
    return io.StringIO("\n".join(lines) + "\n")


class GafOneZeroTests(unittest.TestCase):
    """Tests for the GAF 1.0 iterator, which has only 15 fields."""

    def test_gaf10_iterator(self):
        """GAF 1.0 records use the 15 GAF 1.0 field names."""
        stream = handle(
            "!gaf-version: 1.0",
            "!comment",
            gaf10_line("P12345", "GO:0003674"),
            "",  # a blank line is skipped
            gaf10_line("P67890", "GO:0005575"),
        )
        records = list(GOA.gafiterator(stream))
        self.assertEqual(len(records), 2)
        self.assertEqual(sorted(records[0]), sorted(GOA.GAF10FIELDS))
        # the 15 field format has no Annotation_Extension
        self.assertNotIn("Annotation_Extension", records[0])
        self.assertEqual(records[0]["DB_Object_ID"], "P12345")
        self.assertEqual(records[0]["GO_ID"], "GO:0003674")
        self.assertEqual(records[0]["Assigned_By"], "SGD")
        # the pipe separated fields are split into lists
        self.assertEqual(records[0]["Qualifier"], ["NOT", "contributes_to"])
        self.assertEqual(records[0]["DB:Reference"], ["PMID:12345", "GO_REF:0000002"])
        self.assertEqual(records[0]["With"], ["UniProtKB:P12345", "InterPro:IPR000001"])
        self.assertEqual(records[0]["Synonym"], ["ABC1", "YAL001C"])
        self.assertEqual(records[0]["Taxon_ID"], ["taxon:559292", "taxon:4932"])
        # ...and the others are not
        self.assertEqual(records[0]["Evidence"], "IDA")
        self.assertEqual(records[1]["DB_Object_ID"], "P67890")

    def test_all_gaf_versions_are_recognised(self):
        """GAF 2.0, 2.1 and 2.2 all use the 17 field parser."""
        for version in ("2.0", "2.1", "2.2"):
            with self.subTest(version=version):
                stream = handle(
                    f"!gaf-version: {version}",
                    "",  # a blank line is skipped
                    gaf20_line("P12345", "GO:0003674"),
                )
                records = list(GOA.gafiterator(stream))
                self.assertEqual(len(records), 1)
                self.assertEqual(sorted(records[0]), sorted(GOA.GAF20FIELDS))

    def test_unknown_gaf_version(self):
        """An unrecognised GAF version is rejected."""
        stream = handle("!gaf-version: 9.9", gaf20_line("P12345", "GO:0003674"))
        with self.assertRaises(ValueError) as context:
            GOA.gafiterator(stream)
        self.assertIn("Unknown GAF version", str(context.exception))


class GafByProteinTests(unittest.TestCase):
    """Tests for the GAF iterators that group records by protein."""

    def records(self, version, line_maker):
        """Return the grouped records of a three protein file."""
        stream = handle(
            f"!gaf-version: {version}",
            "!comment",
            "",  # a blank line is skipped
            line_maker("P00001", "GO:0000001"),
            line_maker("P00001", "GO:0000002"),
            line_maker("P00002", "GO:0000003"),
            line_maker("P00003", "GO:0000004"),
            line_maker("P00003", "GO:0000005"),
        )
        return list(GOA.gafbyproteiniterator(stream))

    def test_consecutive_records_are_grouped(self):
        """Consecutive records for one protein come back as one list."""
        for version, line_maker in (
            ("1.0", gaf10_line),
            ("2.0", gaf20_line),
            ("2.1", gaf20_line),
            ("2.2", gaf20_line),
        ):
            with self.subTest(version=version):
                groups = self.records(version, line_maker)
                self.assertEqual(
                    [[record["DB_Object_ID"] for record in group] for group in groups],
                    [["P00001", "P00001"], ["P00002"]],
                )
                self.assertEqual(
                    [record["GO_ID"] for record in groups[0]],
                    ["GO:0000001", "GO:0000002"],
                )

    @unittest.expectedFailure
    def test_final_group_is_yielded(self):
        """The records of the last protein in the file should be returned.

        This currently fails: neither _gaf10byproteiniterator nor
        _gaf20byproteiniterator flushes the accumulated list when the input
        is exhausted, so the annotations of the last protein in the file are
        silently dropped.  Remove the expectedFailure decorator once that is
        fixed.
        """
        groups = self.records("2.1", gaf20_line)
        self.assertEqual(
            [[record["DB_Object_ID"] for record in group] for group in groups],
            [["P00001", "P00001"], ["P00002"], ["P00003", "P00003"]],
        )

    def test_unknown_gaf_version(self):
        """An unrecognised GAF version is rejected."""
        stream = handle("!gaf-version: 9.9", gaf20_line("P12345", "GO:0003674"))
        with self.assertRaises(ValueError) as context:
            GOA.gafbyproteiniterator(stream)
        self.assertIn("Unknown GAF version", str(context.exception))

    def test_writebyproteinrec_round_trip(self):
        """A group of records survives being written and read back."""
        groups = self.records("2.1", gaf20_line)
        f_number, filename = tempfile.mkstemp()
        os.close(f_number)
        try:
            with open(filename, "w") as output:
                output.write("!gaf-version: 2.1\n")
                for group in groups:
                    GOA.writebyproteinrec(group, output)
            with open(filename) as input_handle:
                read_back = list(GOA.gafbyproteiniterator(input_handle))
        finally:
            os.remove(filename)
        # the last group is lost on the way back in, see test_final_group
        self.assertEqual(read_back, groups[:-1])
        self.assertEqual(
            [record["DB_Object_ID"] for record in read_back[0]],
            ["P00001", "P00001"],
        )


class GpaTests(unittest.TestCase):
    """Tests for the GPA iterators."""

    def test_gpa10_iterator(self):
        """GPA 1.0 records use the GPA 1.0 field names."""
        stream = handle(
            "!gpa-version: 1.0",
            "!comment",
            gpa10_line("P12345", "GO:0003674"),
            "",  # a blank line is skipped
            gpa10_line("P67890", "GO:0005575"),
        )
        records = list(GOA.gpa_iterator(stream))
        self.assertEqual(len(records), 2)
        self.assertEqual(sorted(records[0]), sorted(GOA.GPA10FIELDS))
        # GPA 1.0 has a Spliceform_ID where 1.1 has Annotation_Properties
        self.assertEqual(records[0]["Spliceform_ID"], "P12345-1")
        self.assertNotIn("Annotation_Properties", records[0])
        self.assertEqual(records[0]["Qualifier"], ["enables", "NOT"])
        self.assertEqual(records[0]["DB:Reference"], ["PMID:12345", "GO_REF:0000002"])
        self.assertEqual(records[0]["With"], ["UniProtKB:P12345", "InterPro:IPR000001"])
        self.assertEqual(
            records[0]["Annotation_Extension"],
            ["part_of(GO:0005634)", "occurs_in(CL:0000001)"],
        )
        self.assertEqual(records[0]["Evidence code"], "IDA")
        self.assertEqual(records[1]["DB_Object_ID"], "P67890")

    def test_unknown_gpa_version(self):
        """An unrecognised GPA version is rejected."""
        stream = handle("!gpa-version: 9.9", gpa10_line("P12345", "GO:0003674"))
        with self.assertRaises(ValueError) as context:
            GOA.gpa_iterator(stream)
        self.assertIn("Unknown GPA version", str(context.exception))


class GpiTests(unittest.TestCase):
    """Tests for the GPI iterators."""

    def test_gpi10_iterator(self):
        """GPI 1.0 records use the GPI 1.0 field names."""
        stream = handle(
            "!gpi-version: 1.0",
            "!comment",
            gpi10_line("P12345"),
            "",  # a blank line is skipped
            gpi10_line("P67890"),
        )
        records = list(GOA.gpi_iterator(stream))
        self.assertEqual(len(records), 2)
        self.assertEqual(sorted(records[0]), sorted(GOA.GPI10FIELDS))
        # DB_subset and Annotation_Target_Set only exist in GPI 1.0
        self.assertEqual(records[0]["DB_subset"], "Swiss-Prot")
        self.assertEqual(records[0]["Annotation_Target_Set"], ["PRO", "Complex Portal"])
        self.assertEqual(records[0]["DB_Object_Synonym"], ["ABC1_YEAST", "YAL001C"])
        # the name is a plain string in 1.0, unlike in 1.1 and 1.2
        self.assertEqual(records[0]["DB_Object_Name"], "Example protein")
        self.assertEqual(records[1]["DB_Object_ID"], "P67890")

    def test_gpi_version_two_not_implemented(self):
        """GPI 2.1 files are recognised but not yet parsed."""
        stream = handle("!gpi-version: 2.1", gpi10_line("P12345"))
        with self.assertRaises(NotImplementedError):
            GOA.gpi_iterator(stream)

    def test_unknown_gpi_version(self):
        """An unrecognised GPI version is rejected."""
        stream = handle("!gpi-version: 9.9", gpi10_line("P12345"))
        with self.assertRaises(ValueError) as context:
            GOA.gpi_iterator(stream)
        self.assertIn("Unknown GPI version", str(context.exception))


class RecordHasTests(unittest.TestCase):
    """Tests for GOA.record_has."""

    record = {
        "DB_Object_ID": "P12345",
        "Evidence": "IDA",
        "Synonym": ["ABC1", "YAL001C"],
    }

    def test_matches_a_string_field(self):
        """A string field matches when its value is in the wanted set."""
        self.assertTrue(GOA.record_has(self.record, {"Evidence": {"IDA", "IPI"}}))
        self.assertFalse(GOA.record_has(self.record, {"Evidence": {"ND"}}))

    def test_matches_a_list_field(self):
        """A list field matches when any of its values is wanted."""
        self.assertTrue(GOA.record_has(self.record, {"Synonym": {"YAL001C"}}))
        self.assertFalse(GOA.record_has(self.record, {"Synonym": {"YAL002C"}}))

    def test_several_fields_are_combined_with_or(self):
        """Any one matching field is enough."""
        self.assertTrue(
            GOA.record_has(self.record, {"Evidence": {"ND"}, "Synonym": {"YAL001C"}})
        )
        self.assertFalse(
            GOA.record_has(self.record, {"Evidence": {"ND"}, "Synonym": {"YAL002C"}})
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
