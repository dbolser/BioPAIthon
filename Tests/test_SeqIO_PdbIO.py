# Copyright 2012 by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO PdbIO module."""

import unittest
import warnings
from io import StringIO

try:
    import numpy as np
    from numpy import dot  # Missing on PyPy's micronumpy

    del dot
    # We don't need this (?) but Bio.PDB imports it automatically :(
    from numpy.linalg import det  # Missing in PyPy 2.0 numpypy
    from numpy.linalg import svd  # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use PDB formats with SeqIO."
    ) from None

import support

from Bio import BiopythonDeprecationWarning
from Bio import BiopythonParserWarning
from Bio import SeqIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqIO.PdbIO import CifSeqresIterator


def SeqresTestGenerator(extension, parser):
    """Test factory for tests reading SEQRES (or similar) records.

    This is a factory returning a parameterised superclass for tests reading
    sequences from the sequence records of structure files.

    Arguments:
        extension:
            The extension of the files to read from the ``PDB`` directory (e.g.
            ``pdb`` or ``cif``).
        parser:
            The name of the SeqIO parser to use (e.g. ``pdb-atom``).

    """

    class SeqresTests(unittest.TestCase):
        """Use "parser" to parse sequence records from a structure file.

        Args:
            parser (str): Name of the parser used by SeqIO.
            extension (str): Extension of the files to parse.

        """

        def test_seqres_parse(self):
            """Parse a multi-chain PDB by SEQRES entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
            """
            chains = list(
                SeqIO.parse(support.DATA / "PDB" / f"2BEG.{extension}", parser)
            )
            self.assertEqual(len(chains), 5)
            actual_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
            for chain, chn_id in zip(chains, "ABCDE"):
                self.assertEqual(chain.id, "2BEG:" + chn_id)
                self.assertEqual(chain.annotations["chain"], chn_id)
                self.assertEqual(chain.seq, actual_seq)

        def test_seqres_read(self):
            """Read a single-chain structure by sequence entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
            """
            chain = SeqIO.read(support.DATA / "PDB" / f"1A8O.{extension}", parser)
            self.assertEqual(chain.id, "1A8O:A")
            self.assertEqual(chain.annotations["chain"], "A")
            self.assertEqual(
                chain.seq,
                "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPD"
                "CKTILKALGPGATLEEMMTACQG",
            )

        def test_seqres_missing(self):
            """Parse a PDB with no SEQRES entries."""
            chains = list(
                SeqIO.parse(support.DATA / "PDB" / f"a_structure.{extension}", parser)
            )
            self.assertEqual(len(chains), 0)

    return SeqresTests


class TestPdbSeqres(SeqresTestGenerator("pdb", "pdb-seqres")):
    """Test pdb-seqres SeqIO driver."""


class TestCifSeqres(SeqresTestGenerator("cif", "cif-seqres")):
    """Test cif-seqres SeqIO driver."""


class TestCifSeqresChainIds(unittest.TestCase):
    """Test the choice between label and author chain ids in cif-seqres."""

    def chains(self, records):
        return [record.annotations["chain"] for record in records]

    def test_default_warns_when_ids_differ(self):
        """Keep the label ids by default, but warn if they are not author ids."""
        # 4ZHL: label chain A is author chain U, label chain B is author chain P
        with self.assertWarns(BiopythonDeprecationWarning):
            records = list(SeqIO.parse(support.DATA / "PDB" / "4ZHL.cif", "cif-seqres"))
        self.assertEqual(self.chains(records), ["A", "B"])

    def test_warning_names_caller(self):
        """Report the warning at the caller's line, however it was called."""
        path = support.DATA / "PDB" / "1A7G.cif"
        with self.assertWarns(BiopythonDeprecationWarning) as direct:
            CifSeqresIterator(path)
        with self.assertWarns(BiopythonDeprecationWarning) as parsed:
            SeqIO.parse(path, "cif-seqres")
        with self.assertWarns(BiopythonDeprecationWarning) as read:
            SeqIO.read(path, "cif-seqres")
        for context in (direct, parsed, read):
            self.assertEqual(context.filename, __file__)

    def test_default_silent_when_ids_agree(self):
        """Do not warn if the label and author ids are the same."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonDeprecationWarning)
            records = list(SeqIO.parse(support.DATA / "PDB" / "2BEG.cif", "cif-seqres"))
        self.assertEqual(self.chains(records), ["A", "B", "C", "D", "E"])

    def test_auth_chains(self):
        """Name chains by author id with auth_chains=True, like cif-atom."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonDeprecationWarning)
            records = list(
                CifSeqresIterator(support.DATA / "PDB" / "4ZHL.cif", auth_chains=True)
            )
        self.assertEqual([r.id for r in records], ["4ZHL:P", "4ZHL:U"])
        self.assertEqual(self.chains(records), ["P", "U"])
        self.assertEqual(records[0].seq, "CPAYSRYIGC")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            atom_records = list(
                SeqIO.parse(support.DATA / "PDB" / "4ZHL.cif", "cif-atom")
            )
        self.assertEqual(self.chains(atom_records), ["P", "U"])

    def test_label_chains(self):
        """Name chains by label id without a warning with auth_chains=False."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonDeprecationWarning)
            records = list(
                CifSeqresIterator(support.DATA / "PDB" / "4ZHL.cif", auth_chains=False)
            )
        self.assertEqual([r.id for r in records], ["4ZHL:A", "4ZHL:B"])
        self.assertEqual(self.chains(records), ["A", "B"])
        self.assertEqual(records[1].seq, "CPAYSRYIGC")

    def test_dbxrefs_follow_chain(self):
        """Attach each cross-reference to the chain it belongs to."""
        # 1LCD: the LacI protein is author chain A but label chain C, and its
        # _struct_ref_seq entry names it by the author id.
        laci = ["UNP:P03023", "UNP:LACI_ECOLI"]
        records = list(
            CifSeqresIterator(support.DATA / "PDB" / "1LCD.cif", auth_chains=False)
        )
        self.assertEqual([r.dbxrefs == laci for r in records], [False, False, True])
        self.assertTrue(records[2].seq.startswith("MKPVTLYDVAEY"))
        records = list(
            CifSeqresIterator(support.DATA / "PDB" / "1LCD.cif", auth_chains=True)
        )
        self.assertEqual([r.dbxrefs == laci for r in records], [True, False, False])
        # With author ids, the records match those read from the PDB file
        pdb_records = list(SeqIO.parse(support.DATA / "PDB" / "1LCD.pdb", "pdb-seqres"))
        for record, pdb_record in zip(records, pdb_records, strict=True):
            self.assertEqual(record.id, pdb_record.id)
            self.assertEqual(record.annotations, pdb_record.annotations)
            self.assertEqual(record.seq, pdb_record.seq)
            self.assertEqual(record.dbxrefs, pdb_record.dbxrefs)

    def test_no_author_ids(self):
        """Fall back to the label ids where the file has no author ids."""
        no_column = (
            "data_TEST\n"
            "loop_\n"
            "_pdbx_poly_seq_scheme.asym_id\n"
            "_pdbx_poly_seq_scheme.mon_id\n"
            "A MET\n"
            "A ALA\n"
            "B GLY\n"
        )
        null_values = (
            "data_TEST\n"
            "loop_\n"
            "_pdbx_poly_seq_scheme.asym_id\n"
            "_pdbx_poly_seq_scheme.mon_id\n"
            "_pdbx_poly_seq_scheme.pdb_strand_id\n"
            "A MET ?\n"
            "A ALA ?\n"
            "B GLY .\n"
        )
        for data in (no_column, null_values):
            for auth_chains in (None, True, False):
                with warnings.catch_warnings():
                    warnings.simplefilter("error", BiopythonDeprecationWarning)
                    records = list(CifSeqresIterator(StringIO(data), auth_chains))
                self.assertEqual(self.chains(records), ["A", "B"])
                self.assertEqual([r.seq for r in records], ["MA", "G"])


def AtomTestGenerator(extension, parser):
    """Test factory for tests reading ATOM (or similar) records.

    See SeqresTestGenerator for more information.
    """

    class AtomTests(unittest.TestCase):
        def test_atom_parse(self):
            """Parse a multi-chain structure by ATOM entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=2BEG
            """
            chains = list(
                SeqIO.parse(support.DATA / "PDB" / f"2BEG.{extension}", parser)
            )
            self.assertEqual(len(chains), 5)
            actual_seq = "LVFFAEDVGSNKGAIIGLMVGGVVIA"
            for chain, chn_id in zip(chains, "ABCDE"):
                self.assertEqual(chain.id, "2BEG:" + chn_id)
                self.assertEqual(chain.annotations["chain"], chn_id)
                self.assertEqual(chain.annotations["model"], 0)
                self.assertEqual(chain.seq, actual_seq)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                chains = list(
                    SeqIO.parse(support.DATA / "PDB" / f"2XHE.{extension}", parser)
                )
            actual_seq = (
                "DRLSRLRQMAAENQXXXXXXXXXXXXXXXXXXXXXXXPEPFMADFFNRVK"
                "RIRDNIEDIEQAIEQVAQLHTESLVAVSKEDRDRLNEKLQDTMARISALG"
                "NKIRADLKQIEKENKRAQQEGTFEDGTVSTDLRIRQSQHSSLSRKFVKVM"
                "TRYNDVQAENKRRYGENVARQCRVVEPSLSDDAIQKVIEHGXXXXXXXXX"
                "XXXXXXXXNEIRDRHKDIQQLERSLLELHEMFTDMSTLVASQGEMIDRIE"
                "FSVEQSHNYV"
            )
            self.assertEqual(chains[1].seq, actual_seq)

        def test_atom_read(self):
            """Read a single-chain structure by ATOM entries.

            Reference:
            http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=1A8O
            """
            chain = SeqIO.read(support.DATA / "PDB" / f"1A8O.{extension}", parser)
            self.assertEqual(chain.id, "1A8O:A")
            self.assertEqual(chain.annotations["chain"], "A")
            self.assertEqual(chain.annotations["model"], 0)
            self.assertEqual(
                chain.seq,
                "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTIL"
                "KALGPGATLEEMMTACQG",
            )

    return AtomTests


class TestPdbAtom(AtomTestGenerator("pdb", "pdb-atom")):
    """Test pdb-atom SeqIO driver."""

    def test_atom_noheader(self):
        """Parse a PDB with no HEADER line."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chains = list(SeqIO.parse(support.DATA / "PDB" / "1LCD.pdb", "pdb-atom"))

        self.assertEqual(len(chains), 1)
        self.assertEqual(
            chains[0].seq, "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR"
        )

    def test_atom_read_noheader(self):
        """Read a single-chain PDB without a header by ATOM entries."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chain = SeqIO.read(support.DATA / "PDB" / "a_structure.pdb", "pdb-atom")
        self.assertEqual(chain.id, "????:A")
        self.assertEqual(chain.annotations["chain"], "A")
        self.assertEqual(chain.seq, "Q")

    def test_atom_with_insertion(self):
        """Read a PDB with residue insertion code."""
        chain = SeqIO.read(support.DATA / "PDB" / "2n0n_M1.pdb", "pdb-atom")
        self.assertEqual(chain.seq, "HAEGKFTSEF")


class TestCifAtom(AtomTestGenerator("cif", "cif-atom")):
    """Test cif-atom SeqIO driver."""

    def test_atom_read_noheader(self):
        """Read a single-chain CIF without a header by ATOM entries."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            warnings.simplefilter("ignore", BiopythonParserWarning)
            chain = SeqIO.read(support.DATA / "PDB" / "a_structure.cif", "cif-atom")
        self.assertEqual(chain.id, "????:A")
        self.assertEqual(chain.annotations["chain"], "A")
        self.assertEqual(
            chain.seq,
            "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
