"""
Tests for BinaryCIF code in the PDB package.
"""

import unittest

import numpy as np

from Bio.PDB import MMCIFParser
from Bio.PDB import _bcif_helper
from Bio.PDB.binary_cif import BinaryCIFParser
from Bio.PDB.binary_cif import _integer_packing_decoder


class TestIntegerUnpack(unittest.TestCase):
    def test_output_too_small(self):
        for input_dtype, output_dtype in [
            (np.uint8, np.uint32),
            (np.uint16, np.uint32),
            (np.int8, np.int32),
            (np.int16, np.int32),
        ]:
            with self.subTest(input_dtype=input_dtype):
                packed = np.array([1, 2], dtype=input_dtype)
                unpacked = np.empty(1, dtype=output_dtype)
                with self.assertRaisesRegex(ValueError, "too small"):
                    _bcif_helper.integer_unpack(packed, unpacked)

    def test_wrong_src_size(self):
        column = {
            "data": {
                "data": np.array([1, 2], dtype=np.uint8),
                "encoding": [
                    {
                        "kind": "IntegerPacking",
                        "byteCount": 1,
                        "srcSize": 1,
                        "isUnsigned": True,
                    }
                ],
            }
        }
        with self.assertRaisesRegex(ValueError, "too small"):
            _integer_packing_decoder(column)

    def test_truncated_packed_integer(self):
        for input_dtype, output_dtype, sentinel in [
            (np.uint8, np.uint32, np.iinfo(np.uint8).max),
            (np.uint16, np.uint32, np.iinfo(np.uint16).max),
            (np.int8, np.int32, np.iinfo(np.int8).min),
            (np.int8, np.int32, np.iinfo(np.int8).max),
            (np.int16, np.int32, np.iinfo(np.int16).min),
            (np.int16, np.int32, np.iinfo(np.int16).max),
        ]:
            with self.subTest(input_dtype=input_dtype, sentinel=sentinel):
                packed = np.array([sentinel], dtype=input_dtype)
                unpacked = np.empty(1, dtype=output_dtype)
                with self.assertRaisesRegex(ValueError, "truncated"):
                    _bcif_helper.integer_unpack(packed, unpacked)

    def test_output_format(self):
        packed = np.array([1], dtype=np.uint8)
        for output_dtype in [np.uint16, np.int32]:
            with self.subTest(output_dtype=output_dtype):
                unpacked = np.empty(1, dtype=output_dtype)
                with self.assertRaisesRegex(ValueError, "32-bit unsigned"):
                    _bcif_helper.integer_unpack(packed, unpacked)

    def test_error_return(self):
        packed = np.array([[1]], dtype=np.uint8)
        unpacked = np.empty(1, dtype=np.uint32)
        with self.assertRaisesRegex(ValueError, "one-dimensional"):
            _bcif_helper.integer_unpack(packed, unpacked)


class TestBinaryCIFParser(unittest.TestCase):
    def test_get_structure(self):
        mmcif_parser = MMCIFParser(auth_chains=False)
        bcif_parser = BinaryCIFParser()

        for entry in ["1GBT", "6WG6", "3JQH"]:
            mmcif_structure = mmcif_parser.get_structure(entry, f"PDB/{entry}.cif")
            bcif_structure = bcif_parser.get_structure(
                entry, f"PDB/{entry.lower()}.bcif.gz"
            )
            self.assertTrue(
                mmcif_structure.strictly_equals(
                    bcif_structure, compare_coordinates=True
                )
            )
