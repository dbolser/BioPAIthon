"""
Tests for BinaryCIF code in the PDB package.
"""

import sys
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.PDB."
    ) from None

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

    def test_output_too_large(self):
        for input_dtype, output_dtype in [
            (np.uint8, np.uint32),
            (np.uint16, np.uint32),
            (np.int8, np.int32),
            (np.int16, np.int32),
        ]:
            with self.subTest(input_dtype=input_dtype):
                packed = np.array([1, 2], dtype=input_dtype)
                unpacked = np.empty(3, dtype=output_dtype)
                with self.assertRaisesRegex(ValueError, "too large"):
                    _bcif_helper.integer_unpack(packed, unpacked)

    def test_wrong_src_size(self):
        for src_size, message in [(1, "too small"), (3, "too large")]:
            column = {
                "data": {
                    "data": np.array([1, 2], dtype=np.uint8),
                    "encoding": [
                        {
                            "kind": "IntegerPacking",
                            "byteCount": 1,
                            "srcSize": src_size,
                            "isUnsigned": True,
                        }
                    ],
                }
            }
            with self.subTest(src_size=src_size):
                with self.assertRaisesRegex(ValueError, message):
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

    def test_byte_swapped_buffers_rejected(self):
        """Test that a buffer in the wrong byte order is refused.

        integer_unpack reads and writes through native-width pointers and
        pays no attention to the declared byte order, so a byte swapped
        buffer would decode to the wrong values. Refusing it means a
        big-endian machine gets an error rather than wrong coordinates.
        """
        swapped = ">" if sys.byteorder == "little" else "<"
        with self.subTest("input"):
            packed = np.array([1, 2], dtype=f"{swapped}u2")
            unpacked = np.empty(2, dtype=np.uint32)
            with self.assertRaisesRegex(ValueError, "native byte order"):
                _bcif_helper.integer_unpack(packed, unpacked)
        with self.subTest("output"):
            packed = np.array([1, 2], dtype=np.uint8)
            unpacked = np.empty(2, dtype=f"{swapped}u4")
            with self.assertRaisesRegex(ValueError, "native byte order"):
                _bcif_helper.integer_unpack(packed, unpacked)

    def test_explicit_native_byte_order_accepted(self):
        """Test that an explicitly native-tagged buffer is not refused.

        NumPy reports "=u2" as a format string carrying a byte order
        character, so parsing only the first character would reject it.
        """
        packed = np.array([1, 2], dtype="=u2")
        unpacked = np.empty(2, dtype=np.uint32)
        _bcif_helper.integer_unpack(packed, unpacked)
        self.assertEqual(unpacked.tolist(), [1, 2])

    def test_decoder_converts_byte_order(self):
        """Test that the decoder hands the helper a native buffer.

        The packed data is little-endian on disk, so on a big-endian machine
        it reaches the decoder byte swapped. On a little-endian machine this
        exercises the no-op path and only checks the conversion is harmless.
        """
        swapped = ">" if sys.byteorder == "little" else "<"
        column = {
            "data": {
                "data": np.array([1, 2], dtype=f"{swapped}u1"),
                "encoding": [
                    {
                        "kind": "IntegerPacking",
                        "byteCount": 1,
                        "srcSize": 2,
                        "isUnsigned": True,
                    }
                ],
            }
        }
        _integer_packing_decoder(column)
        self.assertEqual(column["data"]["data"].tolist(), [1, 2])
        self.assertEqual(
            column["data"]["data"].dtype.byteorder in ("=", "|"),
            True,
            "decoded output should be in native byte order",
        )


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
