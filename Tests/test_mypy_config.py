# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Lint the two per-module typing ratchets in the repository's .mypy.ini.

The file ends with two blocks of one-module sections:

- the check_untyped_defs baseline, which modules may only leave, and
- the disallow_untyped_defs allowlist, which modules may only join.

A slip in either can switch mypy off without anything going red. On a
duplicate section mypy prints "section ... already exists", ignores the
whole file and exits 0, so every global option is lost. These tests parse
the file strictly and check the shape of both blocks. That entries are only
ever removed from one block and added to the other is left to review, as a
test cannot see history.
"""

import configparser
import os
import re
import unittest

MYPY_INI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".mypy.ini")

# The first line of each block's banner comment.
BASELINE_BANNER = "# check_untyped_defs ratchet baseline."
ALLOWLIST_BANNER = "# disallow_untyped_defs allowlist"


def _read():
    """Return the text of .mypy.ini."""
    with open(MYPY_INI, encoding="utf-8") as handle:
        return handle.read()


def _parse(text, strict=True):
    """Parse the text as mypy does; if strict, a duplicate raises an error."""
    # mypy itself uses a RawConfigParser, with the default strict=True.
    parser = configparser.RawConfigParser(strict=strict)
    parser.read_string(text, source=MYPY_INI)
    return parser


def _blocks(text):
    """Return the section names of the baseline and allowlist blocks."""
    baseline_start = text.index(BASELINE_BANNER)
    allowlist_start = text.index(ALLOWLIST_BANNER)
    baseline = []
    allowlist = []
    for match in re.finditer(r"^\[(.+)\]", text, re.MULTILINE):
        if match.start() > allowlist_start:
            allowlist.append(match.group(1))
        elif match.start() > baseline_start:
            baseline.append(match.group(1))
    return baseline, allowlist


@unittest.skipUnless(os.path.isfile(MYPY_INI), ".mypy.ini is not in the sdist")
class MypyConfigTests(unittest.TestCase):
    """Check .mypy.ini parses strictly and its ratchet blocks keep their shape."""

    def setUp(self):
        self.text = _read()
        self.baseline, self.allowlist = _blocks(self.text)
        # Leniently, so that a duplicate fails test_parses_strictly alone.
        self.parser = _parse(self.text, strict=False)

    def test_parses_strictly(self):
        """mypy must not be able to ignore the file over a duplicate."""
        try:
            _parse(self.text)
        except configparser.Error as err:
            self.fail(f"mypy would ignore .mypy.ini entirely and exit 0: {err}")

    def test_block_order(self):
        """The allowlist comes after the baseline, and both are populated."""
        self.assertLess(
            self.text.index(BASELINE_BANNER), self.text.index(ALLOWLIST_BANNER)
        )
        self.assertTrue(self.baseline)
        self.assertTrue(self.allowlist)

    def test_no_module_in_both(self):
        """A module joining the allowlist must leave the baseline."""
        both = sorted(set(self.baseline) & set(self.allowlist))
        self.assertEqual(both, [], "delete the baseline section of a module that joins")

    def test_baseline_sections(self):
        """Each baseline section only switches check_untyped_defs off."""
        for name in self.baseline:
            with self.subTest(section=name):
                self.assertEqual(
                    dict(self.parser[name]), {"check_untyped_defs": "False"}
                )

    def test_allowlist_sections(self):
        """Each allowlist section only switches disallow_untyped_defs on."""
        for name in self.allowlist:
            with self.subTest(section=name):
                self.assertEqual(
                    dict(self.parser[name]), {"disallow_untyped_defs": "True"}
                )

    def test_blocks_sorted(self):
        """Both blocks are in plain string order, so additions do not collide."""
        self.assertEqual(self.baseline, sorted(self.baseline))
        self.assertEqual(self.allowlist, sorted(self.allowlist))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
