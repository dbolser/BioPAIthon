# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Type stubs for the C extension Bio.Nexus.cnexus."""

# Return the text without its comments and with each command-ending ";"
# replaced by chr(7), or just "[" or "]" if that bracket is unmatched.
def scanfile(text: str, /) -> str: ...
