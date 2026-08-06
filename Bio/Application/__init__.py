# Copyright 2001-2004 Brad Chapman.
# Revisions copyright 2009-2013 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Framework for wrapping command line tools (REMOVED).

Bio.Application has been removed from Biopython. It provided AbstractCommandline and the option classes that the Bio.*.Applications wrappers were built on. It was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Build the command line yourself and run it with the standard library subprocess module, which gives the same control over arguments, stdin, stdout and the return code. See DEPRECATED.rst.
"""

raise ImportError(
    "Bio.Application has been removed from Biopython. It provided AbstractCommandline and the option classes that the Bio.*.Applications wrappers were built on. It was declared obsolete in release 1.79, deprecated in 1.82, and removed in 1.86. Build the command line yourself and run it with the standard library subprocess module, which gives the same control over arguments, stdin, stdout and the return code. See DEPRECATED.rst."
)
