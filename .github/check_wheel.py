#!/usr/bin/env python3
"""Check an installed BioPAIthon wheel; cibuildwheel's test-command runs this.

It imports every compiled extension by name, so a wheel that imports ``Bio``
but silently shipped no C extensions fails. On a free-threaded build of
Python it then checks that the GIL is still disabled: importing an extension
that has not declared it can run without the GIL turns the GIL back on for
the whole process, with a ``RuntimeWarning`` naming the extension.

It lives in ``.github/``, not the repository root, because Python puts the
script's own directory first on ``sys.path``, and the source tree's ``Bio``
there would shadow the installed one.
"""

import importlib
import sys
import sysconfig

import Bio

EXTENSIONS = [
    "Bio.Align._aligncore",
    "Bio.Align._alignmentcounts",
    "Bio.Align._codonaligner",
    "Bio.Align._pairwisealigner",
    "Bio.Align.substitution_matrices._arraycore",
    "Bio.Cluster._cluster",
    "Bio.Nexus.cnexus",
    "Bio.PDB._bcif_helper",
    "Bio.PDB.ccealign",
    "Bio.PDB.kdtrees",
    "Bio.SeqIO._twoBitIO",
    "Bio.motifs._pwm",
]

for name in EXTENSIONS:
    importlib.import_module(name)

if sysconfig.get_config_var("Py_GIL_DISABLED") and sys._is_gil_enabled():
    sys.exit("Importing BioPAIthon's C extensions turned the GIL on.")

print(Bio.__version__)
