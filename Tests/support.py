"""Shared helpers for the test suite.

The test data lives in subdirectories of ``Tests/``, and historically the
tests addressed it with paths relative to the current working directory,
which forced the whole suite to be run from inside ``Tests/``.  New and
migrated tests should instead build data paths from :data:`DATA`, so that
they pass from any working directory::

    import support

    with open(support.DATA / "GenBank" / "cor6_6.gb") as handle:
        ...

Tests which deliberately exercise relative-path behaviour (for example
``Bio.SeqIO.index_db`` storing filenames relative to the index) may still
change directory, but should change to an absolute location derived from
:data:`DATA` in ``setUp`` and restore the original directory afterwards,
rather than assume the process already started in ``Tests/``.
"""

import pathlib

DATA = pathlib.Path(__file__).parent
