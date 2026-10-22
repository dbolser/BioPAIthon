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

import importlib
import pathlib
import sys
import warnings

DATA = pathlib.Path(__file__).parent

_MISSING = object()


def import_fresh_with_warnings(module_name):
    """Re-execute a module and return it with the warnings it raised.

    Module-level warnings (deprecations, say) fire only when the module is
    first executed, so an earlier import in the same process would hide
    them. This drops the module from ``sys.modules``, imports it afresh
    with every warning recorded, and returns ``(module, caught)`` where
    ``caught`` is the list of :class:`warnings.WarningMessage` objects.

    The interpreter is left as it was found: the original module object
    goes back into ``sys.modules`` and onto its parent package, so other
    tests holding references to it keep seeing one module, not two.
    """
    parent_name, _, child_name = module_name.rpartition(".")
    parent = sys.modules.get(parent_name)
    old_module = sys.modules.pop(module_name, None)
    old_attribute = getattr(parent, child_name, _MISSING)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            module = importlib.import_module(module_name)
    finally:
        if old_module is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = old_module
        if parent is not None:
            if old_attribute is _MISSING:
                if hasattr(parent, child_name):
                    delattr(parent, child_name)
            else:
                setattr(parent, child_name, old_attribute)
    return module, caught
