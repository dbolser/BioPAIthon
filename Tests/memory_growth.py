# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Measure how much memory a callable retains per call.

The C extensions are the part of this package where a leak is easy to write
and hard to notice, so they want a regression test that a leak will fail.

LeakSanitizer does not serve here, which is worth recording so it is not
tried again. It reports what is still allocated when the process exits, and
for CPython that is dominated by the interpreter's own import machinery: over
one test module it reported 1.16 MB in 1039 allocations, not one of which had
a Biopython frame anywhere in its stack. Suppressions cannot separate the two
either, because a suppression matches if *any* frame matches, and every call
into our C arrives through ``_PyEval_EvalFrameDefault``. Suppress that and the
leaks we care about go with it. This is why the test_linux CI job runs with
``ASAN_OPTIONS=detect_leaks=0``.

Measuring the growth *across repeated calls* answers the question we actually
have. Interpreter and library start-up costs are paid before the baseline is
taken, so they cancel out, and what remains is attributable to the callable.
It also needs no sanitizer build, so it runs everywhere the suite does.

What it does not see is memory obtained from libc directly. ``tracemalloc``
hooks the ``PyMem_*`` allocator domains, raw included, which is why a leak in
``ccealign`` is caught despite it using ``PyMem_RawMalloc`` throughout; a bare
``malloc`` bypasses all three and never appears in the total. Of the shipped
extensions only ``Bio/cpairwise2module.c`` allocates that way, for its score
matrix, trace matrix and column cache. It is not covered here and cannot be
until those become ``PyMem_Malloc`` - which is not worth doing to a module
already deprecated for removal.
"""

import gc
import platform
import sysconfig
import unittest

try:
    import tracemalloc
except ImportError:
    # PyPy has no _tracemalloc.
    tracemalloc = None


# Reference counting is only predictable on a GIL-enabled CPython, and the
# free-threaded build defers reference release in ways that blur the measure.
can_measure_growth = (
    platform.python_implementation() == "CPython"
    and not sysconfig.get_config_var("Py_GIL_DISABLED")
    and tracemalloc is not None
)

requires_growth_measurement = unittest.skipUnless(
    can_measure_growth, "GIL-enabled CPython memory tracing is required"
)


def bytes_retained_per_call(function, warmup=5, calls=80):
    """Return the mean bytes still traced per call of function.

    The callable is run warmup times first so that any one-off allocation it
    triggers - a lookup table, an imported module, a cached conversion - is
    already paid for when the baseline is taken.
    """
    if tracemalloc is None:
        raise unittest.SkipTest("tracemalloc is not available")
    was_tracing = tracemalloc.is_tracing()
    if not was_tracing:
        tracemalloc.start()
    try:
        for _ in range(warmup):
            function()
        gc.collect()
        baseline = tracemalloc.get_traced_memory()[0]

        for _ in range(calls):
            function()
        gc.collect()
        retained = tracemalloc.get_traced_memory()[0] - baseline
    finally:
        if not was_tracing:
            tracemalloc.stop()
    return retained / calls


def assert_bounded_growth(test, function, limit=16, warmup=5, calls=80):
    """Fail if function retains more than limit bytes per call.

    The default is calibrated rather than guessed, from both directions.
    Across 40 runs of the five callables in test_C_extension_memory, a
    function that does not leak reported at most 0.40 bytes per call, and
    usually 0.00. Deleting the free in Parser_dealloc, so that the printed
    alignment parser leaks one row array per parse, reports 64.4.

    16 sits between the two with a factor of 40 above the noise and 4 below
    that leak. It was not the first choice: 64 also caught that mutant, but
    only at 64.4 against a limit of 64, which is not a margin worth shipping.
    A limit loose enough to look obviously safe - a kilobyte, say - would be a
    test that cannot fail, which is worse than having no test at all.
    """
    per_call = bytes_retained_per_call(function, warmup=warmup, calls=calls)
    test.assertLess(
        per_call,
        limit,
        f"{per_call:.0f} bytes retained per call, over the {limit} byte limit",
    )
