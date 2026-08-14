# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This module attempts to detect whether the internet is available.
# To use it, import requires_internet into your Python code, and call
# requires_internet.check().  If the internet is available, then the
# call returns.  If it is not, then it raises a
# MissingExternalDependencyError exception, and run_tests.py reports
# the calling test module as skipped.

"""Common code to check if the internet is available."""

from Bio import MissingExternalDependencyError

RELIABLE_DOMAIN = "biopython.org"
TIMEOUT = 10.0  # seconds


def _internet_available(timeout=TIMEOUT):
    """Probe DNS resolution in a daemon thread, giving up after timeout (PRIVATE).

    socket.getaddrinfo() itself takes no timeout argument, and on a
    firewalled network that silently drops packets it can block for
    minutes.  Running it in a daemon thread bounds the wait: if the
    probe has not answered within the timeout, report the internet as
    unavailable so the caller skips instead of stalling the test run.
    """
    import socket
    import threading

    results = []

    def probe():
        try:
            socket.getaddrinfo(
                RELIABLE_DOMAIN, 80, socket.AF_UNSPEC, socket.SOCK_STREAM
            )
        except OSError:
            results.append(False)
        else:
            results.append(True)

    thread = threading.Thread(target=probe, daemon=True)
    thread.start()
    thread.join(timeout)
    return bool(results) and results[0]


def check():
    try:
        check.available
    except AttributeError:
        # First call: probe for internet availability (time-boxed).
        check.available = _internet_available()
    if not check.available:
        raise MissingExternalDependencyError("internet not available")
