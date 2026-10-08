# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.File's indexes.

Bio.SeqIO.index and Bio.SearchIO.index return these dictionaries, generic in
the record type. Type-checked with the mypy.ini beside it, never run.
"""

import os
from collections.abc import Iterator
from collections.abc import Mapping
from typing import Any
from typing import overload

from typing_extensions import assert_type

from Bio.File import _IndexedSeqFileDict
from Bio.File import _IndexedSeqFileProxy
from Bio.File import _SQLiteManySeqFilesDict
from Bio.SearchIO._model import QueryResult
from Bio.SeqRecord import SeqRecord


def in_memory(d: _IndexedSeqFileDict[SeqRecord]) -> None:
    """Use an index as Bio.SeqIO.index returns it."""
    assert_type(d["x"], SeqRecord)
    assert_type(d.get("x"), SeqRecord | None)
    assert_type(d.get_raw("x"), bytes)
    assert_type(len(d), int)
    d.close()

    records: Mapping[str, SeqRecord] = d

    # A key is a str unless a key_function says otherwise; see below.
    assert_type(list(d), list[str | Any])
    for key in d:
        key.upper()
        key.frobnicate()  # type: ignore[union-attr]


def search_results(d: _IndexedSeqFileDict[QueryResult]) -> None:
    """Use an index as Bio.SearchIO.index returns it."""
    assert_type(d["x"], QueryResult)


def in_sqlite(d: _SQLiteManySeqFilesDict[SeqRecord]) -> None:
    """Use an index as Bio.SeqIO.index_db returns it."""
    assert_type(d["x"], SeqRecord)
    assert_type(d.get_raw("x"), bytes)
    d.close()


# A record type needs an id, as the index checks each record against its key.
_IndexedSeqFileDict[int]  # type: ignore[type-var]


class Proxy(_IndexedSeqFileProxy[SeqRecord]):
    """A format's random access proxy, as Bio.SeqIO._index defines them."""

    def __iter__(self) -> Iterator[tuple[str, int, int]]:
        """Return (identifier, offset, length) tuples."""
        return iter([])

    def get(self, offset: int) -> SeqRecord:
        """Return the record at this offset."""
        raise NotImplementedError


# The dictionary takes its record type from the proxy.
assert_type(
    _IndexedSeqFileDict(Proxy(), None, "repr", "SeqRecord"),
    _IndexedSeqFileDict[SeqRecord],
)


def make_tuple(identifier: str) -> tuple[int, int]:
    """Make a key from an id, as in the Bio.SeqIO.index docstring."""
    parts = identifier.split("_")
    return int(parts[-2]), int(parts[-1])


# A key_function may return any hashable key, not just a str.
by_tuple = _IndexedSeqFileDict(Proxy(), make_tuple, "repr", "SeqRecord")
assert_type(by_tuple[(540, 792)], SeqRecord)
_IndexedSeqFileDict(Proxy(), str.split, "repr", "SeqRecord")  # type: ignore[arg-type]


# Bio.SeqIO.index_db and Bio.SearchIO.index_db pass a factory that is called
# two ways: with a format alone, and with a format and a filename. Written
# with an overload for each, it keeps the record type.
@overload
def proxy_factory(format: str) -> bool: ...


@overload
def proxy_factory(format: str, filename: str | os.PathLike[str]) -> Proxy: ...


def proxy_factory(
    format: str, filename: str | os.PathLike[str] | None = None
) -> Proxy | bool:
    """Say if a format is supported, or make a proxy for one file."""
    if filename:
        return Proxy()
    return format == "fasta"


assert_type(
    _SQLiteManySeqFilesDict(
        ":memory:", None, proxy_factory, "fasta", None, "repr", "SeqRecord"
    ),
    _SQLiteManySeqFilesDict[SeqRecord],
)
