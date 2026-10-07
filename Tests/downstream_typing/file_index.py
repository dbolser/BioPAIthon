# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.File's indexes.

Bio.SeqIO.index and Bio.SearchIO.index return these dictionaries, generic in
the record type. Type-checked with the mypy.ini beside it, never run.
"""

from collections.abc import Iterator
from collections.abc import Mapping

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
    assert_type(list(d), list[str])
    assert_type(len(d), int)
    d.close()

    records: Mapping[str, SeqRecord] = d

    d[1]  # type: ignore[index]


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
