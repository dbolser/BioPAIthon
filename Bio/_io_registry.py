# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Lazy tables mapping file format names to their handlers (PRIVATE).

Bio.SeqIO and Bio.Align look format names up in these tables, so that a
format module is imported only when someone uses that format.  Their tables
also hold the formats that installed distributions declare as entry points.
"""

import importlib
import threading
import warnings
from typing import Any

from Bio import BiopythonWarning

_ABSENT = object()

# The installed distributions are scanned for entry points once per process,
# when a table first needs them.  _installed holds what
# importlib.metadata.entry_points() returned, and _found the (name, spec,
# distribution) of each entry point in a group.
_scan_lock = threading.Lock()
_scan_failed = False
_installed: Any = None
_found: dict[str, list[tuple[str, str, str]]] = {}


def _entry_points(group):
    """Return (name, spec, distribution) for each entry point in a group (PRIVATE).

    Each spec is a "module" or "module:attr" string, as _resolve takes.  It
    is built from the entry point's module and attr, so an "[extra]" suffix
    on its value does not matter.  If the scan raises, warn once and report
    no entry points.
    """
    global _installed, _scan_failed
    failure = None
    with _scan_lock:
        found = _found.get(group)
        if found is None:
            found = []
            if not _scan_failed:
                try:
                    if _installed is None:
                        import importlib.metadata

                        _installed = importlib.metadata.entry_points()
                    for entry_point in _installed.select(group=group):
                        spec = entry_point.module
                        if entry_point.attr:
                            spec += ":" + entry_point.attr
                        distribution = entry_point.dist.name
                        found.append((entry_point.name, spec, distribution))
                except Exception as exception:
                    _scan_failed = True
                    found = []
                    failure = exception
            _found[group] = found
    if failure is not None:
        warnings.warn(
            f"Could not look for file format plugins, so none are used: {failure!r}",
            BiopythonWarning,
        )
    return found


def _resolve(spec):
    """Import and return the object a spec string names (PRIVATE).

    A spec follows ``importlib.metadata.EntryPoint.load``: "pkg.mod" imports
    and returns the module, and "pkg.mod:attr.path" imports pkg.mod and then
    follows the dotted attribute path.
    """
    module_name, _, attributes = spec.partition(":")
    value = importlib.import_module(module_name)
    for attribute in attributes.split("."):
        if attribute:
            value = getattr(value, attribute)
    return value


def _same_handler(stored, value):
    """Return whether a stored entry and a value are the same handler (PRIVATE).

    A spec string and an object are the same handler when the spec resolves
    to that object, so the answer does not depend on whether the entry has
    been used yet.  Comparing the two imports the spec.  Two spec strings
    are compared as strings.
    """
    if stored is value:
        return True
    if isinstance(stored, str) and isinstance(value, str):
        return stored == value
    if isinstance(stored, str):
        return _resolve(stored) is value
    if isinstance(value, str):
        return _resolve(value) is stored
    return False


def _qualify(value, package):
    """Turn a short "Module.attr" value into "package.Module:attr" (PRIVATE).

    Spec strings that already contain a colon, and values that are not
    strings, are returned unchanged.
    """
    if isinstance(value, str) and ":" not in value:
        module_name, _, attributes = value.partition(".")
        return f"{package}.{module_name}:{attributes}"
    return value


class FormatRegistry(dict):
    """Map format names to handlers, importing each handler on first use (PRIVATE).

    Each value is one of:

     - a spec string, "pkg.mod" or "pkg.mod:attr.path", imported on first
       access as importlib.metadata.EntryPoint.load would;
     - None, resolved on first access by calling factory(name);
     - anything else, returned as it is.

    The resolved value then replaces the spec or None, so the import cost of
    each format is paid once, and only by callers who use that format.

    With package given, a colon-free "Module.attr" value is read as
    "package.Module:attr" when the table is built.  This is the short form in
    which Bio.SeqIO writes its tables.

    Membership tests, len() and iterating over the names work as for any
    dict, and import no format.  get(name, default) returns default only when
    the name is absent; an error raised while resolving a name that is present
    propagates.  values() and items() resolve every entry, and return lists.

    The names present when the table is built are its built-in names, kept in
    the builtin attribute.  register() adds a name, or replaces one on request;
    the packages check that the name and value suit them before calling it.

    With group given, the table also holds the formats which installed
    distributions declare in that entry-point group.  It looks for them once:
    on its first miss (a lookup or membership test of an absent name), or the
    first time its names are listed (iteration, keys(), len(), values() or
    items()), but never on a hit, so using a built-in format does not import
    importlib.metadata.  name_rule(name) returns the format name for an entry
    point's name, or raises ValueError if the package cannot use it.  Such an
    entry point is skipped with a warning, as is one naming a built-in format,
    and a name which two distributions give different objects.  A name that
    register() stored first is skipped silently.  The names added from entry
    points are kept in the plugins attribute; register() replaces them
    without replace=True.  An entry point's object is imported when its format
    is first used, and an error importing it then propagates.
    """

    def __init__(
        self, specs, factory=None, *, package=None, group=None, name_rule=None
    ):
        """Initialize from a mapping of format name to value."""
        if package is not None:
            specs = {name: _qualify(value, package) for name, value in specs.items()}
        super().__init__(specs)
        self.builtin = frozenset(super().keys())
        self.plugins = set()
        self._factory = factory
        self._group = group
        self._name_rule = name_rule
        self._discovered = group is None
        self._lock = threading.Lock()

    def _discover(self):
        """Add the formats declared in this table's entry-point group (PRIVATE).

        This happens once.  Return whether it had not yet happened when
        called, as only then can a name missing before be present now.
        """
        if self._discovered:
            return False
        group = self._group
        offers: dict[str, dict[str, set[str]]] = {}  # name -> spec -> distributions
        messages = []
        for name, spec, distribution in _entry_points(group):
            if self._name_rule is not None:
                try:
                    name = self._name_rule(name)
                except ValueError as exception:
                    messages.append(
                        f"Ignoring entry point {name!r} of {distribution} in group"
                        f" {group!r}: {exception}"
                    )
                    continue
            offers.setdefault(name, {}).setdefault(spec, set()).add(distribution)
        with self._lock:
            if self._discovered:  # done meanwhile, by another thread
                return True
            for name, targets in offers.items():
                distributions = " and ".join(sorted(set().union(*targets.values())))
                if name in self.builtin:
                    messages.append(
                        f"Ignoring entry point {name!r} of {distributions} in group"
                        f" {group!r}: {name!r} is a built-in format"
                    )
                elif len(targets) > 1:
                    messages.append(
                        f"Ignoring entry point {name!r} in group {group!r}:"
                        f" {distributions} give it different objects"
                    )
                elif not super().__contains__(name):
                    [spec] = targets
                    super().__setitem__(name, spec)
                    self.plugins.add(name)
            self._discovered = True
        # Warn outside the lock, as a warning may run arbitrary code:
        for message in messages:
            warnings.warn(message, BiopythonWarning)
        return True

    def _handler(self, name, value):
        """Return the handler for name, given its stored value (PRIVATE)."""
        if isinstance(value, str):
            resolved = _resolve(value)
        elif value is None:
            resolved = self._factory(name)
        else:
            return value
        # Resolving happens outside the lock, since imports take their own
        # locks.  Storing happens under it, and only if the entry is still the
        # value read above, so that concurrent first accesses agree on one
        # resolved value (in particular, on one wrapper class), and a value
        # assigned in the meantime is not overwritten.  An entry deleted in
        # the meantime (SeqIO.register_format drops convert shortcuts) was
        # present when looked up, so the lookup still returns it.
        with self._lock:
            current = super().get(name, _ABSENT)
            if current is value:
                super().__setitem__(name, resolved)
                return resolved
            if current is _ABSENT:
                return resolved
        return self[name]

    def _stored(self, name):
        """Return the value stored for name, or _ABSENT (PRIVATE).

        On a miss, add the entry-point formats if not yet done and look again.
        """
        value = super().get(name, _ABSENT)
        if value is _ABSENT and self._discover():
            value = super().get(name, _ABSENT)
        return value

    def __getitem__(self, name):
        """Return the handler for this format, importing it if needed."""
        value = self._stored(name)
        if value is _ABSENT:
            raise KeyError(name)
        return self._handler(name, value)

    def get(self, name, default=None):
        """Return the handler for this format, or default if it is absent."""
        value = self._stored(name)
        if value is _ABSENT:
            return default
        return self._handler(name, value)

    def __contains__(self, name):
        """Return whether the format is present."""
        return self._stored(name) is not _ABSENT

    def __iter__(self):
        """Iterate over the format names."""
        self._discover()
        return super().__iter__()

    def keys(self):
        """Return a view of the format names."""
        self._discover()
        return super().keys()

    def __len__(self):
        """Return the number of formats."""
        self._discover()
        return super().__len__()

    def values(self):
        """Return a list of all handlers, importing any not yet imported."""
        return [self[name] for name in self]

    def items(self):
        """Return a list of all (format, handler) pairs, importing as needed."""
        return [(name, self[name]) for name in self]

    def register(self, name, value, *, replace=False):
        """Store value under name, or check it is already there.

        An absent name is added, and a name added from an entry point is
        replaced.  Any other name is replaced only if replace is true.
        Otherwise registering the handler it already holds does nothing, and
        any other handler raises ValueError.  This never looks for entry
        points.
        """
        with self._lock:
            stored = super().get(name, _ABSENT)
            if replace or stored is _ABSENT or name in self.plugins:
                super().__setitem__(name, value)
                self.plugins.discard(name)
                return
        # Outside the lock, as comparing may import:
        if not _same_handler(stored, value):
            raise ValueError(
                f"Format {name!r} already exists; use replace=True to replace it"
            )
