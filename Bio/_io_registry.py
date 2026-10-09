# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Lazy tables mapping file format names to their handlers (PRIVATE).

Bio.SeqIO and Bio.Align look format names up in these tables, so that a
format module is imported only when someone uses that format.
"""

import importlib
import threading

_ABSENT = object()


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
    dict, and import nothing.  get(name, default) returns default only when
    the name is absent; an error raised while resolving a name that is present
    propagates.  values() and items() resolve every entry, and return lists.

    The names present when the table is built are its built-in names, kept in
    the builtin attribute.  register() adds a name, or replaces one on request;
    the packages check that the name and value suit them before calling it.
    """

    def __init__(self, specs, factory=None, *, package=None):
        """Initialize from a mapping of format name to value."""
        if package is not None:
            specs = {name: _qualify(value, package) for name, value in specs.items()}
        super().__init__(specs)
        self.builtin = frozenset(self)
        self._factory = factory
        self._lock = threading.Lock()

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
        # assigned in the meantime is not overwritten.
        with self._lock:
            if super().__getitem__(name) is value:
                super().__setitem__(name, resolved)
                return resolved
        return self[name]

    def __getitem__(self, name):
        """Return the handler for this format, importing it if needed."""
        return self._handler(name, super().__getitem__(name))

    def get(self, name, default=None):
        """Return the handler for this format, or default if it is absent."""
        try:
            value = super().__getitem__(name)
        except KeyError:
            return default
        return self._handler(name, value)

    def values(self):
        """Return a list of all handlers, importing any not yet imported."""
        return [self[name] for name in self]

    def items(self):
        """Return a list of all (format, handler) pairs, importing as needed."""
        return [(name, self[name]) for name in self]

    def conflicts(self, name, value):
        """Return whether name is present and holds a handler other than value.

        A spec string and an object are the same handler when the spec
        resolves to that object, so the answer does not depend on whether
        the entry has been used yet.  Comparing the two imports the spec.
        """
        stored = super().get(name, _ABSENT)
        if stored is _ABSENT or stored is value:
            return False
        if isinstance(stored, str) and isinstance(value, str):
            return stored != value
        if isinstance(stored, str):
            return _resolve(stored) is not value
        if isinstance(value, str):
            return _resolve(value) is not stored
        return True

    def register(self, name, value, *, replace=False):
        """Store value under name, or check it is already there.

        An absent name is added.  A present name is replaced only if replace
        is true.  Otherwise registering the handler it already holds does
        nothing, and any other handler raises ValueError.
        """
        with self._lock:
            if replace or not super().__contains__(name):
                super().__setitem__(name, value)
                return
        if self.conflicts(name, value):
            raise ValueError(
                f"Format {name!r} already exists; use replace=True to replace it"
            )
