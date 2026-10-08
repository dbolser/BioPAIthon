# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""I/O function wrappers for ``Bio.Nexus`` trees."""

from itertools import chain

from Bio.Nexus import Nexus
from Bio.Nexus.StandardData import NexusError
from Bio.Phylo import Newick
from Bio.Phylo import NewickIO

# Structure of a Nexus tree-only file
NEX_TEMPLATE = """\
#NEXUS
Begin Taxa;
 Dimensions NTax=%(count)d;
 TaxLabels %(labels)s;
End;
Begin Trees;
 %(trees)s
End;
"""

# 'index' starts from 1; 'tree' is the Newick tree string
TREE_TEMPLATE = "Tree tree%(index)d=%(tree)s"

# NewickIO drops the quotes around a label, then reads an internal label such
# as '95' as a number. Bio.Nexus.Trees kept a quoted label as a name, so each
# one is tagged with a NUL character, which has no place in NEXUS text, and
# the tag is removed after parsing.
_QUOTED = "\0"


def _tag_quoted_label(match):
    """Tag a quoted label found by NewickIO's tokenizer (PRIVATE)."""
    token = match.group()
    return f"'{_QUOTED}{token[1:]}" if token.startswith("'") else token


class _NexusWithNewickTrees(Nexus.Nexus):
    """Nexus reader that parses each tree with NewickIO (PRIVATE).

    ``Bio.Nexus.Nexus`` builds a ``Bio.Nexus.Trees.Tree`` for every tree, and
    that parser recurses once per level of nesting, so deep trees raise
    RecursionError. ``NewickIO.Parser`` does not recurse. Everything else in
    the file is still read by ``Bio.Nexus``.
    """

    def _tree(self, options):
        name, weight, rooted, newick = Nexus._split_tree_command(options)
        newick = NewickIO.tokenizer.sub(_tag_quoted_label, newick)
        parsed = next(NewickIO.Parser.from_string(newick).parse(), None)
        # An empty tree description gives a tree of one empty clade, as before
        root = Newick.Clade() if parsed is None else parsed.root
        _match_bio_nexus_trees(root, self.translate)
        self.trees.append(
            Newick.Tree(root=root, rooted=rooted, name=name, weight=weight)
        )


def _match_bio_nexus_trees(root, translate):
    """Give NewickIO's clades the values Bio.Nexus.Trees gave them (PRIVATE).

    Nexus trees used to be parsed by ``Bio.Nexus.Trees``, and callers see
    these differences from plain Newick parsing:

    - a missing branch length is 0.0, not None;
    - a bare number after a clade, with no branch length, is the branch
      length, not the confidence;
    - confidences are floats;
    - comments keep their square brackets, e.g. "[&rate=1.0]";
    - a quoted label is a name, even one that looks like a number, e.g. '95';
    - terminal names come from the TRANSLATE table, if any, quoted by
      ``Bio.Nexus.Nexus.safename``.

    Walks the tree in preorder with a list rather than by recursion, so deep
    trees work and a failed TRANSLATE lookup names the same taxon as before.
    """
    stack = [root]
    while stack:
        clade = stack.pop()
        stack.extend(reversed(clade.clades))
        if clade.branch_length is None:
            # Not "confidence or 0", which would turn -0.0 into 0.0
            clade.branch_length = (
                0.0 if clade.confidence is None else float(clade.confidence)
            )
            clade.confidence = None
        elif clade.confidence is not None:
            clade.confidence = float(clade.confidence)
        if clade.comment is not None:
            clade.comment = f"[{clade.comment}]"
        if clade.name:
            clade.name = clade.name.replace(_QUOTED, "")
        if translate and not clade.clades:
            try:
                clade.name = Nexus.safename(translate[int(clade.name)])
            except (TypeError, ValueError, KeyError):
                raise NexusError(
                    f"Unable to substitute {clade.name} using 'translate' data."
                ) from None


def parse(handle):
    """Parse the trees in a Nexus file.

    ``Bio.Nexus`` reads the file, including any TRANSLATE table, and
    ``NewickIO`` parses each tree description. The values are those the older
    ``Bio.Nexus.Trees`` parser gave: for example, a missing branch length is
    0.0 rather than None, and comments keep their square brackets.
    """
    yield from _NexusWithNewickTrees(handle).trees


def write(obj, handle, **kwargs):
    """Write a new Nexus file containing the given trees.

    Uses a simple Nexus template and the NewickIO writer to serialize just the
    trees and minimal supporting info needed for a valid Nexus file.
    """
    trees = list(obj)
    writer = NewickIO.Writer(trees)
    nexus_trees = [
        TREE_TEMPLATE % {"index": idx + 1, "tree": nwk}
        for idx, nwk in enumerate(
            writer.to_strings(plain=False, plain_newick=True, **kwargs)
        )
    ]
    tax_labels = [str(x.name) for x in chain(*(t.get_terminals() for t in trees))]
    text = NEX_TEMPLATE % {
        "count": len(tax_labels),
        "labels": " ".join(tax_labels),
        "trees": "\n".join(nexus_trees),
    }
    handle.write(text)
    return len(nexus_trees)
