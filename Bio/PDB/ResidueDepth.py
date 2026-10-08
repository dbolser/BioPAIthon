# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# Copyright (C) 2017, Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Calculation of residue depth using command line tool MSMS.

This module uses Michel Sanner's MSMS program for the surface calculation.
See: http://mgltools.scripps.edu/packages/MSMS

Residue depth is the average distance of the atoms of a residue from
the solvent accessible surface.

Residue Depth::

    from Bio.PDB.ResidueDepth import ResidueDepth
    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser()
    structure = parser.get_structure("1a8o", "Tests/PDB/1A8O.pdb")
    model = structure[0]
    rd = ResidueDepth(model)
    print(rd['A',(' ', 152, ' ')])

Direct MSMS interface, typical use::

    from Bio.PDB.ResidueDepth import get_surface
    surface = get_surface(model)

The surface is a NumPy array with all the surface vertices.

Distance to surface::

    from Bio.PDB.ResidueDepth import min_dist
    coord = (1.113, 35.393,  9.268)
    dist = min_dist(coord, surface)

where coord is the coord of an atom within the volume bound by
the surface (ie. atom depth).

To calculate the residue depth (average atom depth of the atoms
in a residue)::

    from Bio.PDB.ResidueDepth import residue_depth
    chain = model['A']
    res152 = chain[152]
    rd = residue_depth(res152, surface)

"""

import functools
import os
import re
import subprocess
import tempfile
import warnings

import numpy as np

from Bio import BiopythonWarning
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import is_aa

# Table 1: Atom Type to radius
_atomic_radii = {
    #   atom num dist  Rexplicit Runited-atom
    1: (0.57, 1.40, 1.40),
    2: (0.66, 1.40, 1.60),
    3: (0.57, 1.40, 1.40),
    4: (0.70, 1.54, 1.70),
    5: (0.70, 1.54, 1.80),
    6: (0.70, 1.54, 2.00),
    7: (0.77, 1.74, 2.00),
    8: (0.77, 1.74, 2.00),
    9: (0.77, 1.74, 2.00),
    10: (0.67, 1.74, 1.74),
    11: (0.70, 1.74, 1.86),
    12: (1.04, 1.80, 1.85),
    13: (1.04, 1.80, 1.80),  # P, S, and LonePairs
    14: (0.70, 1.54, 1.54),  # non-protonated nitrogens
    15: (0.37, 1.20, 1.20),  # H, D  hydrogen and deuterium
    16: (0.70, 0.00, 1.50),  # obsolete entry, purpose unknown
    17: (3.50, 5.00, 5.00),  # pseudoatom - big ball
    18: (1.74, 1.97, 1.97),  # Ca calcium
    19: (1.25, 1.40, 1.40),  # Zn zinc    (traditional radius)
    20: (1.17, 1.40, 1.40),  # Cu copper  (traditional radius)
    21: (1.45, 1.30, 1.30),  # Fe heme iron
    22: (1.41, 1.49, 1.49),  # Cd cadmium
    23: (0.01, 0.01, 0.01),  # pseudoatom - tiny dot
    24: (0.37, 1.20, 0.00),  # hydrogen vanishing if united-atoms
    25: (1.16, 1.24, 1.24),  # Fe not in heme
    26: (1.36, 1.60, 1.60),  # Mg magnesium
    27: (1.17, 1.24, 1.24),  # Mn manganese
    28: (1.16, 1.25, 1.25),  # Co cobalt
    29: (1.17, 2.15, 2.15),  # Se selenium
    30: (3.00, 3.00, 3.00),  # obsolete entry, original purpose unknown
    31: (1.15, 1.15, 1.15),  # Yb ytterbium +3 ion --- wild guess only
    38: (0.95, 1.80, 1.80),  # obsolete entry, original purpose unknown
}

# Table 2: Residue name, atom name and element to Atom Type
# MSMS's pdb_to_xyzr matches residue and atom names against the patterns in
# its atmtypenumbers file, in order, and the first match wins. This table
# follows the same scheme. Each rule is (set of residue names, atom name
# pattern, element, atom type), where None matches anything and the atom name
# pattern is a regular expression that must match the whole name. Order
# matters. Hydrogens and HETATM ions are dealt with before this table, in
# _get_atom_radius.
# fmt: off
_atom_types = tuple(
    (
        residues,
        None if names is None else re.compile(names, re.DOTALL),
        element,
        atom_type,
    )
    for residues, names, element, atom_type in (
        ({"ACE"}, "CA", None, 9),
        # Main chain atoms
        (None, "N", None, 4),
        (None, "CA", None, 7),
        (None, "C", None, 10),
        (None, "O", None, 1),
        (None, "P", None, 13),
        # CB atoms
        ({"ALA"}, "CB", None, 9),
        ({"ILE", "THR", "VAL"}, "CB", None, 7),
        (None, "CB", None, 8),
        # CG atoms
        ({"ASN", "ASP", "ASX", "HIS", "HIP", "HIE", "HID", "HISN", "HISL",
          "LEU", "PHE", "TRP", "TYR"}, "CG", None, 10),
        ({"LEU"}, "CG", None, 7),
        (None, "CG", None, 8),
        # General amino acids in alphabetical order
        ({"GLN"}, None, "O", 3),
        ({"ACE"}, "CH3", None, 9),
        ({"ARG"}, "CD", None, 8),
        ({"ARG"}, "NE|RE", None, 4),
        ({"ARG"}, "CZ", None, 10),
        ({"ARG"}, "(NH|RH).*", None, 5),
        ({"ASN"}, "OD1", None, 1),
        ({"ASN"}, "ND2", None, 5),
        ({"ASN"}, "AD.*", None, 3),
        ({"ASP"}, "(OD|ED).*", None, 3),
        ({"ASX"}, "OD1.*", None, 1),
        ({"ASX"}, "ND2", None, 3),
        ({"ASX"}, "(OD|AD).*", None, 3),
        ({"CYS", "CYX", "CYM"}, "SG", None, 13),
        ({"CYS", "MET"}, "LP.*", None, 13),
        ({"CUH"}, "SG", None, 12),
        ({"GLU"}, "(OE|EE).*", None, 3),
        ({"GLU", "GLN", "GLX"}, "CD", None, 10),
        ({"GLN"}, "OE1", None, 1),
        ({"GLN"}, "NE2", None, 5),
        ({"GLN", "GLX"}, "AE.*", None, 3),
        # Histidines and friends
        # There are 4 kinds of HIS rings: HIS (no protons), HID (proton on Delta),
        #   HIE (proton on epsilon), and HIP (protons on both)
        # Protonated nitrogens are numbered 4, else 14
        # HIS is treated here as the same as HIE
        #
        # HISL is a deprotonated HIS (the L means liganded)
        ({"HIS", "HID", "HIE", "HIP", "HISL"}, "CE1|CD2", None, 11),
        ({"HIS", "HID", "HIE", "HISL"}, "ND1", None, 14),
        ({"HID", "HIP"}, "ND1|RD1", None, 4),
        ({"HIS", "HIE", "HIP"}, "NE2|RE2", None, 4),
        ({"HID", "HISL"}, "NE2|RE2", None, 14),
        ({"HIS", "HID", "HIP", "HISL"}, "(AD|AE).*", None, 4),
        # More amino acids
        ({"ILE"}, "CG1", None, 8),
        ({"ILE"}, "CG2", None, 9),
        ({"ILE"}, "CD|CD1", None, 9),
        ({"LEU"}, "CD.*", None, 9),
        ({"LYS"}, "CG|CD|CE", None, 8),
        ({"LYS"}, "NZ|KZ", None, 6),
        ({"MET"}, "SD", None, 13),
        ({"MET"}, "CE", None, 9),
        ({"PHE"}, "(CD|CE|CZ).*", None, 11),
        ({"PRO"}, "CG|CD", None, 8),
        ({"CSO"}, "SE|SEG", None, 9),
        ({"CSO"}, "OD.*", None, 3),
        ({"SER"}, "OG", None, 2),
        ({"THR"}, "OG1", None, 2),
        ({"THR"}, "CG2", None, 9),
        ({"TRP"}, "CD1", None, 11),
        ({"TRP"}, "CD2|CE2", None, 10),
        ({"TRP"}, "NE1", None, 4),
        ({"TRP"}, "CE3|CZ2|CZ3|CH2", None, 11),
        ({"TYR"}, "CD1|CD2|CE1|CE2", None, 11),
        ({"TYR"}, "CZ", None, 10),
        ({"TYR"}, "OH", None, 2),
        ({"VAL"}, "CG1|CG2", None, 9),
        (None, "CD", None, 8),
        # Co-factors, and other weirdos
        ({"FS3", "FS4"}, "FE.*[1-7]", None, 21),
        ({"FS3", "FS4"}, "S.*[1-7]", None, 13),
        ({"FS3"}, "OXO", None, 1),
        ({"FEO"}, "FE1|FE2", None, 21),
        ({"HEM"}, "O1|O2", None, 1),
        ({"HEM"}, "FE", None, 21),
        ({"HEM"}, "CH[A-D]|CA[BC]|CB[BC]", None, 11),
        ({"HEM"}, "N ?[A-D]", None, 14),
        ({"HEM"}, "C[1-4][A-D]|CG[AD]", None, 10),
        ({"HEM"}, "CM[A-D]", None, 9),
        ({"HEM"}, "OH2", None, 2),
        ({"AZI"}, "N[1-3]", None, 14),
        ({"MPD"}, "C1|C5|C6", None, 9),
        ({"MPD"}, "C2", None, 10),
        ({"MPD"}, "C3", None, 8),
        ({"MPD"}, "C4", None, 7),
        ({"MPD"}, "O7|O8", None, 2),
        ({"SO4", "SUL"}, "S", None, 13),
        ({"SO4", "SUL", "PO4", "PHO"}, "O[1-4]", None, 3),
        ({"PC "}, "O[1-4]", None, 3),
        ({"PC "}, "P1", None, 13),
        ({"PC "}, "C1|C2", None, 8),
        ({"PC "}, "C3|C4|C5", None, 9),
        ({"PC "}, "N1", None, 14),
        ({"BIG"}, "BAL", None, 17),
        ({"POI", "DOT"}, "POI|DOT", None, 23),
        ({"FMN"}, "N1|N5|N10", None, 4),
        ({"FMN"}, "C2|C4|C7|C8|C10|C4A|C5A|C9A", None, 10),
        ({"FMN"}, "O2|O4", None, 1),
        ({"FMN"}, "N3", None, 14),
        ({"FMN"}, "C6|C9", None, 11),
        ({"FMN"}, "C7M|C8M", None, 9),
        ({"FMN"}, "C[1-5].*", None, 8),
        ({"FMN"}, "O[2-4].*", None, 2),
        ({"FMN"}, "O5.*", None, 3),
        ({"FMN"}, "OP[1-3]", None, 3),
        ({"ALK", "MYR"}, "OT1", None, 3),
        ({"ALK", "MYR"}, "C01", None, 10),
        ({"ALK"}, "C16", None, 9),
        ({"MYR"}, "C14", None, 9),
        ({"ALK", "MYR"}, "C.*", None, 8),
        # Metals
        (None, None, "CU", 20),
        (None, None, "ZN", 19),
        (None, None, "MN", 27),
        (None, None, "FE", 25),
        (None, None, "MG", 26),
        (None, None, "CO", 28),
        (None, None, "SE", 29),
        (None, None, "YB", 31),
        # Others
        (None, "SEG", None, 9),
        (None, "OXT", None, 3),
        # Catch-alls
        (None, "(OT|E).*", None, 3),
        (None, "S.*", None, 13),
        (None, "C.*", None, 7),
        (None, "A.*", None, 11),
        (None, "O.*", None, 1),
        (None, "(N|R).*", None, 4),
        (None, "K.*", None, 6),
        (None, "P[A-D]", None, 13),
        (None, "P.*", None, 13),
        ({"FAD", "NAD", "AMX", "APU"}, "O.*", None, 1),
        ({"FAD", "NAD", "AMX", "APU"}, "N.*", None, 4),
        ({"FAD", "NAD", "AMX", "APU"}, "C.*", None, 7),
        ({"FAD", "NAD", "AMX", "APU"}, "P.*", None, 13),
        ({"FAD", "NAD", "AMX", "APU"}, "H.*", None, 15),
    )
)
# fmt: on


# Cached: scanning the table is several times slower per atom than the old
# elif chain, but a structure has only a few hundred distinct (residue, atom,
# element) keys.
@functools.lru_cache(maxsize=4096)
def _get_atom_type(resname, at_name, at_elem):
    """Return the atom type of the first rule in _atom_types to match (PRIVATE).

    Returns None if no rule matches.
    """
    for residues, names, element, atom_type in _atom_types:
        if (
            (residues is None or resname in residues)
            and (names is None or names.fullmatch(at_name))
            and (element is None or at_elem == element)
        ):
            return atom_type
    return None


def _get_atom_radius(atom, rtype="united"):
    """Translate an atom object to an atomic radius defined in MSMS (PRIVATE).

    Uses information from the parent residue and the atom object to define
    the atom type.

    Returns the radius (float) according to the selected type:
     - explicit (reads hydrogens)
     - united (default)

    """
    if rtype == "explicit":
        typekey = 1
    elif rtype == "united":
        typekey = 2
    else:
        raise ValueError(
            f"Radius type ({rtype!r}) not understood. Must be 'explicit' or 'united'"
        )

    resname = atom.parent.resname
    het_atm = atom.parent.id[0]

    at_name = atom.name
    at_elem = atom.element

    # Hydrogens
    if at_elem in ("H", "D"):
        atom_type = 15
    # HETATMs
    elif het_atm == "W" and at_elem == "O":
        atom_type = 2
    elif het_atm != " " and at_elem == "CA":
        atom_type = 18
    elif het_atm != " " and at_elem == "CD":
        atom_type = 22
    else:
        atom_type = _get_atom_type(resname, at_name, at_elem)
        if atom_type is None:
            warnings.warn(
                f"{at_name}:{resname} not in radii library.", BiopythonWarning
            )
            return 0.01
    return _atomic_radii[atom_type][typekey]


def _read_vertex_array(filename):
    """Read the vertex list into a NumPy array (PRIVATE)."""
    with open(filename) as fp:
        vertex_list = []
        for line in fp:
            sl = line.split()
            if len(sl) != 9:
                # skip header
                continue
            vl = [float(x) for x in sl[0:3]]
            vertex_list.append(vl)
    return np.array(vertex_list)


def get_surface(model, MSMS="msms"):
    """Represent molecular surface as a vertex list array.

    Return a NumPy array that represents the vertex list of the
    molecular surface.

    Arguments:
     - model - BioPython PDB model object (used to get atoms for input model)
     - MSMS - msms executable (used as argument to subprocess.call)

    """
    # Replace pdb_to_xyzr
    # Make x,y,z,radius file
    atom_list = Selection.unfold_entities(model, "A")

    xyz_tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(xyz_tmp, "w") as pdb_to_xyzr:
        for atom in atom_list:
            x, y, z = atom.coord
            radius = _get_atom_radius(atom, rtype="united")
            pdb_to_xyzr.write(f"{x:6.3f}\t{y:6.3f}\t{z:6.3f}\t{radius:1.2f}\n")

    # Make surface
    surface_tmp = tempfile.NamedTemporaryFile(delete=False).name
    msms_tmp = tempfile.NamedTemporaryFile(delete=False).name
    MSMS = MSMS + " -probe_radius 1.5 -if %s -of %s > " + msms_tmp
    make_surface = MSMS % (xyz_tmp, surface_tmp)
    subprocess.call(make_surface, shell=True)
    face_file = surface_tmp + ".face"
    surface_file = surface_tmp + ".vert"
    if not os.path.isfile(surface_file):
        raise RuntimeError(
            f"Failed to generate surface file using command:\n{make_surface}"
        )

    # Read surface vertices from vertex file
    surface = _read_vertex_array(surface_file)

    # Remove temporary files
    for fn in [xyz_tmp, surface_tmp, msms_tmp, face_file, surface_file]:
        try:
            os.remove(fn)
        except OSError:
            pass

    return surface


def min_dist(coord, surface):
    """Return minimum distance between coord and surface."""
    d = surface - coord
    d2 = np.sum(d * d, 1)
    return np.sqrt(min(d2))


def residue_depth(residue, surface):
    """Residue depth as average depth of all its atoms.

    Return average distance to surface for all atoms in a residue,
    ie. the residue depth.
    """
    atom_list = residue.get_unpacked_list()
    length = len(atom_list)
    d = 0
    for atom in atom_list:
        coord = atom.get_coord()
        d = d + min_dist(coord, surface)
    return d / length


def ca_depth(residue, surface):
    """Return CA depth."""
    if not residue.has_id("CA"):
        return None
    ca = residue["CA"]
    coord = ca.get_coord()
    return min_dist(coord, surface)


class ResidueDepth(AbstractPropertyMap):
    """Calculate residue and CA depth for all residues."""

    def __init__(self, model, msms_exec=None):
        """Initialize the class."""
        if msms_exec is None:
            msms_exec = "msms"

        depth_dict = {}
        depth_list = []
        depth_keys = []
        # get_residue
        residue_list = Selection.unfold_entities(model, "R")
        # make surface from PDB file using MSMS
        surface = get_surface(model, MSMS=msms_exec)
        # calculate rdepth for each residue
        for residue in residue_list:
            if not is_aa(residue):
                continue
            rd = residue_depth(residue, surface)
            ca_rd = ca_depth(residue, surface)
            # Get the key
            res_id = residue.get_id()
            chain_id = residue.get_parent().get_id()
            depth_dict[(chain_id, res_id)] = (rd, ca_rd)
            depth_list.append((residue, (rd, ca_rd)))
            depth_keys.append((chain_id, res_id))
            # Update xtra information
            residue.xtra["EXP_RD"] = rd
            residue.xtra["EXP_RD_CA"] = ca_rd
        AbstractPropertyMap.__init__(self, depth_dict, depth_keys, depth_list)
