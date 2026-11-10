# Copyright 2014 by Kevin Wu.
# Revisions copyright 2014 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Provides code to access the REST-style KEGG online API.

This module aims to make the KEGG online REST-style API easier to use. See:
https://www.kegg.jp/kegg/rest/keggapi.html

The KEGG REST-style API provides simple access to a range of KEGG databases.
This works using simple URLs (which this module will construct for you),
with any errors indicated via HTTP error levels.

The functionality is somewhat similar to Biopython's Bio.TogoWS and Bio.Entrez
modules.

Currently KEGG does not provide any usage guidelines (unlike the NCBI whose
requirements are reasonably clear). To avoid risking overloading the service,
Biopython will only allow three calls per second.

References:
Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes.
Nucleic Acids Res. 28, 29-34 (2000).

"""

import io
import re
import time
from urllib.request import urlopen

from Bio._utils import function_with_previous

# An unprefixed KEGG compound, drug or glycan entry, such as C00001
_CHEM_ENTRY = re.compile(r"[CDG]\d{5}")


@function_with_previous
def _q(op, arg1, arg2=None, arg3=None):
    delay = 0.333333333  # one third of a second
    current = time.time()
    wait = _q.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _q.previous = current + wait
    else:
        _q.previous = current

    URL = "https://rest.kegg.jp/%s"
    if arg2 and arg3:
        args = f"{op}/{arg1}/{arg2}/{arg3}"
    elif arg2:
        args = f"{op}/{arg1}/{arg2}"
    else:
        args = f"{op}/{arg1}"
    resp = urlopen(URL % (args))

    if "image" == arg2:
        return resp

    handle = io.TextIOWrapper(resp, encoding="UTF-8")
    handle.url = resp.url
    return handle


_q.previous = 0


# https://www.kegg.jp/kegg/rest/keggapi.html
def kegg_info(database):
    """KEGG info - Displays the current statistics of a given database.

    db - database or organism (string)

    The argument db can be a KEGG database name (e.g. 'pathway' or its
    official abbreviation, 'path'), or a KEGG organism code or T number
    (e.g. 'hsa' or 'T01001' for human).

    A valid list of organism codes and their T numbers can be obtained
    via kegg_list('genome') or https://rest.kegg.jp/list/genome

    """
    # TODO - return a string (rather than the handle?)
    # TODO - cache and validate the organism code / T numbers?
    # TODO - can we parse the somewhat formatted output?
    #
    # https://rest.kegg.jp/info/<database>
    #
    # <database> = kegg | pathway | brite | module | ko | genes | <org> | ag |
    #              vg | vp | genome | vtax | vgenome | compound | glycan |
    #              reaction | rclass | rmodule | enzyme | network | ntmap |
    #              variant | disease | drug | dgroup
    # <org> = KEGG organism code or T number
    return _q("info", database)


def kegg_list(database, org=None, option=None):
    """KEGG list - Entry list for database, or specified database entries.

    db - database or organism (string), or database entries (string, or
    list of strings)
    org - optional organism (string), see below.
    option - optional subset of the brite or genome database (string),
    see below.

    For the pathway database the optional organism can be used to
    restrict the results, as in ``kegg_list("pathway", "hsa")``.

    For the brite database, option can be "br", "jp", "ko" or a KEGG
    organism code; for the genome database, it can be a KEGG organism
    group name such as "animals", or the NCBI taxonomy ID of a phylum,
    class, order, family, genus or species, such as "562".

    At most 100 database entries can be listed at once, whether given as
    a list or joined with "+" in a string: KEGG silently drops any beyond
    the first 100.

    """
    # TODO - split into two functions (dbentries seems separate)?
    #
    # https://rest.kegg.jp/list/pathway/<org>
    #
    # <org> = KEGG organism code
    #
    # (KEGG no longer offers /list/module/<org>, nor /list/organism.)
    if database == "pathway" and org and not option:
        return _q("list", database, org)

    # https://rest.kegg.jp/list/brite/<option>
    #
    # <option> = br | jp | ko | <org>
    #
    # https://rest.kegg.jp/list/genome/<option>
    #
    # <option> = <group_name> | <rank_id>
    # <group_name> = KEGG organism group name
    # <rank_id> = Taxonomy ID for phylum, class, order, family, genus and species
    if database in ("brite", "genome") and option and not org:
        return _q("list", database, option)

    if org and database in ("brite", "genome") and not option:
        raise ValueError(
            f"Invalid org arg for kegg list request: use "
            f"kegg_list({database!r}, option={org!r}) for /list/{database}/{org}."
        )
    if org:
        raise ValueError("Invalid database arg for kegg list request.")
    if option:
        raise ValueError("Invalid option arg for kegg list request.")

    # https://rest.kegg.jp/list/<database>
    #
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup
    # <org> = KEGG organism code or T number
    #
    # https://rest.kegg.jp/list/<dbentries>
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup
    # <org> = KEGG organism code or T number
    if isinstance(database, list):
        database = "+".join(database)
    if database.count("+") >= 100:
        raise ValueError("Maximum number of databases is 100 for kegg list query")
    return _q("list", database)


def kegg_find(database, query, option=None):
    """KEGG find - Data search.

    Finds entries with matching query keywords or other query data in
    a given database.

    db - database or organism (string)
    query - search terms (string)
    option - search option (string), see below.

    For the compound and drug database, set option to the string 'formula',
    'exact_mass', 'mol_weight' or 'nop' to search on that field only. The
    chemical formula search is a partial match irrespective of the order
    of atoms given. The exact mass (or molecular weight) is checked by
    rounding off to the same decimal place as the query data. A range of
    values may also be specified with the minus(-) sign.

    """
    # TODO - return list of tuples?
    #
    # https://rest.kegg.jp/find/<database>/<query>/<option>
    #
    # <database> = compound | drug
    # <option> = formula | exact_mass | mol_weight | nop
    if database in ["compound", "drug"] and option in [
        "formula",
        "exact_mass",
        "mol_weight",
        "nop",
    ]:
        resp = _q("find", database, query, option)
    elif option:
        raise ValueError("Invalid option arg for kegg find request.")

    # https://rest.kegg.jp/find/<database>/<query>
    #
    #
    # <database> = pathway | brite | module | ko | genes | <org> | ag | vg |
    #              vp | genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | enzyme | network | ntmap | variant | disease |
    #              drug | dgroup
    # <org> = KEGG organism code or T number
    else:
        if isinstance(query, list):
            query = "+".join(query)
        resp = _q("find", database, query)

    return resp


def kegg_get(dbentries, option=None):
    """KEGG get - Data retrieval.

    dbentries - Identifiers (single string, or list of strings), see below.
    option - One of "aaseq", "ntseq", "mol", "kcf", "image", "conf",
    "kgml", "json" (string)

    The input is limited up to 10 entries.
    The input is limited to one pathway entry with the image or kgml option.
    The input is limited to one compound/glycan/drug entry with the image option.

    Returns a handle.
    """
    if isinstance(dbentries, list) and len(dbentries) <= 10:
        dbentries = "+".join(dbentries)
    elif isinstance(dbentries, list) and len(dbentries) > 10:
        raise ValueError("Maximum number of dbentries is 10 for kegg get query")

    # https://rest.kegg.jp/get/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | disease_ja | drug_ja | dgroup_ja |
    #              compound_ja
    # <org> = KEGG organism code or T number
    #
    # <option> = aaseq | ntseq | mol | kcf | image | conf | kgml | json
    if option in ["aaseq", "ntseq", "mol", "kcf", "image", "conf", "kgml", "json"]:
        resp = _q("get", dbentries, option)
    elif option:
        raise ValueError("Invalid option arg for kegg get request.")
    else:
        resp = _q("get", dbentries)

    return resp


def kegg_conv(target_db, source_db, option=None):
    """KEGG conv - convert KEGG identifiers to/from outside identifiers.

    Arguments:
     - target_db - Target database
     - source_db - Source database, or database entries (string, or list
       of strings)
     - option - No longer supported: KEGG withdrew the "turtle" and
       "n-triple" output of conv, so any option raises ValueError.

    Convert a whole database, as in ``kegg_conv("eco", "ncbi-geneid")``,
    or given entries, as in ``kegg_conv("ncbi-proteinid",
    "hsa:10458+ece:Z5100")`` or ``kegg_conv("eco", "ncbi-geneid:948364")``.
    The outside databases are "ncbi-geneid", "ncbi-proteinid" and "uniprot"
    for genes, and "pubchem" and "chebi" for chemical substances; KEGG no
    longer offers "ncbi-gi".

    """
    # https://rest.kegg.jp/conv/<target_db>/<source_db>
    #
    # (<target_db> <source_db>) = (<kegg_db> <outside_db>) |
    #                             (<outside_db> <kegg_db>)
    #
    # For gene identifiers:
    # <kegg_db> = <org>
    # <org> = KEGG organism code or T number
    # <outside_db> = ncbi-geneid | ncbi-proteinid | uniprot
    #
    # For chemical substance identifiers:
    # <kegg_db> = drug | compound | glycan
    # <outside_db> = pubchem | chebi
    #
    # https://rest.kegg.jp/conv/<target_db>/<dbentries>
    #
    # For gene identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = <org> | genes | ncbi-geneid | ncbi-proteinid | uniprot
    # <org> = KEGG organism code or T number
    #
    # For chemical substance identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = drug | compound | glycan | pubchem | chebi
    #
    # KEGG used to offer an RDF <option> (turtle | n-triple) on both forms,
    # and ncbi-gi as an outside database; both are gone.
    if option:
        raise ValueError(
            "Invalid option arg for kegg conv request: "
            "KEGG no longer offers RDF output for conv."
        )

    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    entries = source_db.split("+") if isinstance(source_db, str) else []
    if target_db == "ncbi-gi" or any(
        entry == "ncbi-gi" or entry.startswith("ncbi-gi:") for entry in entries
    ):
        # KEGG answers the database form with HTTP 400, and the entries
        # form, even mixed with other entries, with an empty result.
        raise ValueError(
            "KEGG no longer offers ncbi-gi for kegg conv requests; "
            "use ncbi-proteinid instead."
        )

    gene_outside_db = ("ncbi-geneid", "ncbi-proteinid", "uniprot")
    chem_kegg_db = ("drug", "compound", "glycan")
    chem_outside_db = ("pubchem", "chebi")
    if (
        # Either form, converting to an outside gene database, but not
        # from a chemical database (HTTP 400)
        (
            target_db in gene_outside_db
            and source_db not in chem_kegg_db + chem_outside_db
        )
        # Converting to an outside chemical database, from a KEGG chemical
        # database, or from entries such as C00001, D00001 or cpd:C00001
        # (an organism gets HTTP 400, or an empty result for a T number)
        or (
            target_db in chem_outside_db
            and (
                source_db in chem_kegg_db
                or all(
                    ":" in entry or _CHEM_ENTRY.fullmatch(entry) for entry in entries
                )
            )
        )
        # Database form, from an outside gene database
        or source_db in gene_outside_db
        or (target_db in chem_kegg_db and source_db in chem_outside_db)
        # Entries form, such as hsa:10458 or ncbi-geneid:948364, which can
        # be converted to any database, including an organism or genes
        or (isinstance(source_db, str) and ":" in source_db)
    ):
        return _q("conv", target_db, source_db)

    raise ValueError("Bad argument target_db or source_db for kegg conv request.")


def kegg_link(target_db, source_db, option=None):
    """KEGG link - find related entries by using database cross-references.

    target_db - Target database
    source_db_or_dbentries - source database
    option - Can be "species", "genus", "family", "order", "class", "phylum".
    If getting data with RDF, can be "turtle" or "n-triple" (string).
    """
    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    # https://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = drug | atc | jtc
    #
    # <option> = turtle | n-triple
    # https://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries of the following <database>
    # <database> = drug | atc | jtc
    #
    # <option> = turtle | n-triple
    if option in ["turtle", "n-triple"] and target_db in ["drug", "atc", "jtc"]:
        return _q("link", target_db, source_db, option)

    # https://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | <outside_db>
    # <outside_db> = pubmed | taxonomy | atc | jtc | ndc | yk
    #
    # <option> = species | genus | family | order | class | phylum
    # https://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | <outside_db>
    # <outside_db> = pubmed | taxonomy | atc | jtc | ndc | yk
    #
    # <option> = species | genus | family | order | class | phylum
    if option in ["species", "genus", "family", "order", "class", "phylum"]:
        return _q("link", target_db, source_db, option)

    if not option:
        return _q("link", target_db, source_db)

    raise ValueError("Invalid option arg for kegg link request.")
