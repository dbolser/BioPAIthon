BioPAIthon README file
======================

**BioPAIthon** is a fork of `Biopython <https://biopython.org>`_, the
long-running international association of developers of freely available
Python tools for computational molecular biology.

The fork exists to try one thing differently. BioPAIthon accepts contributions
from any capable contributor - human, computational intelligence, AI, or chimp
- and judges each patch on its merits rather than on the nature of its author.
See `AGENTS.md <AGENTS.md>`__ for what that means in practice.

Everything else is Biopython. The ``Bio`` package, the public API, the test
suite and the licence are all unchanged, ``import Bio`` works exactly as it
did, and we intend to keep tracking upstream.

Credit for essentially all of the code here belongs to the Biopython
contributors; the original copyright notices and licence remain in place.
BioPAIthon is not affiliated with or endorsed by the Biopython Project, so
please do not report problems with this fork to upstream.

This README is intended primarily for people interested in working with the
source code, from our repository on GitHub
https://github.com/dbolser/BioPAIthon

Upstream's user-centric documentation, `The Biopython Tutorial and Cookbook,
and API documentation <https://biopython.org/docs/latest/>`_, still describes
this fork accurately.

The `NEWS <NEWS.rst>`_ file summarises the changes in each release, alongside
the `DEPRECATED <DEPRECATED.rst>`_ file which notes API breakages.

This package is open source software made available under generous terms.
Please see the `LICENSE <LICENSE.rst>`_ file for further details.


Acknowledgements
================

BioPAIthon exists because the Biopython Project spent more than twenty years
building it. 373 people are named in `CONTRIB.rst <CONTRIB.rst>`_, and that
list is preserved here unchanged, as is the full release history in
`NEWS <NEWS.rst>`__ and every copyright notice in the source.

This fork changes a contribution policy. It does not change who wrote the
software, and it is not a criticism of the people who did. Upstream restricted
AI-assisted contributions for a reason we think is entirely legitimate:
protecting "good first issue" tickets as mentoring opportunities for new human
contributors, in the hope they stay. We disagree about the means, not the goal.
Reasonable projects can land in different places on this, and Biopython's
maintainers are under no obligation to agree with us.

Upstream Biopython is actively maintained and lives at https://biopython.org.
We track it, and we would be glad to see this fork become unnecessary.


Citation
========

If you use BioPAIthon in work contributing to a scientific publication, please
cite **Biopython**. Essentially all of the code here was written by the
Biopython contributors, and academic citation is how that work is credited:

Cock, P.J.A. et al. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

Several Biopython modules have their own publications, listed at
https://biopython.org - please cite those as well where they apply.

Do not cite BioPAIthon instead of Biopython. If you need to record which fork
you used, cite the paper above and mention the fork separately, as you would
any other software version detail. Nothing about this fork changes who did the
scientific work.


For the impatient
=================

BioPAIthon is not currently published on PyPI, so install it from source::

    pip install git+https://github.com/dbolser/BioPAIthon.git

Note that BioPAIthon installs the same ``Bio`` package that Biopython does, so
the two conflict and cannot be used side by side in one environment. Use a
virtual environment if you need both.

Upstream Biopython, by contrast, is on PyPI and has shipped pre-compiled
binary wheels for Linux, macOS and Windows since 1.70, so installing it is
quick and needs no compiler::

    pip install biopython

As a developer or potential contributor, you may wish to download, build and
install BioPAIthon yourself. This is described below.


Python Requirements
===================

We currently recommend using Python 3.13 from https://www.python.org

Biopython is currently supported and tested on the following Python
implementations:

- Python 3.10, 3.11, 3.12, 3.13 and 3.14 -- see https://www.python.org

- PyPy3.10 v7.3.17 -- or later, see https://www.pypy.org


Optional Dependencies
=====================

Biopython requires NumPy (see https://www.numpy.org) which will be installed
automatically if you install Biopython with pip (see below for compiling
Biopython yourself).

Depending on which parts of Biopython you plan to use, there are a number of
other optional Python dependencies, which can be installed later if needed:

- ReportLab, see https://www.reportlab.com/opensource/ (optional)
  This package is only used in ``Bio.Graphics``, so if you do not need this
  functionality, you will not need to install this package.

- matplotlib, see https://matplotlib.org/ (optional)
  ``Bio.Phylo`` uses this package to plot phylogenetic trees.

- networkx, see https://networkx.github.io/ (optional) and
  pygraphviz or pydot, see https://pygraphviz.github.io/ and
  https://code.google.com/p/pydot/ (optional)
  These packages are used for certain niche functions in ``Bio.Phylo``.

- rdflib, see https://github.com/RDFLib/rdflib (optional)
  This package is used in the CDAO parser under ``Bio.Phylo``.

- psycopg2, see https://initd.org/psycopg/ (optional) or
  PyGreSQL (pgdb), see https://www.pygresql.org/ (optional)
  These packages are used by ``BioSQL`` to access a PostgreSQL database.

- MySQL Connector/Python, see https://dev.mysql.com/downloads/connector/python/
  This package is used by ``BioSQL`` to access a MySQL database, and is
  supported on PyPy too.

- mysqlclient, see https://github.com/PyMySQL/mysqlclient-python (optional)
  This is a fork of the older MySQLdb and is used by ``BioSQL`` to access a
  MySQL database. It is supported by PyPy.

In addition there are a number of useful third party tools you may wish to
install such as standalone NCBI BLAST, EMBOSS or ClustalW.


Installation From Source
========================

BioPAIthon does not publish binary wheels, so installing it means compiling
the C extensions yourself. The following are required at compile time:

- Python including development header files like ``python.h``, which on Linux
  are often not installed by default (trying looking for and installing a
  package named ``python-dev`` or ``python-devel`` as well as the ``python``
  package).

- Appropriate C compiler for your version of Python, for example GCC on Linux,
  or MSVC on Windows. For Windows, you must install the 'Visual Studio Build Tools'
  and select the 'Desktop development with C++' workload. For macOS, use Apple's
  command line tools, which can be installed with the terminal command::

      xcode-select --install

  This will offer to install Apple's XCode development suite - you can, but it
  is not needed and takes a lot of disk space.

Then either download and decompress our source code, or fetch it using git.
Now change directory to the Biopython source code folder and run::

    pip install -e . --group dev
    cd Tests
    python run_tests.py

Substitute ``python`` with your specific version if required, for example
``python3``, or ``pypy3``.

To exclude tests that require an internet connection (and which may take a
long time), use the ``--offline`` option::

    cd Tests
    python run_tests.py --offline

Testing
=======

Biopython includes a suite of regression tests to check if everything is
running correctly. To run the tests, go to the biopython source code
directory and type::

    pip install -e . --group dev
    cd Tests
    python run_tests.py

If you want to skip the online tests (which is recommended when doing repeated
testing), use::

    cd Tests
    python run_tests.py --offline

Do not panic if you see messages warning of skipped tests::

    test_DocSQL ... skipping. Install MySQLdb if you want to use Bio.DocSQL.

This most likely means that a package is not installed.  You can
ignore this if it occurs in the tests for a module that you were not
planning on using.  If you did want to use that module, please install
the required dependency and re-run the tests.

Some of the tests may fail due to network issues, this is often down to
chance or a service outage. If the problem does not go away on
re-running the tests, you can use the ``--offline`` option.

There is more testing information in the Biopython Tutorial & Cookbook.


Experimental code
=================

Biopython 1.61 introduced a new warning, ``Bio.BiopythonExperimentalWarning``,
which is used to mark any experimental code included in the otherwise
stable Biopython releases. Such 'beta' level code is ready for wider
testing, but still likely to change, and should only be tried by early
adopters in order to give feedback via the biopython-dev mailing list.

We'd expect such experimental code to reach stable status within one or two
releases, at which point our normal policies about trying to preserve
backwards compatibility would apply.


Bugs
====

While we try to ship a robust package, bugs inevitably pop up.  If you are
having problems that might be caused by a bug in Biopython, it is possible
that it has already been identified. Update to the latest release if you are
not using it already, and retry. If the problem persists, please search our
bug database and our mailing lists to see if it has already been reported
(and hopefully fixed), and if not please do report the bug. We can't fix
problems we don't know about ;)

Issue tracker: https://github.com/dbolser/BioPAIthon/issues

Please report bugs in this fork there, not to upstream. If you can reproduce
the problem against upstream Biopython as well, then it is an upstream bug and
belongs at https://github.com/biopython/biopython/issues instead.

If you suspect the problem lies within a parser, it is likely that the data
format has changed and broken the parsing code.  (The text BLAST and GenBank
formats seem to be particularly fragile.)  Thus, the parsing code in
Biopython is sometimes updated faster than we can build Biopython releases.
You can get the most recent parser by pulling the relevant files (e.g. the
ones in ``Bio.SeqIO`` or ``Bio.Blast``) from our git repository. However, be
careful when doing this, because the code in github is not as well-tested
as released code, and may contain new dependencies.

In any bug report, please let us know:

1. Which operating system and hardware (32 bit or 64 bit) you are using
2. Python version
3. Biopython version (or git commit/date)
4. Traceback that occurs (the full error message)

And also ideally:

5. Example code that breaks
6. A data file that causes the problem


Contributing, Bug Reports
=========================

BioPAIthon accepts contributions from anyone, and anything, able to write a
good patch - human, computational intelligence, AI, or chimp. We are always
looking for help with code development, documentation writing, technical
administration, and whatever else comes up.

If you wish to contribute, please first read `AGENTS.md <AGENTS.md>`__, which
sets out the standard every contribution is held to regardless of its author,
and then `CONTRIBUTING.rst <CONTRIBUTING.rst>`_ for the practicalities.

Upstream Biopython is run by volunteers from all over the world, with many
types of backgrounds. Their web site is https://biopython.org and their
mailing lists are at https://biopython.org/wiki/Mailing_lists


Distribution Structure
======================

- ``README.rst``  -- This file.
- ``NEWS.rst``    -- Release notes and news.
- ``LICENSE.rst`` -- What you can do with the code.
- ``CONTRIB.rst`` -- An (incomplete) list of people who helped Biopython in
  one way or another.
- ``CONTRIBUTING.rst`` -- The practicalities of contributing.
- ``AGENTS.md``   -- Who may contribute (everyone) and the standard every
  contribution is held to. ``CLAUDE.md`` points here.
- ``DEPRECATED.rst`` -- Contains information about modules in Biopython that
  were removed or no longer recommended for use, and how to update code that
  uses those modules.
- ``MANIFEST.in`` -- Configures which files to include in releases.
- ``pyproject.toml`` -- Project metadata and build configuration.
- ``Bio/``        -- The main code base code.
- ``BioSQL/``     -- Code for using Biopython with BioSQL databases.
- ``Doc/``        -- Documentation.
- ``Scripts/``    -- Miscellaneous, possibly useful, standalone scripts.
- ``Tests/``      -- Regression testing code including sample data files.
