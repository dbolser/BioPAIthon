Contributing to BioPAIthon
==========================

We welcome pull requests to fix bugs or add new features, from anyone - and
anything - able to write a good patch. See `AGENTS.md <AGENTS.md>`__ for what
that means in practice. This file is primarily concerned with the
practicalities.

BioPAIthon is a fork of Biopython, and most of the guidance below is inherited
from that project. Upstream's general overview is still worth reading:
http://biopython.org/wiki/Contributing


Licensing
---------

BioPAIthon inherits Biopython's licensing unchanged. Contributions are offered
under *both* the original "Biopython License Agreement" and the more widely
used "3-Clause BSD License" (see our LICENSE file for more details). Opening a
pull request against this repository is taken as offering your work under
both, on the usual inbound-equals-outbound basis; you do not need to restate it
in each commit or pull request.

Nothing is assigned to anyone by this. You keep the copyright in what you
write, and the project gains no rights you have not equally granted to
everybody else. This is not a contributor licence agreement and there is
nothing to sign.

Both licences are named because upstream Biopython is part-way through
offering the whole library under the pair, and a contribution offered under
only one of them could not follow.

Do not remove or alter existing copyright notices. Essentially all of this
code was written by the Biopython contributors, and that attribution stays.


Tagged issues
-------------

In general, unless there is an open pull request associated with an issue,
or the issue discussion makes it clear that someone is working on it, any
issue can be worked on.

We have been using the `good first issue
<https://github.com/dbolser/BioPAIthon/issues?q=is%3Aissue%20state%3Aopen%20label%3A%22good%20first%20issue%22>`_
label specifically for beginners. These are minor issues where the team knows
what is needed, and could easily do it themselves (albeit some are tedious),
but we have chosen instead to use these as teachable moments to help potential
new contributors get involved with a little mentoring as needed. This is with
the hope they become long term open source contributors.

Use whatever tools help you do that well, AI included. We ask only that you
understand the change you are proposing, can explain why it is correct, and
have run the tests. Those conditions apply identically to every contributor;
see `AGENTS.md <AGENTS.md>`__.

The `help wanted
<https://github.com/dbolser/BioPAIthon/issues?q=is%3Aissue%20state%3Aopen%20label%3A%22help%20wanted%22>`_
label is typically used for issues where none of the regular contributors have the
expertise or setup to work on an issue (eg specific platforms, or restricted datasets).
The `needs info
<https://github.com/dbolser/BioPAIthon/issues?q=is%3Aissue%20state%3Aopen%20label%3A%22Needs%20info%22>`_
label is for unresloved issues where we lack actionable information. In some cases
these may be best closed, but perhaps another person will run into the same problem
and can help solve it.

Git Usage
---------

We have a git introduction online at http://biopython.org/wiki/GitUsage

If you are planning to make a pull request, start by creating a new branch
with a short but descriptive name (rather than using your main branch).

Failing to follow our pull result template places additional burden on the
maintainers in reviewing your changes and if they can be accepted. We may
simply close it without further discussion.

Coding Conventions
------------------

Biopython tries to follow the coding conventions laid out in PEP8 and PEP257,
with the notable exception of existing module names which are not lower case.

- http://www.python.org/dev/peps/pep-0008/
- http://www.python.org/dev/peps/pep-0257/

For docstrings (Python's in-code documentation), in addition to PEP257 we are
using reStructuredText (RST) markup language which allows basic formatting
like *italics* and **bold** once rendered into HTML webpages for our online
API documentation. We also use reStructuredText for files like ``README.rst``.

To facilitate style checking, we make use of ``pre-commit`` to automatically
run ``black``, ``ruff``, ``flake8`` and ``mypy`` on the Python code, and
``rstcheck`` and ``doc8`` on the RST files, among other checks. We also use
continuous integration on GitHub to run these checks, but we strongly suggest
you install and run ``pre-commit`` in your local machine (see the Setup
section of `AGENTS.md <AGENTS.md>`__).


Testing
-------

Any new feature or functionality will not be accepted without tests. Likewise
for any bug fix, we encourage including an additional test. See the testing
chapter in the Biopython Tutorial for more information on our test framework:
http://biopython.org/DIST/docs/tutorial/Tutorial.html


Local Development and Testing
-----------------------------

Please run the full test suite locally before submitting a pull request, and
run ``pre-commit`` so the style checks pass before CI sees your change. How to
install BioPAIthon for development, enable the ``pre-commit`` hooks, and run
the test suite is described once, in the Setup and Testing sections of
`AGENTS.md <AGENTS.md>`__; the README carries the same commands for users.

Continuous Integration
----------------------

Once you submit your pull request on GitHub, several automated tests should be
run via continuous integration services, and their results reported on the pull
request. These will run most of the Biopython tests (although not with all the
optional dependencies included), plus also style checks using ``pre-commit``
(also used for git pre-commit checks, see above).

The continuous integration tests collect test coverage information via
CodeCov: https://codecov.io/github/dbolser/BioPAIthon/

**The continuous integration checks must pass before your pull request will be
merged.** Ideally the CodeCov checks will also pass, but we currently do not insist on
this when reviewing pull requests.

Contributing to the Biopython Tutorial
--------------------------------------

For instructions, see `the Biopython Tutorial README <Doc/README.rst>`_
