.. _`chapter:testing`:

The Biopython testing framework
===============================

Biopython has a regression test suite written with
`unittest <https://docs.python.org/3/library/unittest.html>`__, the
standard unit testing framework for Python, and run by
`pytest <https://docs.pytest.org/>`__. Providing comprehensive
tests for modules is one of the most important aspects of making sure
that the Biopython code is as bug-free as possible before going out. It
also tends to be one of the most undervalued aspects of contributing.
This chapter is designed to make running the Biopython tests and writing
good test code as easy as possible. Ideally, every module that goes into
Biopython should have a test (and should also have documentation!). All
our developers, and anyone installing Biopython from source, are
strongly encouraged to run the unit tests.

Running the tests
-----------------

When you download the Biopython source code, or check it out from our
source code repository, you should find a subdirectory called ``Tests``.
This contains lots of individual scripts named ``test_XXX.py``, the
script ``run_tests.py``, the pytest configuration file ``conftest.py``,
and lots of other subdirectories which contain input files for the test
suite.

The tests need pytest 9 or later. Installing Biopython from its source
folder with the ``test`` extra brings pytest and the optional packages
the tests use:

.. code:: console

   $ pip install -e ".[test]"

As part of building and installing Biopython you will typically run the
full test suite at the command line from the Biopython source top level
directory using the following:

.. code:: console

   $ cd Tests
   $ python run_tests.py

You’ll often want to run just some of the tests, and this is done like
this:

.. code:: console

   $ python run_tests.py test_SeqIO.py test_AlignIO.py

When giving the list of tests, the ``.py`` extension is optional, so you
can also just type:

.. code:: console

   $ python run_tests.py test_SeqIO test_AlignIO

To run the docstring tests (see section :ref:`sec:doctest`), you can
use

.. code:: console

   $ python run_tests.py doctest

You can also skip any tests which have been setup with an explicit
online component by adding ``--offline``, e.g.

.. code:: console

   $ python run_tests.py --offline

``--offline`` also makes any ``socket.connect`` call to another machine
raise an error. ``urllib`` connects that way, so a test which goes
online by mistake will usually fail rather than quietly use the
network. This is a safety net, not a sandbox: DNS lookups and
subprocesses, for example, are not blocked.

By default, ``run_tests.py`` runs all tests, including the docstring
tests.

If an individual test is failing, you can also try running it directly
(e.g. ``python test_SeqIO.py``), which may give you more information.
This works from any directory, so from the top level source directory
it is ``python Tests/test_SeqIO.py``.

Tests based on Python’s standard ``unittest`` framework will
``import unittest`` and then define ``unittest.TestCase`` classes, each
with one or more sub-tests as methods starting with ``test_`` which
check some specific aspect of the code.

Running the tests with pytest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``run_tests.py`` translates its arguments into a pytest command line
and runs pytest, so you can also call pytest yourself and use any of
its options. In the ``Tests`` directory,

.. code:: console

   $ python -m pytest --offline

is the same as ``python run_tests.py --offline``. pytest wants test
files with their ``.py`` extension, and picks out tests whose names
match an expression with ``-k``:

.. code:: console

   $ python -m pytest --offline test_Seq_objs.py -k translate

Add ``-x`` to stop at the first failure, or ``--lf`` to rerun only the
tests which failed last time:

.. code:: console

   $ python -m pytest --offline -x
   $ python -m pytest --offline --lf

pytest names each test by a node ID, such as
``test_Seq_objs.py::StringMethodTests::test_str_count``, which it
prints when the test fails. Give one to run just that test, or leave
off the method name to run the whole ``TestCase`` class:

.. code:: console

   $ python -m pytest test_Seq_objs.py::StringMethodTests::test_str_count
   $ python -m pytest test_Seq_objs.py::StringMethodTests

The docstring tests are collected through the file
``test_docstrings.py``, one test for each docstring, so you can run all
of them, or just those of one module:

.. code:: console

   $ python -m pytest --offline test_docstrings.py
   $ python -m pytest test_docstrings.py::Bio.Seq

With the editable install shown above, pytest also works from any other
directory, given the path to the test files, such as from the top level
source directory:

.. code:: console

   $ python -m pytest --offline Tests/test_Seq_objs.py

The tests run in the directory pytest was started in. Only the docstring
tests change into ``Tests``, as their examples name files relative to
it.

Running the tests using Tox
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like most Python projects, you can also use
`Tox <https://tox.readthedocs.org/en/latest/>`__ to run the tests on
multiple Python versions, provided they are already installed in your
system.

We do not provide the configuration ``tox.ini`` file in our code base
because of difficulties pinning down user-specific settings (e.g.
executable names of the Python versions). You may also only be
interested in testing Biopython only against a subset of the Python
versions that we support.

If you are interested in using Tox, you could start with the example
``tox.ini`` shown below:

.. code:: text

   [tox]
   env_list = py310, py311, py312, py313, py314, pypy3

   [testenv]
   extras = !pypy3: test
   deps = pypy3: pytest>=9
   change_dir = Tests
   commands = python -m pytest --offline {posargs}

Using the template above, executing ``tox`` will test your Biopython
code against Python 3.10 to 3.14 and PyPy. It assumes that those
Pythons’ executables are named ``python3.12`` for Python 3.12, and so
on. Each environment installs Biopython with the ``test`` extra, except
PyPy, which gets just pytest because SciPy has no wheels for PyPy.
Anything after ``--`` is passed on to pytest, so
``tox -e py312 -- -x test_SeqIO.py`` runs one module under Python 3.12
and stops at its first failure.

Writing tests
-------------

Let’s say you want to write some tests for a module called ``Biospam``.
This can be a module you wrote, or an existing module that doesn’t have
any tests yet. In the examples below, we assume that ``Biospam`` is a
module that does simple math.

Each Biopython test consists of a script containing the test itself, and
optionally a directory with input files used by the test:

#. ``test_Biospam.py`` – The actual test code for your module.

#. ``Biospam`` [optional]– A directory where any necessary input files
   will be located. If you have any output files that should be manually
   reviewed, output them here (but this is discouraged) to prevent
   clogging up the main Tests directory. In general, use a temporary
   file/folder.

Build the paths to your input files from ``support.DATA``, the
``Tests`` directory (defined in ``Tests/support.py``), rather than
relative to the current directory:

.. code:: python

   import support

   with open(support.DATA / "Biospam" / "example.txt") as handle:
       ...

That way your test passes however it is run, not only from inside
``Tests``. A relative path such as ``"Biospam/example.txt"`` fails
whenever the tests are run from another directory, as our continuous
integration does.

Any script with a ``test_`` prefix in the ``Tests`` directory will be
found and run by pytest. Only the tests in ``unittest.TestCase``
subclasses are collected, so write your tests as those: a plain
pytest-style ``def test_...()`` function is ignored. Below, we show an
example test script ``test_Biospam.py``. If you put this script in the
Biopython ``Tests`` directory, then pytest will find it and execute the
tests contained in it:

.. code:: console

   $ python run_tests.py --offline
   ============================= test session starts ==============================
   platform linux -- Python 3.12.3, pytest-9.1.1, pluggy-1.6.0
   rootdir: /home/user/biopython
   configfile: pyproject.toml
   collected 4216 items / 22 skipped

   test_Ace.py ..........                                                   [  0%]
   test_Affy.py ......                                                      [  0%]
   test_AlignIO.py .....................................................    [  1%]
   ...
   test_BioSQL_sqlite3.py ................................................. [ 16%]
   ........................                                                 [ 17%]
   test_Biospam.py ....                                                     [ 17%]
   ...
   test_docstrings.py ..................................................... [ 82%]
   ...
   ========================== modules skipped at import ===========================
   ...
   test_BioSQL_sqlite3_online -- internet not available
   ...
   ========== 4180 passed, 58 skipped, 29 warnings in 326.42s (0:05:26) ===========

Each dot is a test which passed, and each ``s`` one which was skipped.
Modules which skip as a whole, for example because they need the
internet or an optional package which is not installed, are listed at
the end with the reason.

Writing a test using ``unittest``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``unittest``-framework has been included with Python since version
2.1, and is documented in the Python Library Reference (which I know you
are keeping under your pillow, as recommended). There is also `online
documentation for
unittest <https://docs.python.org/3/library/unittest.html>`__. If you
are familiar with the ``unittest`` system (or something similar like the
nose test framework), you shouldn’t have any trouble. You may find
looking at the existing examples within Biopython helpful too.

Here’s a minimal ``unittest``-style test script for ``Biospam``, which
you can copy and paste to get started:

.. code:: python

   import unittest
   from Bio import Biospam


   class BiospamTestAddition(unittest.TestCase):
       def test_addition1(self):
           result = Biospam.addition(2, 3)
           self.assertEqual(result, 5)

       def test_addition2(self):
           result = Biospam.addition(9, -1)
           self.assertEqual(result, 8)


   class BiospamTestDivision(unittest.TestCase):
       def test_division1(self):
           result = Biospam.division(3.0, 2.0)
           self.assertAlmostEqual(result, 1.5)

       def test_division2(self):
           result = Biospam.division(10.0, -2.0)
           self.assertAlmostEqual(result, -5.0)


   if __name__ == "__main__":
       runner = unittest.TextTestRunner(verbosity=2)
       unittest.main(testRunner=runner)

In the division tests, we use ``assertAlmostEqual`` instead of
``assertEqual`` to avoid tests failing due to roundoff errors; see the
``unittest`` chapter in the Python documentation for details and for
other functionality available in ``unittest`` (`online
reference <https://docs.python.org/3/library/unittest.html>`__).

These are the key points of ``unittest``-based tests:

-  Test cases are stored in classes that derive from
   ``unittest.TestCase`` and cover one basic aspect of your code

-  You can use methods ``setUp`` and ``tearDown`` for any repeated code
   which should be run before and after each test method. For example,
   the ``setUp`` method might be used to create an instance of the
   object you are testing, or open a file handle. The ``tearDown``
   should do any “tidying up”, for example closing the file handle.

-  The tests are prefixed with ``test_`` and each test should cover one
   specific part of what you are trying to test. You can have as many
   tests as you want in a class.

-  At the end of the test script, you can use

   .. code:: python

      if __name__ == "__main__":
          runner = unittest.TextTestRunner(verbosity=2)
          unittest.main(testRunner=runner)

   to execute the tests when the script is run by itself (rather than
   collected by pytest). If you run this script, then you’ll
   see something like the following:

   .. code:: console

      $ python test_Biospam.py
      test_addition1 (__main__.BiospamTestAddition.test_addition1) ... ok
      test_addition2 (__main__.BiospamTestAddition.test_addition2) ... ok
      test_division1 (__main__.BiospamTestDivision.test_division1) ... ok
      test_division2 (__main__.BiospamTestDivision.test_division2) ... ok

      ----------------------------------------------------------------------
      Ran 4 tests in 0.059s

      OK

-  To indicate more clearly what each test is doing, you can add
   docstrings to each test. These are shown when running the tests,
   which can be useful information if a test is failing.

   .. code:: python

      import unittest
      from Bio import Biospam


      class BiospamTestAddition(unittest.TestCase):
          def test_addition1(self):
              """An addition test"""
              result = Biospam.addition(2, 3)
              self.assertEqual(result, 5)

          def test_addition2(self):
              """A second addition test"""
              result = Biospam.addition(9, -1)
              self.assertEqual(result, 8)


      class BiospamTestDivision(unittest.TestCase):
          def test_division1(self):
              """Now let's check division"""
              result = Biospam.division(3.0, 2.0)
              self.assertAlmostEqual(result, 1.5)

          def test_division2(self):
              """A second division test"""
              result = Biospam.division(10.0, -2.0)
              self.assertAlmostEqual(result, -5.0)


      if __name__ == "__main__":
          runner = unittest.TextTestRunner(verbosity=2)
          unittest.main(testRunner=runner)

   Running the script will now show you:

   .. code:: console

      $ python test_Biospam.py
      test_addition1 (__main__.BiospamTestAddition.test_addition1)
      An addition test ... ok
      test_addition2 (__main__.BiospamTestAddition.test_addition2)
      A second addition test ... ok
      test_division1 (__main__.BiospamTestDivision.test_division1)
      Now let's check division ... ok
      test_division2 (__main__.BiospamTestDivision.test_division2)
      A second division test ... ok

      ----------------------------------------------------------------------
      Ran 4 tests in 0.001s

      OK

If your module contains docstring tests (see section
:ref:`sec:doctest`), you *may* want to include those in the tests to
be run. You can do so as follows by adding ``import doctest`` and
``import sys`` to the imports and modifying the code under
``if __name__ == "__main__":`` to look like this:

.. code:: python

   if __name__ == "__main__":
       unittest_suite = unittest.TestLoader().loadTestsFromName("test_Biospam")
       doctest_suite = doctest.DocTestSuite(Biospam)
       suite = unittest.TestSuite((unittest_suite, doctest_suite))
       runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
       runner.run(suite)

This is only relevant if you want to run the docstring tests when you
execute ``python test_Biospam.py`` if it has some complex run-time
dependency checking.

In general you do not need this, as the test suite runs the docstring
tests of every Biopython module itself, as explained below.

.. _`sec:doctest`:

Writing doctests
----------------

Python modules, classes and functions support built-in documentation
using docstrings. The `doctest
framework <https://docs.python.org/3/library/doctest.html>`__ (included
with Python) allows the developer to embed working examples in the
docstrings, and have these examples automatically tested.

Currently only part of Biopython includes doctests. The test suite
takes care of running them: ``Tests/conftest.py`` imports every module
in the ``Bio`` and ``BioSQL`` packages and collects each docstring
which has examples as one test, under the file ``test_docstrings.py``.
The tests are named after the module and then the object the docstring
belongs to, such as
``test_docstrings.py::Bio.Seq::Bio.Seq.translate``. So, if you’ve added
some doctests to the docstrings in a Biopython module, they will be run
without any change to the test suite. If an example fails, pytest shows
its expected and actual output.

The exceptions are named in two lists near the top of
``Tests/conftest.py``.
``EXCLUDE_DOCTEST_MODULES`` holds modules whose doctests are never run,
such as the stubs left behind for modules removed from Biopython, which
deliberately raise an ``ImportError`` naming their replacement.
``ONLINE_DOCTEST_MODULES`` holds modules whose doctests need internet
access, such as ``Bio.Entrez``. These are only left out when you use
``--offline``, so add your module there if its doctests go online.

Modules needing an optional dependency, such as ReportLab for
``Bio.Graphics``, are not listed. Instead such a module raises
``MissingPythonDependencyError`` (from ``Bio``) when it is imported
without the dependency, and the test suite reports it as skipped rather
than failed, in a "modules skipped at import" section at the end of the
run. A ``test_XXX.py`` script can skip itself the same way, by raising
``MissingExternalDependencyError`` while it is being imported. Only an
error raised during the import counts: the same exceptions raised later,
from the code being tested, are failures. The file
``Tests/expected_skips.txt`` lists the modules which may skip, one per
line, as ``test_XXX`` for a test script or ``Bio.XXX`` for a module's
doctests. ``python run_tests.py --check-skips`` (or
``python -m pytest --check-skips``) fails if any other module skips, so
add a line there if your module has a legitimate new reason to skip.

Note that we regard doctests primarily as documentation, so you should
stick to typical usage. Generally complicated examples dealing with
error conditions and the like would be best left to a dedicated unit
test.

Note that if you want to write doctests involving file parsing, defining
the file location complicates matters. Ideally use relative paths
assuming the code will be run from the ``Tests`` directory, where the
test suite runs each docstring's examples; see the ``Bio.SeqIO``
doctests for an example of this.

To run the docstring tests only, use

.. code:: console

   $ python run_tests.py doctest

Note that the doctest system is fragile and care is needed to ensure
your output will match on all the different versions of Python that
Biopython supports (e.g. differences in floating point numbers).

.. _`sec:doctest-tutorial`:

Writing doctests in the Tutorial
--------------------------------

This Tutorial you are reading has a lot of code snippets, which are
often formatted like a doctest. We have our own system in file
``test_Tutorial.py`` to allow tagging code snippets in the Tutorial
source to be run as Python doctests. This works by adding special
``.. doctest`` comment lines before each Python Console (pycon) block,
e.g.

.. code:: rst

   .. doctest

   .. code:: pycon

      >>> from Bio.Seq import Seq
      >>> s = Seq("ACGT")
      >>> len(s)
      4

Often code examples are not self-contained, but continue from the
previous Python block. Here we use the magic comment ``.. cont-doctest``
as shown here:

.. code:: rst

   .. cont-doctest

   .. code:: pycon

      >>> s == "ACGT"
      True

The special ``.. doctest`` comment line can take a working directory
(relative to the ``Doc/`` folder) to use if you have any example data
files, e.g. ``.. doctest examples`` will use the ``Doc/examples`` folder,
while ``.. doctest ../Tests/GenBank`` will use the ``Tests/GenBank``
folder.

After the directory argument, you can specify any Python dependencies
which must be present in order to run the test by adding ``lib:XXX`` to
indicate ``import XXX`` must work, e.g. ``.. doctest examples lib:numpy``

You can run the Tutorial doctests via:

.. code:: console

   $ python test_Tutorial.py

or:

.. code:: console

   $ python run_tests.py test_Tutorial.py
