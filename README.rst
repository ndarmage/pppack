PPPACK
======

:Author: Daniele Tomatis
:Date: 02/05/2017

Piece-wise Polynomial Package for multivariate interpolation of real-valued scalar functions.

This package offers representations of multi-dimensional functions using piece-wise polynomials. Python handles memory management, whereas external FORTRAN libraries perform all the numerics in order to ensure high computational performances with big amount of data.

The main FORTRAN library is an enhanced version of pppack.f90, whose original version by John Burkardt is available on `netlib <http://www.netlib.org/pppack>`_. Indeed, the routines of this library were first developed by Carl de Boor in F77. External FORTRAN libraries are in the folder ``src``. The folder ``pppack/lib`` contains the shared objects of the numerical FORTRAN libraries. Use ``bin/wrap.sh [libname]`` with a ``libname`` from ``src/.`` to compile the requested so lib.

The documentation is produced by Sphinx, using non-standard ReST directives (as math) in the Python docstrings. Additional sphinx-extensions are used too, see ``docs/doc-requirements.txt``. Enter the folder ``docs`` and run ``make html`` to build the documentation, then open the page ``_build/html/index.html`` with any browser.

This package is hosted on the `PyPI <https://pypi.python.org/pypi/pppack/>`_, with the documentation published on `ReadTheDocs.org <http://pppack.readthedocs.io>`_.


Contact
=======

Please send bug reports, patches and other feedback to:

Daniele dot Tomatis at gmail dot com


Contributions
=============

02/05/2017 - Esteban Szames, implementation of a few tests from Carl De Boor's book.
