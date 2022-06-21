PPPACK
======

:Author: Daniele Tomatis
:Date: 02/05/2017

Introduction
------------

Piece-wise Polynomial Package for multivariate interpolation of real-valued scalar functions.

This package offers representations of multi-dimensional functions using piece-wise polynomials. Python handles memory management, whereas external FORTRAN libraries perform all the numerics in order to ensure high computational performances with a large amount of data.

The main FORTRAN library is an improved version of the original Fortran77 routines first developed by Carl de Boor and still available on `netlib <http://www.netlib.org/pppack>`_. They are also contained in this package. Some of the routines were later translated to FORTRAN90 by John Burkardt. Finally, Daniele Tomatis modified the sources to fully comply with standard FORTRAN90; code upgrade with the latest FORTRAN 2018 norm is expected in future versions. Other numerical libraries addressing interpolation problems and approximation theory are available with this package.

All external FORTRAN libraries are in the folder ``src``. The folder ``pppack/lib`` contains the Python extensions of the numerical FORTRAN libraries (as shared objects). Use ``tools/wrap.sh [libname]`` where ``libname`` is a numerical library available in ``src`` to obtain the corresponding Python extension.

The documentation is produced by Sphinx using the math ReST directive in Python docstrings. Additional sphinx-extensions are used, see the list of required packages in the file ``setup.py``. Enter the folder ``docs`` and execute ``make html`` to build the documentation in HTML format, then open the page ``_build/html/index.html`` with the browser.

This package is hosted on `PyPI <https://pypi.python.org/pypi/pppack/>`_, with the documentation published on `ReadTheDocs.org <http://pppack.readthedocs.io>`_.

Installation
------------

To retrieve the development version:

.. code-block:: console

    $ git clone https://github.com/ndarmage/pppack.git

To install the package:

.. code-block:: console

    $ pip install .

or get it directly from PyPI,

.. code-block:: console

    $ pip install pppack

To compile only the Fortran extensions:

.. code-block:: console

    $ python setup.py build_src build_ext --inplace

Users can select different Fortran compiler by adding ``config --fcompiler=gnu95`` before ``build_ext`` (replace ``gnu95`` with your compiler after checking the available choices with ``f2py``).


Contact
-------

Please send bug reports, patches and other feedback to `Daniele Tomatis <mailto:daniele.tomatis@gmail.com>`_.


Contributions
-------------

The main contributions featuring the released versions are displayed in the following table.

.. csv-table:: List of major contributions in the released versions.
   :header: "Version", "Date", "Author(s)", "Description"
   :widths: 10, 10, 30, 50
   :delim: &

   v1.2.0 & 21/06/2022 & Dinh Nguyen, Daniele Tomatis & upgrade of the package using only distutils on Python v3.9
   v1.1.0.post1 & 09/04/2021 & Daniele Tomatis & port pppack to Win10; enforce PEP8 norm in py modules.
   v1.0.0 & 02/05/2017 & Esteban Szames & implementation of a few tests from Carl De Boor's book.
