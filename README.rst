PPPACK
======

:Author: Daniele Tomatis
:Date: 02/05/2017

Piece-wise Polynomial Package for multivariate interpolation of real-valued scalar functions.

This package offers representations of multi-dimensional functions using piece-wise polynomials. Python handles memory management, whereas external FORTRAN libraries perform all the numerics in order to ensure high computational performances with big amount of data.

The main FORTRAN library is an enhanced version of pppack.f90, whose original version by J. Burkardt is available on www.netlib.org/pppack. External FORTRAN libraries are available in the folder src. To get the shared object to import in a python module just enter the pppack/lib folder and run `.../bin/wrap.sh libname` where `libname` can be pppack, chebyshev_inter or divdif.

The documentation is produced by Sphinx, so additional directives (as math) are used in the docstrings. Please see also the additional sphinx.extension-s used before making the doc.


Contact
=======

Please send bug reports, patches and other feedback to:

Daniele . Tomatis ( at ) gmail . com


Contributions
=============

02/05/2017 - Esteban Szames, implementation of tests.
