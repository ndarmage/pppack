Requirements
============

* Numpy + f2py (>= 1.11.0)
* itertools (for product yielding cartesian products of elements
  belonging to multiple lists)
* pppack.f90 (see the lib folder for more information)

Former theory background is necessary to master the functionalities 
of this package; previous reading of C. De Boor's book `A practical 
guide to Splines` is recommended 
:cite:`deBoor1978prguide,deBoor2001prguide`.


Acknowledgments
===============

This package follows closely the notation and the theory from 
C. De Boor's book: A practical guide to Splines 
:cite:`deBoor1978prguide,deBoor2001prguide`. The revised edition of 
2001 is recommended.

J. Burkardt provided a first translation of De Boor's F77 routines to F90.
The version of the library used in this project was refactored by D. 
Tomatis in order to increase readability (e.g. GOTO removal) and to
comply with more standard FORTRAN code. E. Szames is responsible of the
testing environment.

I am greatly thankful to `Pearu Peterson`_, the main developer of f2py_,
for his precious suggestions about the directives and the tests with
wrapped FORTRAN modules.

.. _`Pearu Peterson`: http://cens.ioc.ee/~pearu/

.. _f2py: https://docs.scipy.org/doc/numpy-dev/f2py/


Bibliography
============

.. bibliography:: pppack.bib
   :style: unsrt
