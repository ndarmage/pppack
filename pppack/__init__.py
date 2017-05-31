"""
::
  
      _____    _____    _____               _____   _  __
     |  __ \  |  __ \  |  __ \     /\      / ____| | |/ /
     | |__) | | |__) | | |__) |   /  \    | |      | ' / 
     |  ___/  |  ___/  |  ___/   / /\ \   | |      |  <  
     | |      | |      | |      / ____ \  | |____  | . \ 
     |_|      |_|      |_|     /_/    \_\  \_____| |_|\_\
  
  
  
The Piece-wise Polynomial Package
=================================

Overview
--------

This python package offers Piece-wise Polynomial (PP) and B-spline
representations of multi-dimensional scalar real functions.

The main module pppack.py contains the definition of different
classes and methods used to build function spaces, functions and
tensors (intended here as multi-dimensional arrays). Multi-
dimensionality is achieved by Cartesian products, which are 
implemented by overloading the multiplication methods.

Memory management is handled by python, whereas all the numerics is
delegated to external FORTRAN libraries to achieve high computational
performances. These libraries are C-wrapped by f2py_ 
:cite:`peterson2009f2py`:

 1. pppack.f90 (main library)

    This F90 library is produced from the original Carl de Boor's
    F77 library_ :cite:`deBoor1978prguide`. The current version 
    used here was modified by D. Tomatis to comply with the latest
    F90 standards. Original files_ comes from J. Burkardt.

 2. newton_interp_1d.f90 (here appended at the end of pppack.f90)
    
    Interpolation of polynomials in Newton form by divided diffe- 
    rences. The code has been upgraded with routines for multi-
    dimensional data.

 3. chebyshev_interp_1d.f90
    
    Interpolation by Chebyshev polynomials.

The implementation of multi-dimensional approximations follows De 
Boor's array indexing as in :cite:`deBoor1979TOMS`.

.. note:: The use of the library DIVDIF to work with polynomials in
          Newton form is on stand-by status at present for the 
          commented lines in data_to_dif which prevents osculatory
          fitting. However, please note that this library contains
          many more functionalities than newton_interp_1d.

.. _f2py: https://sysbio.ioc.ee/projects/f2py2e/

.. _library: http://pages.cs.wisc.edu/~deboor/pgs/

.. _files: http://people.sc.fsu.edu/~jburkardt/f_src/pppack/pppack.f90
"""
# __all__ is here commented to allow importing the content of the module 
# pppack.py without the need of pppack.item_name. This may change in future
# according to packaging constraints.
#__all__ = ['pppack']
from pppack import *
