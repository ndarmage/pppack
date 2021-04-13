#!/bin/bash

# concatenate all fortran sources and wrap the whole into a py module
export lib=$(basename $PWD)

# add directives before calling f2py
python3 -m numpy.f2py -m $lib -h $lib.pyf --overwrite-signature \
   f90/$lib.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : \
   ../qr_solve/f90/qr_solve.f90 only: qr_solve :
python3 -m numpy.f2py -c --fcompiler=gnu95 $lib.pyf \
   f90/$lib.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : \
   ../qr_solve/f90/qr_solve.f90 only: qr_solve : \
   ../r8lib/f90/r8lib.f90
