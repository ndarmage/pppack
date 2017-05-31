#!/bin/bash

# concatenate all fortran sources and wrap the whole into a py module
export lib=$(basename $PWD)

# add directives before calling f2py
f2py -m $lib -h $lib.pyf --overwrite-signature \
   src/f90/$lib.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : \
   ../qr_solve/src/f90/qr_solve.f90 only: qr_solve :
f2py -c --fcompiler=gfortran $lib.pyf \
   src/f90/$lib.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : \
   ../qr_solve/src/f90/qr_solve.f90 only: qr_solve : \
   ../r8lib/src/f90/r8lib.f90
