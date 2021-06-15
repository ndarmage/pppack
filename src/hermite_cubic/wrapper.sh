#!/bin/bash

# concatenate all fortran sources and wrap the whole into a py module
export lib=$(basename $PWD)

# add directives before calling f2py
python3 -m numpy.f2py f90/$lib.f90 -m $lib -h $lib.pyf --overwrite-signature
python3 -m numpy.f2py -c --fcompiler=gnu95 -DF2PY_REPORT_ON_ARRAY_COPY $lib.pyf f90/$lib.f90
