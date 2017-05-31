#!/bin/bash

# concatenate all fortran sources and wrap the whole into a py module
export lib=$(basename $PWD)

# add directives before calling f2py
f2py src/f90/$lib.f90 -m $lib -h $lib.pyf --overwrite-signature
f2py -c $lib.pyf src/f90/$lib.f90
