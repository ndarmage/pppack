#!/bin/bash

# concatenate all fortran sources and wrap the whole into a py module
export lib=$(basename $PWD)

# add directives before calling f2py
f2py f90/$lib.f90 -m $lib -h $lib.pyf --overwrite-signature
f2py -c --fcompiler=gnu95 $lib.pyf f90/$lib.f90
