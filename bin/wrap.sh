#!/bin/bash

# this means the script will run successfully only if executed at pppack/lib
dir=../../src/
hme=$(pwd)

# concatenate all fortran sources and wrap the whole into a py module
if [ -z "$1" ]; then
    echo "ERROR, missing input libname"
    exit 1
fi
lib=$dir$1
if [ ! -d $lib ]; then
    echo "ERROR, invalid libname $lib"
    exit 1
fi

libopt=''
# apply special options according to input lib
if [[ $1 == *"chebyshev"* ]]; then
  libopt="only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d :  $dir/qr_solve/f90/qr_solve.f90 only: qr_solve :  $dir/r8lib/f90/r8lib.f90"
fi

# local compilation may be necessary before wrapping
echo "compile the lib $1 on the local machine"
cd $lib/f90
./$1.sh
cd $hme

# add directives before calling f2py
echo "produce the signature file"
f2py $lib/f90/$1.f90 -m $1 -h $lib/$1.pyf --overwrite-signature $libopt
#
echo "produce the so lib"
f2py -c --fcompiler=gnu95 -DF2PY_REPORT_ON_ARRAY_COPY --build-dir $lib/cc_$1/ $lib/$1.pyf $lib/f90/$1.f90 $libopt
