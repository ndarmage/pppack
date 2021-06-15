#! /bin/bash
#
export libname=$(basename $(dirname $PWD))
gfortran -c -Wall ${libname}.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ${libname}.o ~/lib/${libname}.o
#
echo "Normal end of execution."
