#! /bin/bash
#
# the c file seems bugged
#
#gcc -c -Wall f90split.c
gfortran -c -Wall f90split.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
# gcc f90split.o
gfortran f90split.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm f90split.o
#
chmod u+x a.out
mv a.out ../../bin/f90split
#
echo "Normal end of execution."
