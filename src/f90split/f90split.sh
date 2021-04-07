#! /bin/bash
#
gcc -c -Wall f90split.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gcc f90split.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm f90split.o
#
chmod u+x a.out
mv a.out ~/binc/f90split
#
echo "Normal end of execution."
