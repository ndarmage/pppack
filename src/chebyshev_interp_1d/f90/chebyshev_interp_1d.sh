#!/bin/bash
export libname=$(basename $(dirname $PWD))
#export libname=$(basename $(dirname $(dirname $PWD)))

#
mkdir temp
cd temp
rm *
f90split ../${libname}.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
#
# include dependencies of qr_solve
ar -x ../../../qr_solve/f90/libqr_solve.a
ar -x ../../../r8lib/f90/libr8lib.a
ar qc lib${libname}.a *.o
rm *.o
#
#mv lib${libname}.a ~/lib/$ARCH
mv lib${libname}.a ../.
cd ..
rmdir temp
#
#echo "Library installed as ~/lib/$ARCH/lib${libname}.a"
echo "Library installed as ./lib${libname}.a"
#
# test the library
echo "gfortran -o ${libname}_prb.exe ${libname}_prb.f90 ./lib${libname}.a"
gfortran -o ${libname}_prb.exe ${libname}_prb.f90 ./lib${libname}.a ../../r8lib/f90/libr8lib.a
./${libname}_prb.exe > ${libname}_prb_output.txt
echo " --- TEST LIB ---"
echo "diff between ref. and current ${libname}_prb_output.txt"
diff ${libname}_prb_output_ref.txt ${libname}_prb_output.txt
#echo " *** without datetime at the beginning and at the end ***"
#diff <(tail -n +1 ${libname}_prb_output.txt | head -n -1) <(tail -n +1 ${libname}_prb_output_tst.txt | head -n -1)
echo " --- NORMAL END ---"
