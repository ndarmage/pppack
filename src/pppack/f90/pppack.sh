#!/bin/bash
#export libname=$(basename $(dirname $(dirname $PWD)))
export libname=$(basename $(dirname $PWD))

fsplit=../../../../bin/f90split

mkdir temp
cd temp

if [ ! -f $fsplit ]; then
  echo "missing $fsplit from `pwd`"
  exit 1
fi
# the original file is pppack_orig.f90
$fsplit ../${libname}.f90

# first compile modules
f90mods=(colloc_data.f90 l2approx_data.f90 l2ir_data.f90 ppcolloc_data.f90)
for FILE in ${f90mods[@]};
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done

# then the remaining f90 files
for FILE in $(ls -1 *.f90 | grep -v _data);
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done

ar qc lib${libname}.a *.o
# rm *.f90 *.mod *.o

#mv lib${libname}.a ~/lib/$ARCH
mv lib${libname}.a ../.
cd ..
rm -R temp

#echo "Library installed as ~/lib/$ARCH/lib${libname}.a"
echo "Library installed as lib${libname}.a"

# test the library
gfortran -o ${libname}_prb.exe ${libname}_prb.f90 ./lib${libname}.a
./${libname}_prb.exe > ${libname}_prb_output.txt
echo " --- TEST LIB ---"
echo "diff between ref. and current ${libname}_prb_output.txt"
diff ${libname}_prb_output_ref.txt ${libname}_prb_output.txt
#echo " *** without datetime at the beginning and at the end ***"
#diff <(tail -n +1 ${libname}_prb_output.txt | head -n -1) <(tail -n +1 ${libname}_prb_output_tst.txt | head -n -1)
echo " --- NORMAL END ---"
