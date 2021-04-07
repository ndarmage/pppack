@echo off

set libname=%~n0

mkdir temp
cd temp
del *
REM ~/bin/$ARCH/f90split ../${libname}.f90
REM the original file is pppack_orig.f90
..\..\..\f90split\f90split.exe ..\%libname%.f90

REM first compile modules
set f90mods[0]=colloc_data.f90
set f90mods[1]=l2approx_data.f90
set f90mods[2]=l2ir_data.f90
set f90mods[3]=ppcolloc_data.f90

set "x=0"

:SymLoop
if defined f90mods[%x%] (
    call gfortran -c %%f90mods[%x%]%%
    set /a "x+=1"
    GOTO :SymLoop
)

REM then the remaining f90 files
for /F %%i in ('dir /b *.f90') do gfortran -c %%i

del *.f90 *.mod

ar qc lib%libname%.a *.o
del *.o

REM move lib${libname}.a ~/lib/$ARCH
move lib%libname%.a ..\.
cd ..
rmdir temp

REM #echo "Library installed as ~/lib/$ARCH/lib${libname}.a"
echo "Library installed as lib%libname%.a"

REM test the library
gfortran -o %libname%_prb.exe %libname%_prb.f90 lib%libname%.a
%libname%_prb.exe > %libname%_prb_output.txt
echo " --- TEST LIB ---"
echo "diff between ref. and current %libname%_prb_output.txt"
fc %libname%_prb_output_ref.txt %libname%_prb_output.txt
REM #echo " *** without datetime at the beginning and at the end ***"
REM #diff <(tail -n +1 ${libname}_prb_output.txt | head -n -1) <(tail -n +1 ${libname}_prb_output_tst.txt | head -n -1)
echo " --- NORMAL END ---"
