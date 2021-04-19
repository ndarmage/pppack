@echo off

set libname=%~n0

set fsplit=..\..\..\..\bin\f90split.exe

mkdir temp
cd temp

if not exist %fsplit% (
    echo missing %fsplit% from %cd%
    exit /b 1
)
%fsplit% ..\%libname%.f90

for /F %%i in ('dir /b *.f90') do gfortran -c %%i
	
del *.f90

ar qc lib%libname%.a *.o
del *.o

REM mv lib${libname}.a ~/lib/$ARCH
move lib%libname%.a ..\.
cd ..
rmdir temp

REM echo "Library installed as ~/lib/$ARCH/lib${libname}.a"
echo "Library installed as ./lib%libname%.a"

REM test the library
gfortran -o %libname%_prb.exe %libname%_prb.f90 lib%libname%.a .\..\..\r8lib\f90\libr8lib.a
%libname%_prb.exe > %libname%_prb_output.txt
echo " --- TEST LIB ---"
echo "diff between ref. and current %libname%_prb_output.txt"
fc %libname%_prb_output_ref.txt %libname%_prb_output.txt
REM echo " *** without datetime at the beginning and at the end ***"
REM fc <(tail -n +1 ${libname}_prb_output.txt | head -n -1) <(tail -n +1 ${libname}_prb_output_tst.txt | head -n -1)
echo " --- NORMAL END ---"
