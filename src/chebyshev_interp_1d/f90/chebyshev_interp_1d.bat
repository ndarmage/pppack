@echo off

REM for /F %%i in ("%cd%") do set "libname=%%~ni"
set libname=%~n0

mkdir temp
cd temp
..\..\..\f90split\f90split.exe ..\%libname%.f90

for /F %%i in ('dir /b *.f90') do gfortran -c %%i
	
del *.f90

REM REM include dependencies of qr_solve
ar -x ..\..\..\qr_solve\f90\libqr_solve.a
ar -x ..\..\..\r8lib\f90\libr8lib.a
ar qc lib%libname%.a *.o
del *.o

move lib%libname%.a ..\.
cd ..
rmdir temp

echo "Library installed as lib%libname%.a"

REM test the library
echo "gfortran -o %libname%_prb.exe %libname%_prb.f90 .\lib%libname%.a"
gfortran -o %libname%_prb.exe %libname%_prb.f90 lib%libname%.a ..\..\r8lib\f90\libr8lib.a
%libname%_prb.exe > %libname%_prb_output.txt
echo " --- TEST LIB ---"
echo "diff between ref. and current %libname%_prb_output.txt"
fc %libname%_prb_output_ref.txt %libname%_prb_output.txt
REM echo " *** without datetime at the beginning and at the end ***"
REM FC <(tail -n +1 ${libname}_prb_output.txt | head -n -1) <(tail -n +1 ${libname}_prb_output_tst.txt | head -n -1)
echo " --- NORMAL END ---"
