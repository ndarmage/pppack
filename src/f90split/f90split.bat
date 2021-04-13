@echo off

REM gcc compilation fails for missing link and unlink functions declaration
REM gcc -c -Wall f90split.c
REM if [ $? -ne 0 ]; then
  REM echo "Compile error."
  REM exit
REM fi

REM gcc f90split.o -o f90split.exe
REM if [ $? -ne 0 ]; then
  REM echo "Load error."
  REM exit
REM fi

REM rm f90split.o

gfortran -o f90split.exe f90split.f90
copy f90split.exe ..\..\bin\f90split.exe

echo "Normal end of execution."
