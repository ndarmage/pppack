@echo off

for /F %%i in ("%cd%") do set "lib=%%~ni"

set PY=python -m numpy.f2py
set OPTS=--verbose --fcompiler=gnu95 --compiler=mingw32

REM add directives before calling f2py
%PY% -m %lib% -h %lib%.pyf %OPTS% --overwrite-signature f90\%lib%.f90
%PY% %OPTS% -DF2PY_REPORT_ON_ARRAY_COPY -c %lib%.pyf f90/%lib%.f90
