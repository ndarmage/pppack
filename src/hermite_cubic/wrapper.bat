@Echo Off

for /F %%i in ("%cd%") do set "lib=%%~ni"

set PY=python -m numpy.f2py
set OPTS=--verbose --debug-capi --fcompiler=gnu95

REM set LDFLAGS="-lpthread"
REM "--disable-threads"
REM set CFLAGS="--disable-threads"
REM set FFLAGS="--disable-threads"
REM set NPY_DISTUTILS_APPEND_FLAGS=1

REM add directives before calling f2py
%PY% -m %lib% -h %lib%.pyf --overwrite-signature f90\%lib%.f90
%PY% %OPTS% -c %lib%.pyf f90\%lib%.f90