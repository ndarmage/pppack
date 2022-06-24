@Echo Off

for /F %%i in ("%cd%") do set "lib=%%~ni"

set PY=python -m numpy.f2py
set OPTS=--verbose --fcompiler=gnu95

set LDFLAGS=-static

REM add directives before calling f2py
%PY% -m %lib% -h %lib%.pyf --overwrite-signature ^
    f90\%lib%.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : ^
    ..\qr_solve\f90\qr_solve.f90 only: qr_solve :
%PY% %OPTS% -c %lib%.pyf ^
    f90/%lib%.f90 only: chebyshev_coef_1d chebyshev_interp_1d chebyshev_value_1d : ^
    ../qr_solve/f90/qr_solve.f90 only: qr_solve : ^
    ../r8lib/f90/r8lib.f90