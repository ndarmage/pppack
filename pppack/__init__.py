"""
.. include:: ../../docs/source/figures/pppack_asciiart.rst

The Fortran libraries in the lib folder must be recompiled and wrapped by F2Py
on every machine.

.. note::
   The folder bin contains scripts to get c-extension and
   fortran-wrapped libs.

.. warning::
    When using gcc compilers on WINOS, users may need to copy the following
    libraries in the folder lib in order to resolve missing dependencies:
    libgcc_s_seh-1.dll, libgfortran-5.dll, libquadmath-0.dll and
    libwinpthread-1.dll. The pyd or so extensions can be symlinked in
    pppack\lib.
"""
from pppack.pppack import *
