import os
import sys

# to use f2py extensions
from numpy.distutils.core import setup, Extension

on_rtd = os.environ.get('READTHEDOCS') == 'True'

# list f2py extensions
extensions = [] if on_rtd else [
    Extension(name='pppack.lib.pppack',
              sources=['src/pppack/f90/pppack.f90'],
              define_macros=[('DF2PY_REPORT_ON_ARRAY_COPY', None)]
              # f2py_options=["--fcompiler=gnu95"],
              ),
    Extension(name='pppack.lib.chebyshev_interp_1d',
              sources=['src/chebyshev_interp_1d/f90/chebyshev_interp_1d.f90',
                       'src/qr_solve/f90/qr_solve.f90',
                       'src/r8lib/f90/r8lib.f90'],
              define_macros=[('DF2PY_REPORT_ON_ARRAY_COPY', None)]
              # f2py_options=["--fcompiler=gnu95"]
              ),
    Extension(name='pppack.lib.divdif',
              sources=['src/divdif/f90/divdif.f90'],
              define_macros=[('DF2PY_REPORT_ON_ARRAY_COPY', None)]
              # f2py_options=["--fcompiler=gnu95"]
              ),
    Extension(name='pppack.lib.hermite_cubic',
              sources=['src/hermite_cubic/f90/hermite_cubic.f90'],
              define_macros=[('DF2PY_REPORT_ON_ARRAY_COPY', None)]
              # f2py_options=["--fcompiler=gnu95"]
              ),
]  # no extension build if on RTD servers


def get_values(ifile, *args):
    "extract value of args keys in module ifile"
    o = []
    with open(ifile, 'r') as fid:
        lines = fid.readlines()

    for v in args:
        if isinstance(v, str):
            for line in lines:
                if v in line:
                    val = line.split(v)[1]
                    if '=' in val:
                        val = val.split('=')[1].strip()
                        if '#' in val:
                            val = val.split('#')[0].strip()
                    else:
                        raise RuntimeError('missing assignation for ' + v)
                    o.append(val.replace('"', '').replace("'",''))
        else:
            raise TypeError('input var ' + str(v) + ' is not of type str')

    return o

# Get the long description from the README file
with open("README.rst", 'r', encoding="utf-8") as f:
    long_description = f.read()

# workaround to avoid the issue of mocking lib.pppack import
mainfile = os.path.join('pppack', 'pppack.py')
pkgname, pkgversion = get_values(mainfile, '__title__', '__version__')
# the following is commented because it fails when build_ext is used
# from pppack.pppack import __version__ as pkgversion
# from pppack.pppack import __title__ as pkgname


with open(os.path.join("docs", "RTD-doc-requirements.txt"), 'r',
          encoding="utf-8") as f:
    required_packages = f.readlines()


setup(
    name=pkgname,
    version=pkgversion,
    packages=[pkgname],
    description='A Python Piecewise Polynomial Package',
    long_description=long_description,

    author='Daniele Tomatis',
    author_email='daniele.tomatis@gmail.com',

    url='https://github.com/ndarmage/pppack',

    ext_modules=extensions,

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],

    keywords=['approximation theory', 'multivariate polynomial', 'splines',
              'fortran'],

    install_requires=required_packages,
)

# use MANIFEST.in instead of include_package_data
