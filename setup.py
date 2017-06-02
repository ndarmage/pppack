from os import path

here = path.abspath(path.dirname(__file__))

import setuptools

# to use f2py extensions
from numpy.distutils.core import setup,  Extension

import setuptools

# list f2py extensions
ext1 = Extension(name = 'pppack.lib.pppack',
                 sources = ['src/pppack/f90/pppack.f90'],
                 define_macros = [('DF2PY_REPORT_ON_ARRAY_COPY',''),],
                 )
ext2 = Extension(name = 'pppack.lib.chebyshev_interp_1d',
                 sources = ['src/chebyshev_interp_1d/f90/chebyshev_interp_1d.f90',
                            'src/qr_solve/f90/qr_solve.f90',
                            'src/r8lib/f90/r8lib.f90'],
                 define_macros = [('DF2PY_REPORT_ON_ARRAY_COPY',''),],
                 )
ext3 = Extension(name = 'pppack.lib.divdif',
                 sources = ['src/divdif/divdif.pyf','src/divdif/f90/divdif.f90'],
                 define_macros = [('DF2PY_REPORT_ON_ARRAY_COPY',''),],
                 )
ext_mods = [ext1, ext2, ext3]

# Get the long description from the README file
# to use a consistent encoding
from codecs import open
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

from pppack.pppack import __version__ as pkgversion
from pppack.pppack import __title__ as pkgname

import sys, os
sys.path.extend('config_fc --fcompiler=gnu95'.split())

on_rtd = os.environ.get('READTHEDOCS') == 'True'

with open('docs/doc-requirements.txt','r') as f:
    requirements = f.readlines()

if on_rtd:
    requirements = [r for r in requirements if not 'numpy' in r]
    ext_mods = []

setup(
    name = pkgname,
    version = pkgversion+'.post1',
    
    description       = 'A Python Piecewise Polynomial Package',
    long_description = long_description,
    
    author            = 'Daniele Tomatis',
    author_email      = 'daniele.tomatis@gmail.com', # google groups ?
    
    url='https://github.com/ndarmage/pppack',
    
    ext_modules = ext_mods,
    
    license='MIT',
    
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    
    keywords=['multivariate polynomial','splines','fortran'],
    
    packages=['pppack'],

    install_requires=requirements,
    
    include_package_data=True,
    
    #scripts=['bin/wrap.sh'],
)






