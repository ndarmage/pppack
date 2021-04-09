import os
import sys
import setuptools

# to use f2py extensions
from numpy.distutils.core import setup, Extension

here = os.path.abspath(os.path.dirname(__file__))

on_rtd = os.environ.get('READTHEDOCS') == 'True'

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

def get_value(ifile, *args):
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
# to use a consistent encoding
from codecs import open
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

if on_rtd:
    # workaround to avoid the issue of mocking lib.pppack import
    mainfile = os.path.join(here, 'pppack', 'pppack.py')
    pkgversion = get_value(mainfile, '__version__')[0]
    pkgname = get_value(mainfile, '__version__')[0]
else:
    from pppack.pppack import __version__ as pkgversion
    from pppack.pppack import __title__ as pkgname


sys.path.extend('config_fc --fcompiler=gnu95'.split())
sys.exit()

with open('docs/doc-requirements.txt','r') as f:
    requirements = f.readlines()
    requirements = [r.replace('\n','') for r in requirements]

if on_rtd:
    #requirements = [r for r in requirements if not 'numpy' in r]
    ext_mods = []

setup(
    name = pkgname,
    version = pkgversion + '.post1',

    description = 'A Python Piecewise Polynomial Package',
    long_description = long_description,

    author = 'Daniele Tomatis',
    author_email = 'daniele.tomatis@gmail.com',

    url = 'https://github.com/ndarmage/pppack',

    ext_modules = ext_mods,

    license = 'MIT',

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

    keywords = ['multivariate polynomial', 'splines', 'fortran'],

    packages = ['pppack'],

    install_requires = requirements,
    
    include_package_data = True,
    
    #scripts=['bin/wrap.sh'],
)
