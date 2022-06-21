# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# PPPACK documentation build configuration file, first created
# on Tue May 2 17:15:50 2017.

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

def extract_value(txt, key):
    string = [line for line in txt if key in line][0]
    return string.split()[-1].strip('"')

with open(os.path.abspath('../../pppack/pppack.py'), 'r') as f:
    lines = f.readlines()
author = extract_value(lines, '__author__')
project = extract_value(lines, '__title__').upper()
release = extract_value(lines, '__version__')
copyright = '2022, ' + author


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinxcontrib.bibtex',
    'sphinx.ext.mathjax',
    'sphinxfortran.fortran_domain',
    'sphinxfortran.fortran_autodoc',
    'sphinx.ext.intersphinx'
]

# set the filepath to the fortran 90 extension libs
fortran_src = [
    '../src/chebyshev_interp_1d/f90/chebyshev_interp_1d.f90',
    '../src/pppack/f90/pppack.f90',
    '../src/divdif/f90/divdif.f90',
]

# link to external documentation of other python modules
intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    # 'f2py' : ('https://numpy.org/doc/stable/f2py/', None)
}

# change the member order used by autodoc
autodoc_member_order = 'bysource'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# main bibtex file
bibtex_bibfiles = ['pppack.bib']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'rainbow_dash'  # 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'  # possible alternative sphinx_rtd_theme
html_css_files = ['morestyle.css']
html_js_files = [('custom.js', {'defer': 'defer'})]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'description': 'Piecewise Polynomials Package',
    'github_user': 'ndarmage',
    'github_repo': 'https://github.com/ndarmage/pppack',
    'github_button': True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
   '**': ['about.html', 'navigation.html', 'searchbox.html','relations.html'],
}


# -- Options for LaTeX output ---------------------------------------------

# latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',

# Latex figure (float) alignment
#'figure_align': 'htbp',
# }

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'PPPACK.tex', u'PPPACK Documentation',
     u'Daniele Tomatis', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'pppack', u'PPPACK Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    ('index', 'PPPACK', u'PPPACK Documentation',
     author, 'PPPACK', 'One line description of project.',
     'Miscellaneous'),
]

# If true, do not generate a @detailmenu in the "Top" node's menu.
texinfo_no_detailmenu = False


# -- READTHEDOCS ----------------------------------------------------------

# Ignore some modules during documentation building on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS') == 'True'

if on_rtd:

    from unittest import mock

    autodoc_mock_imports = [
        'lib.pppack',
        'lib.chebyshev_interp_1d',
        'lib.divdif',
        # 'numpy',
    ]

    for mod_name in autodoc_mock_imports:
        sys.modules[mod_name] = mock.Mock()
