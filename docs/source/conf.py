# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
# print(sys.path)
sys.path.insert(0, os.path.abspath('../../'))

# from pkg_resources import DistributionNotFound, get_distribution
# # try:
# import asteca
# __version__ = asteca.__version__
# # __version__ = get_distribution("asteca").version
# # except DistributionNotFound:
# #     __version__ = "unknown version"
# print(__version__)

__version__ = '0.5.1'


# -- Project information -----------------------------------------------------

project = 'ASteCA'
copyright = '2024, Gabriel I Perren'
author = 'Gabriel I Perren'

# The full version, including alpha/beta/rc tags
version = __version__
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    'sphinx.ext.githubpages',
    "myst_nb",
    "IPython.sphinxext.ipython_console_highlighting",
]

# https://myst-nb.readthedocs.io/en/v0.12.2/use/execute.html
nb_execution_mode = "cache"

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_favicon = "_static/favicon.ico"
html_title = "[ ASteCA ]"
html_theme_options = {
    "repository_url": "https://github.com/asteca/ASteCA",
    "use_repository_button": True,
    "use_fullscreen_button": False,
    "show_toc_level": 1,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
