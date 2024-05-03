# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# import os
# import sys
# sys.path.insert(0, os.path.abspath("../../"))
# import asteca
# __version__ = asteca.__version__

# Read version from pyproject.toml
with open("../pyproject.toml", encoding="utf-8") as pyproject_toml:
    __version__ = (
        next(line for line in pyproject_toml if line.startswith("version"))
        .split("=")[1]
        .strip("'\"\n ")
    )


# -- Project information -----------------------------------------------------

project = "ASteCA"
copyright = "2024, Gabriel I Perren"
# author = 'Gabriel I Perren'

# The full version, including alpha/beta/rc tags
version = __version__
release = __version__

rst_epilog = """.. |ProjectVersion| replace:: {versionnum}""".format(
    versionnum = version,
)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # "sphinx.ext.autodoc",
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "myst_nb",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_math_dollar",
    "sphinx.ext.mathjax",
]

autoapi_dirs = ['../asteca']
# autoapi_options = ["imported-members", "show-inheritance", "inherited-members"]
# autoapi_ignore = ["_*"]
autoapi_add_toctree_entry = False
autoapi_keep_files = True
# autoapi_own_page_level = "attribute"


# Hide private files
def skip_submodules(app, what, name, obj, skip, options):
    if what in ("package", "function", "attribute"):
        skip = True
    if what == "method":  # Skip private methods
        if name.split('.')[-1].startswith("_"):
            skip = True
    # if skip is False:
    #     print(what, ",", name, ",", obj)
    #     # breakpoint()
    return skip
def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_submodules)


# https://www.sympy.org/sphinx-math-dollar/
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"]],
        "displayMath": [["\\[", "\\]"]],
    }
}

# https://myst-nb.readthedocs.io/en/v0.12.2/use/execute.html
# nb_execution_mode = "auto"
nb_execution_mode = "off"
# nb_execution_timeout = -1

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
]

# # Don't re-order methods alphabetically
# autodoc_member_order = "bysource"

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
html_theme = "sphinx_book_theme"
html_favicon = "_static/favicon.ico"
html_title = "[ ASteCA ]"
html_theme_options = {
    "repository_url": "https://github.com/asteca/ASteCA",
    "use_repository_button": True,
    "use_issues_button": True,
    "path_to_docs": "docs",
    "use_edit_page_button": True,
    "use_fullscreen_button": False,
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "repository_branch": "main",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "notebook_interface": "classic",
    },
}


# Hide parent class name in right sidebar TOC for methods
toc_object_entries_show_parents = 'hide'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
