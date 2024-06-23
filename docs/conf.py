# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# import os
# import sys
# sys.path.insert(0, os.path.abspath("../../"))
# import asteca
# __version__ = asteca.__version__


# -- Project information -----------------------------------------------------
project = "ASteCA"
copyright = "2024, Gabriel I Perren"
# author = 'Gabriel I Perren'

# Read version from pyproject.toml
with open("../pyproject.toml", encoding="utf-8") as pyproject_toml:
    __version__ = (
        next(line for line in pyproject_toml if line.startswith("version"))
        .split("=")[1]
        .strip("'\"\n ")
    )

# The full version, including alpha/beta/rc tags
version = __version__
release = __version__

rst_epilog = """.. |ProjectVersion| replace:: {versionnum}""".format(
    versionnum=version,
)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = [
    "autodoc2",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "myst_nb",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_math_dollar",
    "sphinx.ext.mathjax",
    # "sphinx.ext.doctest"
]

autodoc2_packages = [
    {
        "path": "../asteca",
        "exclude_dirs": [
            "__pycache__",
            "modules",
        ],
    }
]
autodoc2_hidden_objects = ["dunder", "private", "inherited"]


################################################################
# This block removes the attributes from the apidocs .rst files.
# I could not find a simpler way to do this 26/05/24
# https://www.sphinx-doc.org/en/master/extdev/appapi.html#sphinx-core-events
def source_read_handler(app, docname, source):
    """'docname, source' not used but required"""
    path = "./apidocs/asteca/"
    for file in os.listdir(path):
        with open(path + file) as f:
            lines = f.readlines()
            idxs = []
            for i, line in enumerate(lines):
                if ".. py:attribute:" in line:
                    idxs += list(range(i, i + 6))
        for i in range(len(lines), 0, -1):
            if i in idxs:
                del lines[i]
        with open(path + file, "w") as f:
            for line in lines:
                f.write(line)


def setup(app):
    app.connect("source-read", source_read_handler)


################################################################


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

# templates_path = ['_templates']
# exclude_patterns = []


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
toc_object_entries_show_parents = "hide"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
