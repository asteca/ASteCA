# -- Project information -----------------------------------------------------
project = "ASteCA"
copyright = "2024, Gabriel I Perren"
author = "Gabriel I Perren"

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

rst_epilog = """.. |ProjectVersion| replace:: **{versionnum}**""".format(
    versionnum=version,
)


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = [
    "sphinx.ext.githubpages",
    "myst_nb",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx.ext.mathjax",
    "sphinx_inline_tabs",
    "sphinx_automodapi.automodapi",
]

#  automodapi
numpydoc_show_class_members = False
automodapi_toctreedirnm = "build/api"

# https://myst-nb.readthedocs.io/en/v0.12.2/use/execute.html
# nb_execution_mode = "auto"
nb_execution_mode = "off"
# nb_execution_excludepatterns = ["*membership"]
# nb_execution_timeout = -1


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
        "colab_url": "https://colab.research.google.com",
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
html_css_files = ["custom.css"]
html_js_files = ["umami.js"]
