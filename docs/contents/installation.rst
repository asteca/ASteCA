.. _installation:

Installation
############

**ASteCA** is a pure Python package. It can be installed locally or in a
`Google Colaboratory <https://colab.google/>`_ notebook. For a local installation
it is recommended to install it in its own environment, so that its dependencies do
not conflict with those of other packages. Any of the tools below can be used for
this; `uv`_ is the recommended one.

All the commands below install ``asteca[all]``, i.e.: **ASteCA** together with its
only optional dependency, `Numba`_.

.. note::

   `Numba`_ is not required but it is **highly** recommended. If it is installed, **ASteCA**
   uses it to compile its most demanding functions (e.g.: the generation of synthetic
   clusters and the likelihood evaluation), which makes the package considerably
   faster. If it is not installed, **ASteCA** falls back to equivalent pure NumPy
   functions. See :ref:`install_no_numba`.


uv (recommended)
================

`uv`_ is a fast Python package and project manager, that also installs Python
itself if needed. After `installing uv <https://docs.astral.sh/uv/getting-started/installation/>`_,
create a project and add **ASteCA** to it:

.. code-block:: shell

   uv init my_project
   cd my_project
   uv add "asteca[all]"

Scripts are then run inside the project's environment with ``uv run``:

.. code-block:: shell

   uv run python my_script.py

To install **ASteCA** in an existing environment instead, use
``uv pip install "asteca[all]"``.


venv
====

`venv`_ is included with Python. Create an environment in the ``.venv`` folder,
activate it, and install **ASteCA** with `pip`_:

.. tab:: Unix/macOS

   .. code-block:: shell

      python -m venv .venv
      source .venv/bin/activate
      python -m pip install "asteca[all]"

.. tab:: Windows

   .. code-block:: shell

      py -m venv .venv
      .venv\Scripts\activate
      py -m pip install "asteca[all]"

The environment must be activated (second line) in every new terminal session
before using **ASteCA**.


conda
=====

**ASteCA** is not distributed through `conda`_ channels, but it can be installed
with `pip`_ inside a conda environment. Create the environment with Python and
`Numba`_ from ``conda-forge``, activate it, and install **ASteCA**:

.. code-block:: shell

   conda create -n asteca -c conda-forge python=3.13 numba
   conda activate asteca
   python -m pip install "asteca[all]"

The environment must be activated (second line) in every new terminal session
before using **ASteCA**.


Google Colab
============

Run in a notebook's code cell:

.. code-block:: bash

   !pip install "asteca[all]"

If Colab asks to restart the session after the installation (this happens when some
of the packages already loaded by Colab, e.g. NumPy, are upgraded), restart it and
continue with the next cell. The installation must be repeated in every new Colab
session.


.. _install_no_numba:

Installing without Numba
========================

If `Numba`_ can not be installed (for example, if it does not yet support your
Python version or platform), install ``asteca`` instead of ``asteca[all]`` in any
of the commands above, e.g.:

.. code-block:: shell

   uv add asteca

`Numba`_ can be added later by installing ``asteca[all]``.




.. _pip: https://pip.pypa.io/en/stable/
.. _conda: https://docs.conda.io/
.. _uv: https://docs.astral.sh/uv/
.. _Numba: https://numba.pydata.org/
.. _venv: https://docs.python.org/3/library/venv.html
