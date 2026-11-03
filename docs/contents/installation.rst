.. _installation:

Installation
############

**ASteCA** can be installed locally or in a
`Google Colaboratory <https://colab.google/>`_ notebook. To install **ASteCA** locally
run the following `pip`_ command in a terminal session:

.. tab:: Unix/macOS

   .. code-block:: shell

      python -m pip install asteca

.. tab:: Windows

   .. code-block:: shell

      py -m pip install asteca

It is recommended to set up a local environment first. There are several tools that
handle this: `conda`_, `uv`_, or Python's own `venv`_.

.. note::

   `Numba`_ is not required but it is highly recommended. If it is installed, **ASteCA**
   uses it to compile its most demanding functions (e.g.: the generation of synthetic
   clusters and the likelihood evaluation), which makes the package considerably
   faster. If it is not installed, **ASteCA** falls back to equivalent pure NumPy
   functions. To install it, run:

   .. code-block:: shell

      python -m pip install numba

To install **ASteCA** in Google Colaboratory, run in a notebook's code cell:

.. code-block:: bash

  !pip install asteca

To verify that the installation was successful, import **ASteCA** and print the
installed version number with:

.. code-block:: bash

  import asteca
  print(asteca.__version__)


.. _pip: https://pip.pypa.io/en/stable/
.. _conda: https://conda.io/docs/index.html
.. _uv: https://docs.astral.sh/uv/
.. _Numba: https://numba.pydata.org/
.. _venv: https://docs.python.org/3/library/venv.html