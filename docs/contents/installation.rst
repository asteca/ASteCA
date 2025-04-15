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
.. _venv: https://docs.python.org/3/library/venv.html