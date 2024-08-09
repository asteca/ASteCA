.. _installation:

Installation
############

**ASteCA** can be installed locally or in a
`Google Colaboratory <https://colab.google/>`_ notebook. To install in Google
Colaboratory use the following command in a new notebook:

.. code-block:: bash

  !pip install asteca

To check that the installation was successful, import **ASteCA** and print the
installed version number with:

.. code-block:: bash

  import asteca
  asteca.__version__

If you want to install the package locally, I recommend using the `conda`_ package and environment manager to install **ASteCA** in an isolated Python environment.
To install with ``conda``, follow these steps:

1. Go to https://conda.io/miniconda.html and download the appropriate version
   for your system. I will assume in what follows that you are running a 64 bit Linux
   system.
2. Open a terminal instance where this file was downloaded, and install with the
   command:

   .. code-block:: bash

     $ bash Miniconda3-latest-Linux-x86_64.sh

   Select yes when asked: *Do you whish the installer to prepend the Miniconda3
   install location to PATH in your ~/path?*
3. Restart your terminal for the changes to take effect.
4. Create a virtual environment for **ASteCA** with the command:

   .. code-block:: bash

     $ conda create --name asteca python=3.10

5. Activate the environment:

   .. code-block:: bash

     $ conda activate asteca


   You can tell that the environment is activated because its name is now
   shown in the terminal before the ``$`` symbol as:

   .. code-block:: bash

     (asteca) $

   You need to activate this environment each time *before* attempting to
   run **ASteCA**.

6. Finally, install **ASteCA** with:

   .. code-block:: bash

     (asteca) $ pip install asteca


If everything went well, you can now access a ``Python`` shell within your ``conda``
environment and import **ASteCA**:

.. code-block:: bash

 (asteca) $ python
 Python 3.10.14 (main, Mar 21 2024, 16:24:04) [GCC 11.2.0] on linux
 Type "help", "copyright", "credits" or "license" for more information.
 >>> import asteca
 >>> asteca.__version__




.. _conda: https://conda.io/docs/index.html
