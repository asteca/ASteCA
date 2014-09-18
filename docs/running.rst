Running
-------

To run the code simply drop a cluster’s photometric data file in the
``/asteca/input`` folder, modify the ``params_input.dat`` file to both
inform the code of the distribution of columns in your photometric file
and set the values for each function within the code.

Photometric files can also be stored in ``/asteca/input`` inside a
sub-folder.

**Every** file inside the ``/asteca/input`` folder (inside or outside a
sub-folder) will be processed by ``ASteCA``, with the exception of the
membership probabilities files that end with a *memb.dat* extension.

Once file(s) are in place and the ``params_input.dat`` file correctly
modified, the code can be executed with the command:

::

    python asteca.py

The *CLUSTER.DAT* file located in the ``/asteca/input`` folder contains
a synthetic open cluster generated via the `MASSCLEAN`_ package with the
following parameter values:

::

    z=
    log(age)=
    E(B-V)=
    (m-M)o=

and serves as an example cluster to be analyzed with ``ASteCA``.

.. _MASSCLEAN: http://www.physics.uc.edu/~bogdan/massclean/