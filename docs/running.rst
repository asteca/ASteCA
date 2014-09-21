Running
=======

To run the code simply position a cluster’s photometric data file in the
``/input`` folder, modify the ``params_input.dat`` file to both
inform the code of the distribution of columns in your photometric file
and set the values for each function within the code.
Photometric files can also be stored in ``/input`` inside a
sub-folder.

**Every** file inside the ``/input`` folder (inside or outside a
sub-folder) will be processed by **ASteCA**, with the exception of the
membership probabilities files that end with a ``.memb.dat`` extension
(see XXX).

Once file(s) are in place and the ``params_input.dat`` file correctly
modified, the code can be executed with the command:

.. code-block:: bash

    $ python asteca.py

The ``CLUSTER.DAT`` file located in the ``/input`` folder contains
a synthetic open cluster generated via the `MASSCLEAN`_ package with the
following parameter values:

::

    z=
    log(age)=
    E(B-V)=
    (m-M)o=

and serves as an example cluster to be analyzed with **ASteCA**.


Theoretical isochrones
----------------------

**ASteCA** needs at least one set of theoretical isochrones stored in a
``/asteca/isochrones`` folder to be able to apply the function that
estimates the clusters’ parameters. `Girardi isochrones`_ are currently
the default but any set can in practice be used (with minor changes made
to the code). Please `contact me <gabrielperren@gmail.com>`_ if you wish
to use a different set of theoretical isochrones.

The isochrones can be downloaded manually or the package `ezPadova-2`_
can be used to automatically fetch them from the site.


.. _MASSCLEAN: http://www.physics.uc.edu/~bogdan/massclean/
.. _Girardi isochrones: http://stev.oapd.inaf.it/cgi-bin/cmd
.. _ezPadova-2: https://github.com/Gabriel-p/ezpadova