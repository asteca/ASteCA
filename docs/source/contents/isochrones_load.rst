.. _isochronesload:

Isochrones
==========


``ASteCA`` requires a set of theoretical isochrones to be able to estimate the
fundamental parameters for the clusters. The supported services are:
`PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd_3.7>`_,
`MIST <https://waps.cfa.harvard.edu/MIST/>`_, and
`BASTI <http://basti-iac.oa-abruzzo.inaf.it/isocs.html>`_.


Each service produces the isochrone files in a different way. A single file
produced with these services will contain:

* PARSEC : multiple metallicities and multiple ages
* MIST   : single metallicity and multiple ages
* BASTI  : single metallicity and single age

After downloading the isochrone files from whichever service chosen, these need
to be stored in a folder named after the service used to produce them. It is important
to use the correct name of the service used to generate the file(s), so that the
:class:`isochrones` class can then load the files properly.

The block below shows how the files should be stored for each service. A single
``>`` character points to a folder, while two ``>>`` characters point to a file.

.. code-block:: console

    > PARSEC/
      |---> phot_syst/
            |-------->> mets_ages.dat

    > MIST/
      |---> phot_syst/
            |-------->> met_age_1.iso.cmd
            |-------->> met_age_2.iso.cmd
            |-------->> ...

    > BASTI/
      |---> phot_syst/
            |--------> met_1/
                       |---->> age_1.isc_xxxx
                       |---->> age_2.isc_xxxx
                       |---->> ...
            |--------> met_2/
                       |---->> age_1.isc_xxxx
                       |---->> ...

In this example ``phot_syst`` is the name of the photometric system that is
being employed (not required but suggested), ``mets_ages`` is the single file
produced by the ``PARSEC`` service, ``met_age_X.iso.cmd`` is the name of a file
produced by the ``MIST`` system, and ``age_X.isc_xxxx`` is the names of a file produced
by the ``BASTI`` system, stored in a sub-folder named after its metallicity
``met_x``.


Please `contact me <gabrielperren@gmail.com>`_ if you have any issues with the loading
process of the theoretical isochrones.


