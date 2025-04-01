.. _cluster_module:

Cluster
#######

The :py:class:`asteca.cluster.Cluster` class in **ASteCA** allows generating
an object with data for either a cluster or an observed field (i.e.: a data file
that contains not only the cluster members but also the surrounding field stars)

.. code-block:: python

    import asteca

    # Generate a `cluster` object
    my_field = asteca.cluster(
        ra=df['RA_ICRS'],
        dec=df['DE_ICRS'],
        magnitude=df["Gmag"],
        e_mag=df["e_Gmag"],
        color=df["BP-RP"],
        e_color=df['e_BP-RP'],
        pmra=df["pmRA"],
        e_pmra=df["e_pmRA"],
        pmde=df["pmDE"],
        e_pmde=df["e_pmDE"],
        plx=df["Plx"],
        e_plx=df["e_Plx"],
        verbose=2            # Print info to screen
    )

    Instantiating cluster...
    Columns read   : RA, DEC, Magnitude, e_mag, Color, e_color, Plx, pmRA, pmDE
    N_stars        : 8683
    N_clust_min    : 25
    N_clust_max    : 5000
    Cluster object generated

This object, called ``my_field`` here, can be processed to extract a given cluster's
properties such as: its center coordinates and radius (:ref:`structure`),
the estimated number of members (:ref:`nmembers`), as well as their membership
probabilities (:ref:`membership_module`).

The following sections explain these processes in more detail.


.. _structure:

Structure analysis
******************

Structural analysis includes a cluster's center coordinates estimation as well as the
estimation of its radius.


Center estimation
=================

The simplest structural analysis (and usually the first one to be required) is the
cluster's center estimation. **ASteCA** provides the
:py:meth:`get_center() <asteca.cluster.Cluster.get_center>` method to perform this
estimation.

There are currently two algorithms available for center estimation. The default method
is called ``knn_5d`` and it requires that your :py:class:`cluster` object contains
the five dimensions of data: ``(ra, dec, pmra, pmde, plx)`` . The other algorithm is
called ``kde_2d`` and it requires either ``(ra, dec)`` or ``(pmra, pmde)`` 2D data,
depending on the center values you want to estimate.





.. code-block:: python

    my_cluster.get_center()

    >> Center coordinates found:
    >> radec_c        : (6.3049, 61.3218)
    >> pms_c          : (-2.811, -1.070)
    >> plx_c          : 0.288

The ``radec_c, pms_c, plx_c`` values containing the center coordinates will be stored in
your :py:class:`my_cluster` object as attributes and can be accessed to use later
on:




Radius estimation
=================

A method to estimate the cluster's radius  will be added in future versions.
In the meantime you can manually add the attribute with:

.. code-block:: python

    my_cluster.radius = 0.1

where the value is always in units of degrees.


.. _nmembers:

Number of members
=================

Estimating the number of members for a given cluster is a crucial step for the
membership analysis. Currently **ASteCA** integrates two methods to perform this
estimation, as shown in
:py:meth:`get_nmembers() <asteca.cluster.Cluster.get_nmembers>`.

If the estimated number is not a proper representation of the believed number of members
for this cluster, the user can easily input this value manually with:

.. code-block:: python

    my_cluster.N_cluster = 300



.. include:: cluster.ipynb
   :parser: myst_nb.docutils_

