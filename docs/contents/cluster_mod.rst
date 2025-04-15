.. _cluster_module:

Cluster module
##############

The :py:class:`asteca.Cluster` class in **ASteCA** generates an object with
data for either an observed cluster (i.e.: a data file that contains only true
cluster members) or an observed field (i.e.: a data file that contains not just the
cluster members but also the surrounding field stars).

Given a data file containing either a cluster or field data, the user can generate
a :py:class:`asteca.Cluster` object as follows:


.. code-block:: python

    import pandas as pd
    import asteca

    # Load the data file as a pandas DataFrame
    df = pd.read_csv("my_data_file.csv")

    # Create a Cluster object using the columns in the loaded data file
    my_field = asteca.Cluster(
        ra=df["RA_ICRS"],
        dec=df["DE_ICRS"],
        pmra=df["pmRA"],
        pmde=df["pmDE"],
        plx=df["Plx"],
        e_pmra=df["e_pmRA"],
        e_pmde=df["e_pmDE"],
        e_plx=df["e_Plx"],
    )

.. note::
    Loading the file as a `pandas.DataFrame`_  is optional, you can also pass ``numpy``
    arrays or lists/tuples for each of the data dimensions.

The resulting ``my_field`` object, assumed here to be an observed field containing a
cluster, can be processed to extract a given cluster's structural properties such as
its center coordinates and radius (see :ref:`structure`).

The :py:class:`asteca.Cluster` class also contains methods to obtain the estimated
number of members (see :ref:`nmembers`), as well as their membership probabilities
(see :ref:`membership_module`).


.. _pandas.DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html



.. _structure:

Center and radius estimation
============================

The simplest structural analysis and usually the first one to be required, is the
cluster's center estimation. **ASteCA** provides the
:py:meth:`asteca.Cluster.get_center` method to perform this estimation.

There are currently two algorithms available for center estimation. The default method
is called ``knn_5d`` and it requires that your :py:class:`asteca.Cluster` object contains
the five dimensions of data: ``(ra, dec, pmra, pmde, plx)``. The other algorithm is
called ``kde_2d`` and it requires either ``(ra, dec)`` or ``(pmra, pmde)`` 2D data,
depending on the center values you want to estimate. Calling either method is as
easy as:

.. code-block:: python

    # Will apply the default `knn_5d` algorithm on `(ra, dec, pmra, pmde, plx)` data
    my_field.get_center()

    # Will apply the `kde_2d` algorithm on `(ra, dec)` data
    my_field.get_center('kde_2d')

    # Will apply the `kde_2d` algorithm on `(pmra, pmde)` data
    my_field.get_center('kde_2d', data_2d='pms')

The estimated values will be stored as attributes of the ``my_field`` object which can
be accessed as:

.. code-block:: python

    my_field.radec_c  # (ra, dec) center coordinates
    my_field.pms_c    # proper motions center coordinates
    my_field.plx_c    # parallax center value

These attributes can also be manually set if necessary, for example:

.. code-block:: python

    my_field.radec_c = (127.3, -3.7)

.. important::

    A method to estimate the cluster's radius  will be added in future versions.
    In the meantime you can manually add the attribute with:

    .. code-block:: python

        my_field.radius = 0.1

    where the value is always assumed to be in **units of degrees**.

See the :ref:`structure_ntbk` tutorial for a step-by-step example of how to estimate the
cluster's center and radius.



.. _nmembers:

Number of members
=================

Estimating the number of members for a given cluster is a crucial previous step before
the membership analysis. Currently **ASteCA** integrates two methods to perform this
estimation, as shown in :py:meth:`asteca.Cluster.get_nmembers`.

There are currently two algorithms available for the estimation of the true number
of members associated to a give cluster. The default algorithm is called ``ripley``, it was originally introduced along with the ``fastMP`` membership method
in `Perren et al. (2023) <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__.
It requires ``(ra, dec, pmra, pmde, plx)`` data and their center estimates
(``pms_c, plx_c``).

The other algorithm is called ``density`` and it is a simple algorithm that counts the
number of stars within the cluster region (center+radius) and subtracts the expected
number of field stars within that region. It requires the ``(ra, dec)`` center and the
radius of the cluster to be defined.

Again, the methods can be accessed from the ``my_field`` object as:

.. code-block:: python

    # Will apply the default `ripley` algorithm on `(ra, dec, pmra, pmde, plx)` data
    # This assumes that the `get_center()` method was already applied or the
    # (pms_c, plx_c) center values were stored as attributes by the user
    my_field.get_nmembers()

    # Will apply the `density` algorithm on `(ra, dec)` data. This also assumes that
    # the `get_center()` method was applied, and it also requires a `radius`
    # attribute to be set in the `my_field` object
    my_field.radius = 0.1
    my_field.get_nmembers('density')

The estimated value will be stored as an attribute of the ``my_field`` object which can
be accessed as:

.. code-block:: python

    my_field.N_cluster

If the estimated number is not a proper representation of the believed number of members
for this cluster, the user can easily modify this attribute manually with:

.. code-block:: python

    my_field.N_cluster = 300

See the :ref:`membership_module` section for more details on the membership estimation
process, and the :ref:`membership_ntbk` tutorial for an example of how this value is
used to estimate the cluster's true members.
