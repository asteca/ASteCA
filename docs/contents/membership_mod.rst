.. _membership_module:

Membership module
#################

The :py:class:`asteca.Membership` class allows estimating the membership probabilities
for all the stars in a given observed field. There are currently two methods included in
this class: :py:meth:`asteca.Membership.bayesian` and
:py:meth:`asteca.Membership.fastmp`.

The :py:meth:`bayesian` method was described in detail in the `article`_ where **ASteCA**
was originally introduced. The method requires ``(ra, dec)``  data and will use any extra
data dimensions stored in the :py:class:`Cluster <asteca.cluster.Cluster>` object, i.e.:
photometry, proper motions, and parallax. A minimum of two data dimensions are required,
in addition to ``(ra, dec)``. This method can produce membership probabilities on
photometric data alone.

The :py:meth:`fastmp` method was described in detail in the
`article <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__
where the `Unified Cluster Catalogue (UCC) <https://ucc.ar/>`__ was introduced. The
method requires proper motions, and parallax data dimensions stored in the
:py:class:`Cluster <asteca.cluster.Cluster>` object. Photometric data is not employed.

.. important::
    The only advantage of the :py:meth:`bayesian` method over the :py:meth:`fastmp`
    method is that the former works with photometric data. Hence it should only be used
    in cases were only photometric data is available, as :py:meth:`fastmp` is not only
    much faster but also more precise in those cases where proper motions and/or
    parallax data is available.

To use these methods we need to estimate the cluster's number of members as described in
the :ref:`nmembers` section, which is done by calling the
:py:meth:`asteca.Cluster.get_nmembers` method.

With the ``N_cluster`` attribute in place in a
:py:class:`Cluster <asteca.cluster.Cluster>` object, here called ``my_field``, you can
define a :py:class:`Membership <asteca.membership.Membership>` object, here called
:py:obj:`memb`:

.. code-block:: python

    # Define a `membership` object
    memb = asteca.Membership(my_field)

and apply either the :py:meth:`bayesian` or the :py:meth:`fastmp` method:


.. code-block:: python

    # Run `fastmp` method
    probs_fastmp = memb.fastmp()

    # Run `bayesian` method
    probs_bayes = memb.bayesian()

The arrays stored in the ``probs_fastmp`` or ``probs_bayes`` variables are the
per-star membership probabilities. The results will naturally not be equivalent as both
algorithms are rather different. The :py:meth:`bayesian` algorithm for example tends to
assign lower probabilities than the :py:meth:`fastmp` algorithm.

A step-by-step example is shown in the :ref:`membership_ntbk` tutorial.


.. _article: https://doi.org/10.1051/0004-6361/201424946