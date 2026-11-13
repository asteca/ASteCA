.. _sampler_module:

Sampler module
##############

The :py:class:`asteca.Sampler` class estimates the fundamental parameters of an
observed cluster. It takes a calibrated :py:class:`asteca.Synthetic` object, which
generates the synthetic clusters, compares them to the observed cluster with the
:ref:`plr` (see :ref:`sampler_lkl`), and samples the posterior distribution of the
parameters given uniform priors.

The method is Sequential Monte Carlo (SMC), which moves a set of particles from the
prior to the posterior by slowly increasing the weight given to the likelihood, with
persistent sampling (`Karamanis & Seljak 2025`_): the particles of all the iterations
are combined with importance weights. The implementation is adapted from the
`pocoMC`_ package (`Karamanis et al. 2022a`_, `2022b`_), keeping only the parts used
here: uniform priors, persistent sampling, and t-preconditioned Crank-Nicolson
(t-pCN) MCMC moves, whose proposals use a Student-t distribution fitted to the
particles in each iteration. The normalizing flow of pocoMC, i.e.: the
preconditioning of Preconditioned Monte Carlo, is not included (on a CPU it is slow to
train, and it gave the same results for this problem), so this is not Preconditioned
Monte Carlo but the SMC sampler that pocoMC runs with ``precondition=False``. No extra
packages are required, and for the same seed and settings the results are identical
to those of pocoMC with ``precondition=False``.

The likelihood is fixed: it is the one the sampler (and its ``temperature="auto"``
rule) was tested and calibrated with. To use a different likelihood or distance, or a
different sampler, see :ref:`own_lkl_sampler`.


Defining and running the sampler
================================

The fitted parameters and their uniform priors are given as a dictionary of
``(min, max)`` values. The parameters not included keep the default values of the
``synthcl`` object:

.. code-block:: python

    # Calibrate the Synthetic object with the observed cluster
    synthcl.calibrate(my_cluster)

    priors = {
        "met": (0.01, 0.02),
        "loga": (7.0, 9.5),
        "dm": (8.0, 10.5),
        "Av": (0.0, 2.0),
    }
    sampler = asteca.Sampler(synthcl, priors, temperature="auto")
    sampler.run(n_total=512)

The run ends when the particles sample the posterior and the effective sample size of
the posterior samples is at least ``n_total``. Calling :py:meth:`asteca.Sampler.run`
again with a larger ``n_total`` continues the run. The ``n_effective`` and
``n_active`` arguments set the number of particles (``256`` and ``128`` by default).
In our tests, these values and ``n_total=512`` took about half the time of twice as
many particles, with the same accuracy and uncertainties, while halving them again
made the uncertainties too small.

The sampler inherits the ``seed`` of the ``synthcl`` object, which also fixes the
random sampling of its synthetic clusters, so setting that seed is enough to make the
results reproducible.


Temperature
-----------

The log-likelihood (see :ref:`sampler_lkl`) is divided by the ``temperature``. A value
of ``1`` (the default) gives the posterior of the Poisson likelihood ratio, while
values larger than ``1`` widen the posterior by a factor of about
``sqrt(temperature)``. The synthetic clusters are only an approximation of a real
cluster, and for clusters with many members a ``temperature`` of ``1`` gives
uncertainties that are too small: sources of uncertainty that are not included in the
likelihood, such as the random sampling of the synthetic clusters or the mismatch
between the models and the observed data, dominate in this case. For clusters with
few members, on the other hand, values larger than ``1`` give uncertainties that are
too large and worse fits.

With ``temperature="auto"`` the value is selected from the number :math:`N` of stars in
the observed cluster:

============================  ===============
:math:`N`                     ``temperature``
============================  ===============
:math:`N \leq 150`            1
:math:`150 < N \leq 400`      2
:math:`N > 400`               4
============================  ===============

This is a **heuristic rule**, not derived from first principles: it was calibrated on
simulated clusters (so that the 16th-84th percentile range contains the true values
for close to 68% of them) and checked on a set of real clusters.


Faster runs
-----------

The particles are evaluated in batches. With `numba <https://numba.pydata.org/>`__
installed, the synthetic clusters of a batch and their log-likelihoods are obtained in
a single compiled call, which removes the Python overhead of each model (the gain is
largest for clusters with few stars). The results are identical in all cases.

.. important::

    **Without numba there is no speedup**: the synthetic clusters are generated and
    evaluated one at a time. The code works the same in both cases.

With ``parallel=True`` the models of each batch are distributed among numba's
threads: each thread generates the synthetic clusters of its models and evaluates
their log-likelihoods, which also speeds up clusters with many stars, where the cost
of each model is dominated by the actual computation. The ``parallel`` argument is
ignored if numba is not installed.

.. code-block:: python

    # Evaluate the models in parallel threads, four threads
    sampler = asteca.Sampler(synthcl, priors, parallel=True, n_threads=4)

By default (``n_threads=None``) numba's own setting is used, which is one thread per
CPU core unless changed with the ``NUMBA_NUM_THREADS`` environment variable (set
before importing ``asteca``) or with ``numba.set_num_threads()``; ``n_threads`` can
not be larger than ``NUMBA_NUM_THREADS``. If several fits are run at the same time,
give each one a fraction of the CPU cores so that they do not compete for them. Each
thread needs several models per batch to be used efficiently, so ``n_active`` (the
size of the batches) should be larger than the number of threads.


Results
=======

The :py:meth:`asteca.Sampler.summary` method returns the median, mean, standard
deviation, and 16th and 84th percentiles of each fitted parameter, as dictionaries
that can be passed directly to other methods, and the effective sample size (ESS) of
the importance weighted posterior samples used to compute them. The ESS is the number
of independent samples they are equivalent to, and it should not be much smaller than
the ``n_total`` value passed to :py:meth:`asteca.Sampler.run`:

.. code-block:: python

    summ = sampler.summary()
    print(summ["median"], summ["std"], summ["ess"])

    # Synthetic clusters sampled from the fitted solution
    synthcl.get_models(summ["median"], summ["std"])

The posterior samples are returned by :py:meth:`asteca.Sampler.posterior`, with their
importance weights, or resampled to equally weighted samples with ``resample=True``
(e.g.: for a corner plot). The columns follow the order of the ``priors`` dictionary,
stored in the ``names`` attribute.

After a run, the object also stores the best solution found (``best_params``) and its
distance (``best_dist``, see :ref:`plr`), the number of iterations (``n_iter``) and
likelihood calls (``n_calls``), the effective sample size of the posterior samples
(``ess``), and the log-evidence estimate (``logz``), which can be used to compare
models fitted to the same observed cluster.

The distance ``best_dist`` can not be compared across observed clusters. To assess
whether the fitted synthetic cluster is a good representation of the observed one,
use :py:meth:`asteca.Synthetic.cmd_overlap` (see :ref:`cmd_overlap`).

See the :ref:`paramfit_sampler_ntbk` tutorial for a complete example.



.. _sampler_lkl:

Likelihood
==========

The likelihood compares the observed and synthetic clusters through their Hess diagrams
(the color-magnitude diagram split into cells). The observed cluster and its Hess
diagram are obtained when the ``synthcl`` object is calibrated with
:py:meth:`asteca.Synthetic.calibrate`.


.. _ref_binning:

Binning
-------

The cells are defined from the observed photometry when the ``synthcl`` object is
calibrated. The number of bins in magnitude and in each color is set with Knuth's rule
(see `Choosing Histogram Bins
<https://docs.astropy.org/en/stable/visualization/histogram.html>`_ in the astropy
documentation), and the bins have equal widths and span the observed range of each
dimension.


.. _plr:

Poisson likelihood ratio
------------------------

This is the Poisson likelihood ratio as defined in `Tremmel et al (2013)`_. Starting
from Eq 10 with :math:`v_{i,j}=1` (i.e.: equal volumes for all cells), we have:

.. math::

    p(d|\theta) = \prod_i^B \frac{\Gamma(n_i+m_i+\frac{1}{2})}
    {2^{n_i+m_i+\frac{1}{2}} n_i!\Gamma(m_i+\frac{1}{2})}

where :math:`n_i` and :math:`m_i` are the number of observed and synthetic stars in
bin :math:`i`, respectively, and :math:`B` is the total number of bins.
Applying the logarithm:

.. math::

    \log(p) = \sum_i^B \left[\log\Gamma(n_i+m_i+\frac{1}{2})
    - (m_i+n_i+\frac{1}{2})\log2 -\log n_i!
    - \log \Gamma(m_i+\frac{1}{2}) \right]

re-arranging:

.. math::

    \log(p) = \sum_i^B \left[\log\Gamma(n_i+m_i+\frac{1}{2})-
    \log \Gamma(m_i+\frac{1}{2}) \right]
    - (M+N+\frac{1}{2}B) \log 2 - \sum_i^B \log n_i!

employing the :math:`LogGamma()` function;

.. math::

    \log(p) = \sum_i^B \left[LogGamma(n_i+m_i+\frac{1}{2}) - LogGamma(m_i+\frac{1}{2}) \right] - \\ (N+\frac{1}{2}B) \log 2 - \sum_i^B \log n_i! -M  \log 2

.. note::

    The :math:`LogGamma()` function is used instead of the :math:`\Gamma()`
    function because it is numerically more stable, especially when
    dealing with large numbers. This is implemented in `scipy.special.gammaln`_.


We can discard the second and third terms as they depend on :math:`N` and :math:`B`
only (i.e.: the number of stars in observed cluster and the number of bins), and hence
do not change for different synthetic data.

The last term depends on :math:`M`, the total number of synthetic stars. Once the
:py:class:`asteca.Synthetic` object is calibrated, every synthetic cluster generated
by :py:meth:`asteca.Synthetic.generate` contains the same number of stars as the
observed cluster (with valid photometry), i.e.: :math:`M=N` (except in the case
described below), and this term is a constant that can also be discarded.

.. important::

    A synthetic cluster can contain fewer stars than the observed one
    (:math:`M<N`) only when the IMF sample (see the ``max_mass`` argument of
    :py:class:`asteca.Synthetic`) can not provide enough stars brighter than the
    maximum observed magnitude. This only happens for models far from the observed
    cluster (e.g.: too distant or too old), which have a large distance with or
    without this term.

Finally, discarding the last three terms, we have the expression for the Poisson
likelihood ratio:

.. math::

    \log(p)\approx \sum_i^B \left[LogGamma(n_i+m_i+\frac{1}{2}) - LogGamma(m_i+\frac{1}{2}) \right]


This likelihood is normalized and inverted into a distance, so that the best value is
``0`` and the worst value is ``1``:

.. math::

    d = \frac{\log(p)_{max} - \log(p)}{\log(p)_{max} - \log(p)_{min}}

where :math:`\log(p)_{max}` is the value obtained comparing the observed cluster with
itself (:math:`m_i=n_i`), and :math:`\log(p)_{min}` is the value obtained when no
synthetic star falls in any of the bins occupied by observed stars (:math:`m_i=0`).
Hence :math:`d=0` means an exact match, and :math:`d=1` means that the synthetic and
observed clusters do not overlap at all in the CMD. This puts the end points of the
distance on the same scale for all clusters, but the values in between are still not
comparable across clusters: the distance obtained for a correct model depends on the
number of stars and on the binning (see :ref:`cmd_overlap`).

When the synthetic and observed clusters do not overlap at all, the likelihood takes
the same value regardless of where the synthetic stars are located in the CMD. To
break this tie, which can stall the parameter fit when most models sampled from the
priors have no overlap (e.g.: for small clusters), the distance for these models is
set to:

.. math::

    d = 1 + \left(\frac{\tilde{G}_{obs} - \tilde{G}_{syn}}{\Delta G}\right)^2 +
    \left(\frac{\tilde{C}_{obs} - \tilde{C}_{syn}}{\Delta C}\right)^2

where :math:`\tilde{G}` and :math:`\tilde{C}` are the median magnitude and color of
the observed and synthetic clusters, and :math:`\Delta G` and :math:`\Delta C` are the
widths of the magnitude and color ranges of the observed Hess diagram. The latter
put both dimensions on the same scale, so that the magnitude (which usually spans a
much larger range) does not dominate over the color. Hence any synthetic cluster with
some overlap (:math:`d<1`) is always considered better than any synthetic cluster with
no overlap (:math:`d\geq 1`), and the latter are pulled towards the region of the
CMD where the observed cluster is located.

Synthetic clusters with no stars (e.g.: parameters that place the whole isochrone below
the maximum magnitude of the observed cluster) are assigned a fixed distance of
:math:`d=10`, so that they are always considered worse than synthetic clusters with
stars (unless these have medians several Hess diagram widths away from the observed
cluster).

The log-likelihood used by the sampler is the one that corresponds to this distance,
relative to its maximum value and divided by the ``temperature`` :math:`T`:

.. math::

    \log L = -d\,\frac{\log(p)_{max} - \log(p)_{min}}{T}

so that ``0`` is again the best value. For synthetic clusters that overlap the
observed one this is simply :math:`(\log(p) - \log(p)_{max})/T`.


.. _pocoMC: https://pocomc.readthedocs.io/
.. _Karamanis et al. 2022a: https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.1644K
.. _2022b: https://arxiv.org/abs/2207.05660
.. _Karamanis & Seljak 2025: https://arxiv.org/abs/2407.20722
.. _Tremmel et al (2013): https://ui.adsabs.harvard.edu/abs/2013ApJ...766...19T/abstract
.. _scipy.special.gammaln: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaln.html
