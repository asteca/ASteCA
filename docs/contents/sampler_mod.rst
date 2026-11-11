.. _sampler_module:

Sampler module
##############

The :py:class:`asteca.Sampler` class estimates the fundamental parameters of an
observed cluster. It combines a calibrated :py:class:`asteca.Synthetic` object, which
generates the synthetic clusters, with a :py:class:`asteca.Likelihood` object, which
compares them to the observed cluster, and samples the posterior distribution of the
parameters given uniform priors.

The method is Preconditioned Monte Carlo (PMC), a Sequential Monte Carlo method that
moves a set of particles from the prior to the posterior by slowly increasing the
weight given to the likelihood. The implementation is adapted from the `pocoMC`_
package (`Karamanis et al. 2022a`_, `2022b`_), keeping only the parts used here:
uniform priors, persistent sampling (`Karamanis & Seljak 2025`_), and t-preconditioned
Crank-Nicolson MCMC moves, without the normalizing flow (which on a CPU is slow to
train and gave the same results for this problem). No extra packages are required.
For the same seed and settings the results are identical to those of pocoMC with
``precondition=False``.


Defining and running the sampler
================================

The fitted parameters and their uniform priors are given as a dictionary of
``(min, max)`` values. The parameters not included keep the default values of the
``synthcl`` object:

.. code-block:: python

    # Calibrated Synthetic object and Likelihood object
    synthcl.calibrate(my_cluster)
    likelihood = asteca.Likelihood(my_cluster)

    priors = {
        "met": (0.01, 0.02),
        "loga": (7.0, 9.5),
        "dm": (8.0, 10.5),
        "Av": (0.0, 2.0),
    }
    sampler = asteca.Sampler(synthcl, likelihood, priors, temperature="auto", seed=1)
    sampler.run(n_total=512)

The run ends when the particles sample the posterior and the effective sample size of
the posterior samples is at least ``n_total``. Calling :py:meth:`asteca.Sampler.run`
again with a larger ``n_total`` continues the run. The ``n_effective`` and
``n_active`` arguments set the number of particles (``256`` and ``128`` by default).
In our tests, these values and ``n_total=512`` took about half the time of twice as
many particles, with the same accuracy and uncertainties, while halving them again
made the uncertainties too small.

The results depend on the ``seed`` of the sampler and also on the seed of the
``synthcl`` object, which fixes the random sampling of its synthetic clusters.


Temperature
-----------

The log-likelihood is obtained with :py:meth:`asteca.Likelihood.log_lkl`, divided by
the ``temperature``. A value of ``1`` (the default) gives the posterior of the Poisson
likelihood ratio, while values larger than ``1`` widen the posterior by a factor of
about ``sqrt(temperature)``. The synthetic clusters are only an approximation of a real
cluster, and for clusters with many members a ``temperature`` of ``1`` gives
uncertainties that are too small. For clusters with few members, on the other hand,
values larger than ``1`` give uncertainties that are too large and worse fits.

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

The particles are evaluated in batches with
:py:meth:`asteca.Likelihood.log_lkl_batch` (see :ref:`ref_log_lkl_batch`). With
`numba <https://numba.pydata.org/>`__ installed this removes the Python overhead of
each synthetic cluster, and with ``parallel=True`` the models of each batch are also
distributed among the CPU cores (the number of threads can be set with
``n_threads``). **Without numba there is no speedup**, and the results are the same in
all cases.


Results
=======

The :py:meth:`asteca.Sampler.summary` method returns the median, mean, standard
deviation, and 16th and 84th percentiles of each fitted parameter, as dictionaries
that can be passed directly to other methods:

.. code-block:: python

    summ = sampler.summary()
    print(summ["median"], summ["std"])

    # Synthetic clusters sampled from the fitted solution
    synthcl.get_models(summ["median"], summ["std"])

The posterior samples are returned by :py:meth:`asteca.Sampler.posterior`, with their
importance weights, or resampled to equally weighted samples with ``resample=True``
(e.g.: for a corner plot). The columns follow the order of the ``priors`` dictionary,
stored in the ``names`` attribute.

After a run, the object also stores the best solution found (``best_params``) and its
distance (``best_dist``, see :py:meth:`asteca.Likelihood.dist`), the number of
iterations (``n_iter``) and likelihood calls (``n_calls``), the effective sample size
of the posterior samples (``ess``), and the log-evidence estimate (``logz``), which can
be used to compare models fitted to the same observed cluster.

See the :ref:`paramfit_pocomc_ntbk` tutorial for a complete example.


.. _pocoMC: https://pocomc.readthedocs.io/
.. _Karamanis et al. 2022a: https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.1644K
.. _2022b: https://arxiv.org/abs/2207.05660
.. _Karamanis & Seljak 2025: https://arxiv.org/abs/2407.20722
