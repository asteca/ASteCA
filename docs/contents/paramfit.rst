.. _paramfit_section:

Parameters fit
##############

**ASteCA** includes its own sampler, :py:class:`asteca.Sampler` (see
:ref:`sampler_module`), which estimates the fundamental parameters with Preconditioned
Monte Carlo. The synthetic clusters and the likelihood can also be combined with any
other sampler of the user's choosing (see this list of Python-based
`MCMC and ABC packages`_). The following tutorials show both options.

.. toctree::
   :titlesonly:

   notebooks/paramfit_pyabc
   notebooks/paramfit_pocomc


The :ref:`paramfit_ntbk` notebook uses the `pyABC`_ library (Approximate Bayesian
Computation), where the **ASteCA** likelihood is used as a distance to be minimized.

The :ref:`paramfit_pocomc_ntbk` notebook uses the :py:class:`asteca.Sampler` class,
adapted from the `pocoMC`_ library (Preconditioned Monte Carlo), where the **ASteCA**
likelihood is turned into a log-likelihood.

Both tutorials show how to set up the fitting process, run it, and extract the
parameter estimates (medians and standard deviations) from the results. They finish by
assessing the goodness of fit of the synthetic cluster generated with the median
parameters via its CMD overlap with the observed cluster (see :ref:`cmd_overlap`), and
a comparison Color-Magnitude Diagram (CMD).

pyABC versus asteca.Sampler
---------------------------

**pyABC** uses the distance returned by :py:meth:`asteca.Likelihood.dist` directly,
with no likelihood scaling to choose. Its uncertainties are usually somewhat wider,
which in our tests with simulated clusters gave better calibrated error bars for
clusters with many members. On the other hand it is slower, its results depend on
when the run is stopped (maximum time, number of populations, or minimum distance),
it needs a database file, and its runs are only reproducible with the (slower) single
core sampler.

**asteca.Sampler** (as pocoMC, from which it is adapted) was several times faster in
our tests (about ten times on real clusters), stops on its own once the posterior is
reached, and is reproducible through its ``seed`` argument. It also avoided some high ``Av`` and ``DR`` solutions where
pyABC occasionally got stuck. On the other hand the width of the posterior depends on
the chosen ``temperature``: with ``temperature=1`` the uncertainties can be too small
for clusters with many members, while large values give too wide uncertainties for
clusters with few members. With ``temperature="auto"`` it is set with a heuristic rule
that depends on the number of members (see :ref:`sampler_module`).


.. _MCMC and ABC packages: https://gabriel-p.github.io/pythonMCMC/
.. _pyABC: https://pyabc.readthedocs.io/
.. _pocoMC: https://pocomc.readthedocs.io/
