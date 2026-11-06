.. _paramfit_section:

Parameters fit
##############

**ASteCA** generates the synthetic clusters and quantifies how similar they are to the
observed cluster, but it does not include its own sampler: the fundamental parameters
are estimated by combining it with any sampler of the user's choosing (see this list of
Python-based `MCMC and ABC packages`_). The following tutorials show how to do this
with two different packages.

.. toctree::
   :titlesonly:

   notebooks/paramfit_pyabc
   notebooks/paramfit_pocomc


The :ref:`paramfit_ntbk` notebook uses the `pyABC`_ library (Approximate Bayesian
Computation), where the **ASteCA** likelihood is used as a distance to be minimized.

The :ref:`paramfit_pocomc_ntbk` notebook uses the `pocoMC`_ library (Preconditioned
Monte Carlo), where the **ASteCA** likelihood is turned into a log-likelihood.

Both tutorials show how to set up the fitting process, run it, and extract the
parameter estimates (medians and standard deviations) from the results. They finish by
assessing the goodness of fit of the synthetic cluster generated with the median
parameters via its CMD overlap with the observed cluster (see :ref:`cmd_overlap`), and
a comparison Color-Magnitude Diagram (CMD).

pyABC versus pocoMC
-------------------

**pyABC** uses the distance returned by :py:meth:`asteca.Likelihood.get` directly,
with no likelihood scaling to choose. Its uncertainties are usually somewhat wider,
which in our tests with simulated clusters gave better calibrated error bars for
clusters with many members. On the other hand it is slower, its results depend on
when the run is stopped (maximum time, number of populations, or minimum distance),
it needs a database file, and its runs are only reproducible with the (slower) single
core sampler.

**pocoMC** was several times faster in our tests (about ten times on real clusters),
stops on its own once the posterior is reached, and is reproducible through its
``random_state`` argument. It also avoided some high ``Av`` and ``DR`` solutions where
pyABC occasionally got stuck. On the other hand the width of the posterior depends on
the chosen ``temperature``: with ``temperature=1`` the uncertainties can be too small
for clusters with many members, while large values give too wide uncertainties for
clusters with few members. The :ref:`paramfit_pocomc_ntbk` tutorial sets it with a
heuristic rule that depends on the number of members.


.. _MCMC and ABC packages: https://gabriel-p.github.io/pythonMCMC/
.. _pyABC: https://pyabc.readthedocs.io/
.. _pocoMC: https://pocomc.readthedocs.io/
