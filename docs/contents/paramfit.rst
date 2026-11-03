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


.. _MCMC and ABC packages: https://gabriel-p.github.io/pythonMCMC/
.. _pyABC: https://pyabc.readthedocs.io/
.. _pocoMC: https://pocomc.readthedocs.io/
