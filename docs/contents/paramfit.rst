.. _paramfit_section:

Parameters fit
##############

**ASteCA** includes its own sampler, :py:class:`asteca.Sampler` (see
:ref:`sampler_module`), which estimates the fundamental parameters with Sequential
Monte Carlo and a built-in likelihood. The synthetic clusters can also be combined
with a likelihood (or distance) and a sampler of the user's choosing (see
:ref:`own_lkl_sampler` and this list of Python-based `MCMC and ABC packages`_). The
following tutorials show both options.

.. toctree::
   :titlesonly:

   notebooks/paramfit_sampler
   notebooks/paramfit_pyabc
   notebooks/paramfit_pocomc
   notebooks/paramfit_gof


The :ref:`paramfit_sampler_ntbk` notebook uses the :py:class:`asteca.Sampler` class
(Sequential Monte Carlo, adapted from the `pocoMC`_ library), with the likelihood
built into **ASteCA** (see :ref:`sampler_lkl`). This is the recommended option.

The :ref:`paramfit_ntbk` notebook uses the `pyABC`_ library (Approximate Bayesian
Computation), with a simple distance defined by the user, which pyABC minimizes.

The :ref:`paramfit_pocomc_ntbk` notebook uses the `pocoMC`_ library itself, with a
Poisson likelihood defined by the user.

All the tutorials show how to set up the fitting process, run it, and extract the
parameter estimates (medians and standard deviations) from the results.

The :ref:`paramfit_gof_ntbk` notebook shows how to assess the goodness of fit of the
synthetic cluster generated with the fitted parameters (whichever method was used to
obtain them) via its CMD overlap with the observed cluster (see :ref:`cmd_overlap`),
and a comparison Color-Magnitude Diagram (CMD).

pyABC versus asteca.Sampler
---------------------------

**pyABC** minimizes a distance directly, with no likelihood scaling to choose. In our
tests, where the distance used was the one that corresponds to the likelihood of
:py:class:`asteca.Sampler` (see :ref:`plr`) instead of the simple one used in the
tutorial, its uncertainties were usually somewhat wider, which with simulated clusters
gave better calibrated error bars for clusters with many members. These results
depend on the distance used. On the other hand it is slower, its results depend on
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
