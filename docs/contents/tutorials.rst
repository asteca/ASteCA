.. _tutorials:

Tutorials
##########

This section contains a set of tutorials designed to guide the user through different
types of analyses, e.g.: analyzing the structure of the cluster, estimating the
cluster's membership probabilities, generating synthetic clusters,
and estimating a cluster's fundamental parameters.

.. toctree::
   :titlesonly:

   notebooks/structure
   notebooks/membership
   notebooks/synthetic
   paramfit
   notebooks/post_process


The :ref:`structure_ntbk` notebook serves as a guide to the initial structural
characterization of a stellar cluster field. It covers loading data,
visualizing distributions, and applying built-in methods to estimate the cluster center
and the number of likely members, while also pointing out parameters like the radius that
currently require manual input or external estimation.


The :ref:`membership_ntbk` notebook serves as a practical guide to applying the
:py:meth:`asteca.Membership.fastmp` and :py:meth:`asteca.Membership.bayesian` membership
probability estimation methods. It highlights the different data requirements (5D for
:py:obj:`fastmp`, photometric for :py:obj:`bayesian`) and necessary preliminary parameter
estimations (center, ``N_cluster``, radius for :py:obj:`bayesian`) for each method,
concluding with a visual comparison of their results.


The :ref:`synth_ntbk` notebook demonstrates the process of generating synthetic star
clusters via the :py:meth:`asteca.Synthetic.generate` method using specified fundamental
parameters. Optionally, the generator can be calibrated with observed cluster data using
the :py:meth:`asteca.Synthetic.calibrate` method to enhance realism in terms of star count, magnitude limits, and errors. The process culminates in visualizing the resulting synthetic Color-Magnitude Diagrams (CMDs), often distinguishing stellar systems and potentially overlaying the underlying isochrone.

The :ref:`paramfit_section` section contains two notebooks that provide a
step-by-step guide to estimating the fundamental parameters of a stellar cluster,
with the `pyABC`_ library and with the :py:class:`asteca.Sampler` class (adapted from
the `pocoMC`_ library). They show how to set
up the fitting process, run it, and interpret the results, including the final
parameter estimates (medians and standard deviations), the goodness of fit (CMD
overlap), and a comparison Color-Magnitude Diagram (CMD).


The :ref:`post_process_ntbk` notebook shows how to apply the following methods:

- :meth:`asteca.Synthetic.get_models`: sample a set of synthetic clusters
  with the same parameters as the observed cluster (previously calibrated and fitted)
- :meth:`asteca.Synthetic.cluster_masses`: determine various associated cluster masses
  (initial, actual, dynamical, etc.)
- :meth:`asteca.Synthetic.stellar_parameters`: estimate per-star masses and
  probabilities of being binary systems and blue straggler stars



.. _pyABC: https://pyabc.readthedocs.io/
.. _pocoMC: https://pocomc.readthedocs.io/