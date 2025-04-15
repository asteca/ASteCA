.. _tutorials:

Tutorials
##########

This section contains a set of tutorials designed to guide the user through different
types of analyses, e.g.: analyzing the structure of the cluster, estimating the
cluster's membership probabilities, generating synthetic clusters,
and estimating a cluster's fundamental parameters.

.. toctree::
   :titlesonly:

   structure
   membership
   synthetic
   paramfit
   masses_bfr


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

The :ref:`paramfit_ntbk` notebook provides a step-by-step guide to estimating the
fundamental parameters of a stellar cluster using **ASteCA** and the `pyABC`_
library. The tutorial shows how to set up the fitting process, run it, and interpret
the results, including the final parameter estimates (medians and standard deviations)
along with ``pyABC``  visualizations and a comparison Color-Magnitude Diagram (CMD).


Given an observed cluster the :ref:`masses_bfr_ntbk` notebook shows how to apply
the following methods:

- :meth:`asteca.Synthetic.get_models`: sample a set of synthetic clusters
  with the same parameters as the observed cluster (previously calibrated and fitted)
- :meth:`asteca.Synthetic.stellar_masses`: estimate its individual stellar masses and
  their probabilities of being binary systems 
- :meth:`asteca.Synthetic.binary_fraction`: calculate the overall binary fraction of
  the cluster
- :meth:`asteca.Synthetic.cluster_masses`: determine various associated cluster masses
  (initial, actual, observed, etc.)



.. _pyABC: https://pyabc.readthedocs.io/