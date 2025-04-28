:orphan:

.. _changelog:


Changelog
#########


`[v0.6.1] <https://github.com/asteca/asteca/releases/tag/v0.6.1>`__ - 2025-04-28
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Isochrones module:
   - Fixed improper MIST isochrone generation (https://github.com/asteca/ASteCA/issues/569) via EEP (Equivalent Evolutionary Point) interpolation
   - Changed to per-isochrone mass based interpolation for other services

- Documentation & Code Structure:
   - Fixed LaTeX rendering issues
   - Added umami tracker (https://github.com/asteca/ASteCA/issues/571)



`[v0.6.0] <https://github.com/asteca/asteca/releases/tag/v0.6.0>`__ - 2025-04-15
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Downgraded the versions of almost all the required packages to play nice with
  Google Colab (again)
- Faster import of `asteca`, using lazy imports in classes
  (https://github.com/asteca/ASteCA/issues/564); removed the `plot` class,
  removed `pandas` and `matplotlib` as dependencies
- Allow generating synthetic clusters without calibrating an observed cluster
  (https://github.com/asteca/ASteCA/issues/565>)
- Remove `fix_params` parameter from `synthetic.calibrate()`
- Major changes to the documentation
- Small performance improvements in the generation of synthetic cluster
  (zaWAverage function)



`[v0.5.9] <https://github.com/asteca/asteca/releases/tag/v0.5.9>`__ - 2025-01-24
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Fix issue with BASTI isochrones that would store ages out of order
- Upgrade versions, Google Colab will most likely not work due to running older
  versions of Python and the required packages.
- Remove new version checker https://github.com/asteca/ASteCA/issues/560
- Many type hints added to the codebase



`[v0.5.8] <https://github.com/asteca/asteca/releases/tag/v0.5.8>`__ - 2024-08-09
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Downgraded the versions of almost all the required packages to play nice with
  Google Colab
- Fix `AttributeError: ptp was removed from the ndarray class in NumPy 2.0. Use
  np.ptp(arr, ...) instead.` for newer versions of `numpy`
- Allow `isochrones` to load single files, not just paths to folders
- Added an update notifier to `__init__.py` (added `requests` as requirement)



`[v0.5.7] <https://github.com/asteca/asteca/releases/tag/v0.5.7>`__ - 2024-07-31
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Faster (and cleaner) :py:meth:`asteca.membership.Membership.fastmp` and
  :py:meth:`asteca.membership.Membership.bayesian`
- Change the normalization constant used in ``fastmp`` (from 2*median to IQR)
- Add ``kde_2d`` method to :py:meth:`asteca.cluster.Cluster.get_center`
- Fix ``knn_5d`` method to better handle large frames
- Minor fixes to center and total masses estimation



`[v0.5.6] <https://github.com/asteca/asteca/releases/tag/v0.5.6>`__ - 2024-06-25
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Simpler error distribution function
- Convert `dataclass` to `class` (simpler API documentation)



`[v0.5.5] <https://github.com/asteca/asteca/releases/tag/v0.5.5>`__ - 2024-06-21
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Add :py:meth:`asteca.membership.Membership.bayesian` method decontamination algorithm
- Move plot functions to its own module
- Removed ``masses_binary_probs`` method from the
  :py:class:`asteca.synthetic.Synthetic` class
  and added the methods: :py:meth:`asteca.synthetic.Synthetic.get_models`,
  :py:meth:`asteca.synthetic.Synthetic.stellar_masses`,
  :py:meth:`asteca.synthetic.Synthetic.binary_fraction`,
  :py:meth:`asteca.synthetic.Synthetic.cluster_masses`, to estimate individual masses,
  total binarity fraction, and total masses (see: :ref:`masses_and_binarity`)
- Added the :py:meth:`asteca.cluster.Cluster.get_nmembers` method to estimate the
  number of members for a cluster embedded in a field.



`[v0.5.4] <https://github.com/asteca/asteca/releases/tag/v0.5.4>`__ - 2024-05-26
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Fix small issue with documentation in v0.5.3



`[v0.5.3] <https://github.com/asteca/asteca/releases/tag/v0.5.3>`__ - 2024-05-26
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Added the :py:meth:`asteca.membership.Membership.fastmp` method
- Added the :py:meth:`asteca.cluster.Cluster.get_center` method



`[v0.5.2] <https://github.com/asteca/asteca/releases/tag/v0.5.2>`__ - 2024-05-07
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Re-wrote isochrones reading module. See :ref:`isochrones_module`.
- Added `Gaia (E)DR3 extinction law <https://www.cosmos.esa.int/web/gaia/edr3-extinction-law>`_.
  See :ref:`synthetic_module`.
- Downgrade requirements to ``python=3.10`` to allow running in default Google Colab
- Fixed wrong masses and binary probabilities being assigned to stars with missing
  photometric data
- Fixed bug in :class:`likelihood` when two colors were defined
- Move ``z_to_FeH`` argument from :class:`synthetic` to :class:`isochrones`
- Remove ``min_max()`` method from :class:`synthetic`, values are now shown when
  instantiating an :class:`isochrones` class
- Allow plotting two colors in :class:`cluster` and :class:`synthetic` plotting modules



`[v0.5.1] <https://github.com/asteca/asteca/releases/tag/v0.5.1>`__ - 2024-04-19
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Skip hidden files when reading isochrones (https://github.com/asteca/ASteCA/issues/551)



`[v0.5.0] <https://github.com/asteca/asteca/releases/tag/v0.5.0>`__ - 2024-04-18
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Major update. **ASteCA** is now a proper Python package that can be installed with
`pip install asteca`. Some of the most important changes are:

- Added support for MIST and Basti isochrones, in addition to PARSEC
- No longer depends on `ptemcee <https://github.com/willvousden/ptemcee>`_, now the
  user can employ any package to perform the parameter inference
- Updated binary systems generation
- Removed structural analysis function (for now at least)


[v0.4.4-0.4.9]
++++++++++++++

These version numbers were skipped due to the major changes introduced in version
``0.5.0``.



`[v0.4.3] <https://github.com/asteca/asteca/releases/tag/v0.4.3>`__ - 2022-06-01
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This is a large update mostly on the synthetic cluster generation process. Two
extra free fundamental parameters are added: differential reddening and the
ratio of total to selective absorption.

- Fit Av instead of E(B-V) (`544 <https://github.com/asteca/ASteCA/issues/544>`__)
- Clean up before releasing 0.4.3 (`541 <https://github.com/asteca/ASteCA/issues/541>`__)
- Change the 'completeness' function (`540 <https://github.com/asteca/ASteCA/issues/540>`__)
- Structural analysis consumes too much memory for large fields (`539 <https://github.com/asteca/ASteCA/issues/539>`__)
- Remove mass parameter (`538 <https://github.com/asteca/ASteCA/issues/538>`__)
- Use pandas to read cluster data (`537 <https://github.com/asteca/ASteCA/issues/537>`__)
- Remove complete/incomplete codebase + simplify get_data() (`533 <https://github.com/asteca/ASteCA/issues/533>`__)
- Change several input parameters (`532 <https://github.com/asteca/ASteCA/issues/532>`__)
- Fix King profile fit (`529 <https://github.com/asteca/ASteCA/issues/529>`__)
- Add a flag to turn completeness on/off (`528 <https://github.com/asteca/ASteCA/issues/528>`__)
- Make the lower mass limits input parameters (`526 <https://github.com/asteca/ASteCA/issues/526>`__)
- Remove (RA, DE) transformation, RV mentions, trim frame, error rejection (`524 <https://github.com/asteca/ASteCA/issues/524>`__)
- Optimal radius estimation fails for some clusters (`521 <https://github.com/asteca/ASteCA/issues/521>`__)
- Gap in the synthetic cluster (`519 <https://github.com/asteca/ASteCA/issues/519>`__)
- DA 'read' mode: read from column (match with pyUPMASK) (`518 <https://github.com/asteca/ASteCA/issues/518>`__)
- Simplify Bayesian parallax inference (`516 <https://github.com/asteca/ASteCA/issues/516>`__)
- Per cluster fundamental parameters range (`514 <https://github.com/asteca/ASteCA/issues/514>`__)
- Revise default parallax offset following eDR3 release (`501 <https://github.com/asteca/ASteCA/issues/501>`__)
- Generalize mass-ratio distribution for binaries using a power law (`496 <https://github.com/asteca/ASteCA/issues/496>`__)
- Output of cluster memberships is very slow for large clusters (`437 <https://github.com/asteca/ASteCA/issues/437>`__)
- Generate finding chart plot (`210 <https://github.com/asteca/ASteCA/issues/210>`__)
- Probability density for binary assignment (`198 <https://github.com/asteca/ASteCA/issues/198>`__)
- Differential reddening (`174 <https://github.com/asteca/ASteCA/issues/174>`__)
- Make extinction parameter Rv a free parameter (`170 <https://github.com/asteca/ASteCA/issues/170>`__)



`[v0.4.2] <https://github.com/asteca/asteca/releases/tag/v0.4.2>`__ - 2021-05-10
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Fixed two issues: don't read hidden files from the ``input/`` folder, remove
  forgotten parameter that was removed in the previous release.


`[v0.4.1] <https://github.com/asteca/asteca/releases/tag/v0.4.1>`__ - 2021-05-05
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Fixed estimated optimal radius that was too large
   (`513 <https://github.com/asteca/ASteCA/issues/513>`__)
-  Deprecate pixel coordinate support
   (`509 <https://github.com/asteca/ASteCA/issues/509>`__)
-  Coordinates density map shows artifact in corners
   (`511 <https://github.com/asteca/ASteCA/issues/511>`__)
-  Split A-D test into one test per feature
   (`477 <https://github.com/asteca/ASteCA/issues/477>`__)


`[v0.4.0] <https://github.com/asteca/asteca/releases/tag/v0.4.0>`__ - 2021-05-03
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Compensate cluster’s mass for binaries masses?
   (`488 <https://github.com/asteca/ASteCA/issues/488>`__)
-  Estimate individual per-star masses
   (`484 <https://github.com/asteca/ASteCA/issues/484>`__)
-  Improve performance of synth cluster generation (3)
   (`506 <https://github.com/asteca/ASteCA/issues/506>`__)
-  Simplify isochrones download/handling
   (`497 <https://github.com/asteca/ASteCA/issues/497>`__)
-  Add CS 37 COLIBRI track + deprecate old versions 10 & 11 of PARSEC
   (`495 <https://github.com/asteca/ASteCA/issues/495>`__)
-  Optimal radius too large for some clusters
   (`510 <https://github.com/asteca/ASteCA/issues/510>`__)
-  Project equatorial coordinates before processing
   (`237 <https://github.com/asteca/ASteCA/issues/237>`__)
-  Add eccentricity parameter to KP fit?
   (`480 <https://github.com/asteca/ASteCA/issues/480>`__)
-  Finish working on enhanced King profile fitting
   (`456 <https://github.com/asteca/ASteCA/issues/456>`__)
-  Remove KDE_stds and mp_flag parameters
   (`500 <https://github.com/asteca/ASteCA/issues/500>`__)
-  Simplify input of structure parameters
   (`512 <https://github.com/asteca/ASteCA/issues/512>`__)
-  Deprecate all likelihoods except Tremmel
   (`507 <https://github.com/asteca/ASteCA/issues/507>`__)
-  Interpolate IMF masses into the isochrones, not the other way around
   (`503 <https://github.com/asteca/ASteCA/issues/503>`__)
-  Add minimum binary mass ratio to fundamental parameters?
   (`504 <https://github.com/asteca/ASteCA/issues/504>`__)
-  Deprecate Anderson-Darling test
   (`499 <https://github.com/asteca/ASteCA/issues/499>`__)
-  Deprecate “Read mode”
   (`498 <https://github.com/asteca/ASteCA/issues/498>`__)
-  Add IMF and PMF curves obtention
   (`96 <https://github.com/asteca/ASteCA/issues/96>`__)
-  Convert pixel coordinates to RA & DEC
   (`203 <https://github.com/asteca/ASteCA/issues/203>`__)
-  Add ZAMS to CMD final plot
   (`160 <https://github.com/asteca/ASteCA/issues/160>`__)
-  Add semi_input.dat checking to checker
   (`214 <https://github.com/asteca/ASteCA/issues/214>`__)
-  Add weighted spatial density map
   (`167 <https://github.com/asteca/ASteCA/issues/167>`__)
-  Generate output CMD-CCD plots for the mean+median+mode
   (`479 <https://github.com/asteca/ASteCA/issues/479>`__)
-  Exact circle area using geometry instead of Monte Carlo
   (`446 <https://github.com/asteca/ASteCA/issues/446>`__)
-  Use the maximum number of members in the optimal radius?
   (`494 <https://github.com/asteca/ASteCA/issues/494>`__)
-  Add 1-sigma region to King profile
   (`478 <https://github.com/asteca/ASteCA/issues/478>`__)
-  Turn off MP coloring in D2 plots for binned likelihoods
   (`473 <https://github.com/asteca/ASteCA/issues/473>`__)


`[v0.3.1] <https://github.com/asteca/asteca/releases/tag/v0.3.1>`__ - 2020-06-19
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Only the ``ptemcee`` method is kept, all others are now deprecated.

-  Corrected an error in the ``CMD_systs.dat`` file
   (`468 <https://github.com/asteca/ASteCA/issues/468>`__)
-  `Fixed
   path <https://github.com/asteca/ASteCA/commit/3ab2b30d3d107972734112e7f0bd8ce12709ebdc>`__
   for ``CMD_systs.dat``, now works in Windows (and Mac?)
-  Control (some) plotting parameters through custom style and allow the
   selection of one of the supported styles
   (`464 <https://github.com/asteca/ASteCA/issues/464>`__)
-  Dump the results of the fundamental parameters analysis to file
   (`467 <https://github.com/asteca/ASteCA/issues/467>`__)
-  Closed several issues related to the deprecated bootstrap(+GA), brute
   force, and emcee methods
   (`265 <https://github.com/asteca/ASteCA/issues/265>`__,
   `280 <https://github.com/asteca/ASteCA/issues/280>`__,
   `284 <https://github.com/asteca/ASteCA/issues/284>`__,
   `324 <https://github.com/asteca/ASteCA/issues/324>`__,
   `341 <https://github.com/asteca/ASteCA/issues/341>`__,
   `347 <https://github.com/asteca/ASteCA/issues/347>`__,
   `418 <https://github.com/asteca/ASteCA/issues/418>`__,
   `442 <https://github.com/asteca/ASteCA/issues/442>`__,
   `447 <https://github.com/asteca/ASteCA/issues/447>`__)
-  Split D1 plots (MCMC convergence diagnostics plots & values)
   (`389 <https://github.com/asteca/ASteCA/issues/389>`__)
-  Explore Zeus as a possible addition to the best fit process
   (`457 <https://github.com/asteca/ASteCA/issues/457>`__)
-  Add mode, median to King’s profile plot
   (`470 <https://github.com/asteca/ASteCA/issues/470>`__)
-  Make “trim frame” option per cluster
   (`474 <https://github.com/asteca/ASteCA/issues/474>`__)
-  Closed due to old or not applicable
   (`209 <https://github.com/asteca/ASteCA/issues/209>`__,
   `293 <https://github.com/asteca/ASteCA/issues/293>`__,
   `399 <https://github.com/asteca/ASteCA/issues/399>`__)


`[v0.3.0] <https://github.com/asteca/asteca/releases/tag/v0.3.0>`__ - 2020-04-22
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Massive changes introduced in this new version. Python 2.7.x is no
longer supported.

-  Port to Python 3
   (`243 <https://github.com/asteca/ASteCA/issues/243>`__)
-  Upgrade to ``emcee`` v3.0.2
   (`423 <https://github.com/asteca/ASteCA/issues/423>`__)
-  Add ``emcee`` to the best fit process
   (`193 <https://github.com/asteca/ASteCA/issues/193>`__)
-  Upgraded to ``astropy`` v0.0.4
-  Remove (z,a) steps
   (`413 <https://github.com/asteca/ASteCA/issues/413>`__)
-  Bug fix: binary probabilities should not be averaged by
   ``zaWAverage``
   (`462 <https://github.com/asteca/ASteCA/issues/462>`__)
-  Add Tremmel’s implementation of the PLR
   (`447 <https://github.com/asteca/ASteCA/issues/447>`__)
-  Improve performance of synthetic cluster generation
   (`445 <https://github.com/asteca/ASteCA/issues/445>`__)
-  Fix Tolstoy likelihood accounting for uncertainties twice
   (`406 <https://github.com/asteca/ASteCA/issues/406>`__)
-  Add option to apply ’pmRA*cos(DE)’ correction
   (`452 <https://github.com/asteca/ASteCA/issues/452>`__)
-  Added ``optm`` method to local removal of stars
   (`432 <https://github.com/asteca/ASteCA/issues/432>`__)
-  Added ``manual`` binning method to likelihood block
   (`325 <https://github.com/asteca/ASteCA/issues/325>`__)
-  New radius estimating method and many improvements to structural
   functions (RDP, field dens, radius)
   (`454 <https://github.com/asteca/ASteCA/issues/454>`__,
   `449 <https://github.com/asteca/ASteCA/issues/449>`__,
   `346 <https://github.com/asteca/ASteCA/issues/346>`__,
   `378 <https://github.com/asteca/ASteCA/issues/378>`__)
-  Added maximum likelihood method for fitting King profiles
   (`268 <https://github.com/asteca/ASteCA/issues/268>`__,
   `298 <https://github.com/asteca/ASteCA/issues/298>`__)
-  Allow seeding the synthetic cluster generation process
   (`196 <https://github.com/asteca/ASteCA/issues/196>`__)
-  Add stopping condition to the plotting line
   (`443 <https://github.com/asteca/ASteCA/issues/443>`__)
-  Add Nsigma region to the best fit synthetic cluster
   (`460 <https://github.com/asteca/ASteCA/issues/460>`__)
-  Fix small bug in radii arrows
   (`182 <https://github.com/asteca/ASteCA/issues/182>`__)


`[v0.2.7] <https://github.com/asteca/asteca/releases/tag/v0.2.7>`__ - 2019-10-03
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Use inverse transform sampling to sample the IMF
   (`434 <https://github.com/asteca/ASteCA/issues/434>`__)
-  Interpolation of (z,a) values uses wrong m_ini index
   (`440 <https://github.com/asteca/ASteCA/issues/439>`__)
-  Interpolation of isochrone fails when (z,a) are both fixed
   (`439 <https://github.com/asteca/ASteCA/issues/440>`__)
-  Mass ‘alignment’ in zaInterp() gives poor result
   (`441 <https://github.com/asteca/ASteCA/issues/441>`__)
-  Select the N_mass_interp number automatically
   (`438 <https://github.com/asteca/ASteCA/issues/438>`__)


`[v0.2.6] <https://github.com/asteca/asteca/releases/tag/v0.2.6>`__ - 2019-09-19
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Fix normalization in Bayesian DA
   (`426 <https://github.com/asteca/ASteCA/issues/426>`__)
-  Fix function to detect X11 that fails in Mac OS (Windows too?)
   (`428 <https://github.com/asteca/ASteCA/issues/428>`__)
-  Merge ``semi_input.dat`` file into ``params_input.dat`` and copy
   input file as output
   (`427 <https://github.com/asteca/ASteCA/issues/427>`__)
-  Remove modes (`429 <https://github.com/asteca/ASteCA/issues/429>`__)
-  Use one photometric systems file instead of two identical ones
   (`421 <https://github.com/asteca/ASteCA/issues/421>`__)
-  Fix Ext/Imm operator causing spurious points in the GA
   (`424 <https://github.com/asteca/ASteCA/issues/424>`__)


`[v0.2.5] <https://github.com/asteca/asteca/releases/tag/v0.2.5>`__ - 2019-08-07
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Added the ``ptemcee`` method, and deprecated (for now) the BF
   (`367 <https://github.com/asteca/ASteCA/issues/367>`__)
-  Accept a CMD/CCD from mixed photometric systems
   (`228 <https://github.com/asteca/ASteCA/issues/228>`__,
   `229 <https://github.com/asteca/ASteCA/issues/229>`__)
-  Add support for the new set of isochrones PARSEC+COLIBRI
   (`322 <https://github.com/asteca/ASteCA/issues/322>`__)
-  Output all information obtained from the bootstrap
   (`279 <https://github.com/asteca/ASteCA/issues/279>`__)
-  Mask stars with photometry outside of reasonable range
   (`414 <https://github.com/asteca/ASteCA/issues/414>`__)
-  Add proper motions, parallax, and radial velocity support to Bayesian
   DA (`220 <https://github.com/asteca/ASteCA/issues/220>`__)
-  Use stars with no complete data in the Bayesian equation
   (`377 <https://github.com/asteca/ASteCA/issues/377>`__).
-  Add dimensional `weights to Bayesian
   DA <https://github.com/asteca/ASteCA/commit/d8a2ba99f6d36cbfb9e09efe08e1f590eb156743>`__.
-  Use all positions for structural functions
   (`107 <https://github.com/asteca/ASteCA/issues/107>`__).
-  Make the bootstrap the actual method (instead of GA)
   (`64 <https://github.com/asteca/ASteCA/issues/64>`__)
-  Make the GA work with floats instead of a grid
   (`412 <https://github.com/asteca/ASteCA/issues/412>`__)
-  Plot the incomplete dataset with MPs information
   (`411 <https://github.com/asteca/ASteCA/issues/411>`__)
-  Use a total number of masses, not a step value
   (`410 <https://github.com/asteca/ASteCA/issues/410>`__)
-  Use stars after error rejection for LF & completeness
   (`390 <https://github.com/asteca/ASteCA/issues/390>`__)
-  Switch to astropy’s read module
   (`327 <https://github.com/asteca/ASteCA/issues/327>`__) and allow
   `reading columns by
   name <https://github.com/asteca/ASteCA/commit/08d2c04ab5a5307aba3d19762bbb7f64df4f1aae>`__.
-  Update check for `installed
   packages <https://github.com/asteca/ASteCA/commit/bb885f9cc9acc311d57e312ac6c4623ec7ff235b>`__
   (newer ``pip`` threw an error).
-  Added a 2D cluster vs field KDE comparison, and the A-D test
   (`255 <https://github.com/asteca/ASteCA/issues/255>`__,
   `356 <https://github.com/asteca/ASteCA/issues/356>`__)
-  Added MAP, median and mode to output parameters.
-  Added R2 normality estimator to distributions
   (`401 <https://github.com/asteca/ASteCA/issues/401>`__)
-  Deprecated `KDE p-value
   function <https://github.com/asteca/ASteCA/commit/f218148e1f2a7abff591816c2271a7c6e2dc61ac>`__.
-  Deprecated ``trim_frame``, and ``manual`` `mode in photometric error
   rejection <https://github.com/asteca/ASteCA/commit/783975b22b8773c4ab08b3f1588e616cd3c858b2>`__.
-  Deprecated `integrated magnitude
   function <https://github.com/asteca/ASteCA/commit/1130c905e82048053267d3fcba41a967a88f77a2>`__.
-  Store input parameters as .json for each cluster
   (`126 <https://github.com/asteca/ASteCA/issues/126>`__)
-  Don’t read hidden files from the ‘isochrones’ folder
   (`403 <https://github.com/asteca/ASteCA/issues/403>`__)
-  Use KDE instead of Gaussian filters
   (`379 <https://github.com/asteca/ASteCA/issues/379>`__)
-  Split C2 plot into C2 and C3


`[v0.2.4] <https://github.com/asteca/asteca/releases/tag/v0.2.4>`__ - 2018-03-16
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Extend support for up to two colors.
-  Improved performance
   (`#357 <https://github.com/asteca/ASteCA/issues/357>`__):

   -  Make mass sampling optional
      (`#373 <https://github.com/asteca/ASteCA/issues/373>`__)
   -  Move binarity assignment outside of the synthetic cluster
      generation.
   -  Move isochrone sorting outside of the synthetic cluster
      generation.
   -  Move random floats for photometric errors outside of the synthetic
      cluster generation.
   -  Move random floats for completeness outside of the synthetic
      cluster generation. Code is now ~3.3X faster


`[v0.2.3] <https://github.com/asteca/asteca/releases/tag/v0.2.3>`__ - 2017-09-23
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Improved performance of synthetic cluster generation
   (`#227 <https://github.com/asteca/ASteCA/issues/227>`__). Code is now
   ~4X faster.
-  Fix excessive use of memory by Rbf interpolation
   (`#350 <https://github.com/asteca/ASteCA/issues/350>`__)
-  Use equal bin widths in LF and completeness function
   (`#300 <https://github.com/asteca/ASteCA/issues/300>`__)
-  Faster star separation by errors
   (`#351 <https://github.com/asteca/ASteCA/issues/351>`__)
-  Generalize Bayesian DA to N-dimensions, fix statistical issues,
   improve performance
   (`#352 <https://github.com/asteca/ASteCA/issues/352>`__)


`[v0.2.2] <https://github.com/asteca/asteca/releases/tag/v0.2.2>`__ - 2017-08-29
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Add weights to binned likelihood
   (`#216 <https://github.com/asteca/ASteCA/issues/216>`__)
-  Fix `bug in progress
   bar <https://github.com/asteca/ASteCA/commit/65d1f89bd0992120c8401c80ef976ba3c3803c38>`__.
-  Identify binaries in `plotted HR
   diagram <https://github.com/asteca/ASteCA/commit/7c650fb9b65090ea54064d385aa28087b3008c80>`__.
-  Modify the information presented by the `2-parameters density
   plots <https://github.com/asteca/ASteCA/commit/ec38070b4bb2c6d48d50c2bbd265f15bcc6347ee>`__.
   Takes care of `#71 <https://github.com/asteca/ASteCA/issues/71>`__.
-  Smarter empty field region around cluster region
   (`#345 <https://github.com/asteca/ASteCA/issues/345>`__).
-  Detect stars with duplicate IDs in data file
   (`#212 <https://github.com/asteca/ASteCA/issues/212>`__).


`[v0.2.1] <https://github.com/asteca/asteca/releases/tag/v0.2.1>`__ - 2017-08-11
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Fix issue with ‘tolstoy’ likelihood estimation
   (`#340 <https://github.com/asteca/ASteCA/issues/340>`__)
-  Fix a couple of issues with the error curve fitting
   (`#338 <https://github.com/asteca/ASteCA/issues/338>`__)
-  Add ‘fixed’ MPs algorithm (useful when no field region is available)
   (`#326 <https://github.com/asteca/ASteCA/issues/326>`__)
-  Fix crash when obtaining error curve
   (`#256 <https://github.com/asteca/ASteCA/issues/256>`__)


`[v0.2.0] <https://github.com/asteca/asteca/releases/tag/v0.2.0>`__ - 2017-08-07
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Generalized code to accept an arbitrary CMD in any *single*
   photometric system supported by the `CMD
   service <http://stev.oapd.inaf.it/cgi-bin/cmd>`__
   (`#24 <https://github.com/asteca/ASteCA/issues/24>`__).
-  Identify binary systems in synthetic clusters
   (`#199 <https://github.com/asteca/ASteCA/issues/199>`__).
-  Plots are now produced per blocks, instead of all together at the end
   (`#271 <https://github.com/asteca/ASteCA/issues/271>`__)
-  Switch dependency requirement from astroML to astropy
   (`#303 <https://github.com/asteca/ASteCA/issues/303>`__).
-  Remove unused error rejection modes
   (`#331 <https://github.com/asteca/ASteCA/issues/331>`__)
-  Simplify params_input.dat file
   (`#217 <https://github.com/asteca/ASteCA/issues/217>`__)
-  Check that all metallicity files contain the same number of age
   values (`#218 <https://github.com/asteca/ASteCA/issues/218>`__)
-  Add density maps analysis for center function
   (`#164 <https://github.com/asteca/ASteCA/issues/164>`__)
-  Remove weight added to the observed cluster CMD’s histogram
   (`#308 <https://github.com/asteca/ASteCA/issues/308>`__)
-  Fix bad parameter rounding
   (`#248 <https://github.com/asteca/ASteCA/issues/248>`__)
-  Add ‘max mag’ cut for synthetic clusters
   (`#302 <https://github.com/asteca/ASteCA/issues/302>`__,
   `#264 <https://github.com/asteca/ASteCA/issues/264>`__)
-  Simplify installation steps
   (`#88 <https://github.com/asteca/ASteCA/issues/88>`__,
   `#315 <https://github.com/asteca/ASteCA/issues/315>`__)
-  Plot results of brute force minimization
   (`#100 <https://github.com/asteca/ASteCA/issues/100>`__)
-  Make extinction parameter Rv a manual input parameter
   (`#314 <https://github.com/asteca/ASteCA/issues/314>`__)
-  Use numpy’s binning methods
   (`#317 <https://github.com/asteca/ASteCA/issues/317>`__)
-  Modify RDP limit
   (`#294 <https://github.com/asteca/ASteCA/issues/294>`__)
-  Store extra data from theoretical isochrones
   (`#201 <https://github.com/asteca/ASteCA/issues/201>`__)


`[v0.1.9.5] <https://github.com/asteca/asteca/releases/tag/v0.1.9.5>`__ - 2016-08-07
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Remove forgotten print line.
-  Print relevant information when data con not be read
   (`#262 <https://github.com/asteca/asteca/issues/262>`__).
-  Fix bad range issue
   (`#226 <https://github.com/asteca/asteca/issues/226>`__).


`[v0.1.9.4] <https://github.com/asteca/asteca/releases/tag/v0.1.9.4>`__ - 2016-07-25
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Add support for five tracks from the CMD service
   (`#276 <https://github.com/asteca/ASteCA/issues/276>`__).
-  Read metallicity files with underscores instead of decimal dots
   (`#277 <https://github.com/asteca/ASteCA/issues/277>`__).
-  Several important structural changes
   (`#273 <https://github.com/asteca/asteca/issues/273>`__): add
   ``first_run`` check, re-arrange and re-name modules, and move almost
   every part of the code into the ``packages/`` folder.


`[v0.1.9.3] <https://github.com/asteca/asteca/releases/tag/v0.1.9.3>`__ - 2016-05-25
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Add support for CMD in the `HST/ACS WFC photometric
   system <http://www.stsci.edu/hst/acs>`__ (requested by Daniel
   Arbelaez).


`[v0.1.9.2] <https://github.com/asteca/asteca/releases/tag/v0.1.9.2>`__ - 2016-04-17
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Add support for three CMDs in the `Strömgren photometric
   system <https://en.wikipedia.org/wiki/Str%C3%B6mgren_photometric_system>`__
   (requested by J. Hughes Clark).
-  Change likelihood density plots to `scatter
   plots <https://github.com/asteca/ASteCA/commit/6bac8749ba9b6b8c0fbaa2b226cca272e110e1cf>`__
   which show more information.
-  Add extra condition for DA break: minimum 10% of the runs `must have
   passed <https://github.com/asteca/ASteCA/commit/7095c0cd043804cce25d27a9e16650ecf8a2f7a5>`__.
-  Fix bug with `‘mag’
   mode <https://github.com/asteca/ASteCA/commit/272ed205d4beaaa8d3a10b2c664550140e238053>`__
   in ‘Reduced membership’, wouldn’t run if the Bayesian DA was skipped.
-  Fix minor bug
   (`#241 <https://github.com/asteca/asteca/issues/241>`__) when
   `printing KP results to
   screen <https://github.com/asteca/ASteCA/commit/62ffe4dad93fd5291900c08aa05af9e1c1cee5f2>`__.


`[v0.1.9.1] <https://github.com/asteca/asteca/releases/tag/v0.1.9.1>`__ - 2015-08-25
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Fixed rounding of errors that returned 0. values if error was larger
   than value (`#213 <https://github.com/asteca/asteca/issues/213>`__).
-  Check if ``pip`` module is installed + search for installed packages
   `globally, not
   locally <https://github.com/asteca/ASteCA/commit/3d04bb5247e001cf033a3df47e9f89e21c9dd2e5>`__.
-  Catch `badly
   formatted <https://github.com/asteca/ASteCA/commit/11ed705d9b23730ef8752d4553139c45700c0074>`__
   input data file.
-  Restructure `King radii
   obtention <https://github.com/asteca/ASteCA/commit/4d201b76edace038d6651b7c43ac997728de1c82>`__.
-  `Correctly plot
   stars <https://github.com/asteca/ASteCA/commit/c3ccc376a5d46415ae45b9f2e4572be50b75847d>`__
   in cluster region, not used in best fit function.


`[v0.1.9] <https://github.com/asteca/asteca/releases/tag/v0.1.9>`__ - 2015-06-18
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(**Warning**: this release breaks compatibility with the previous
version of the ``params_input.dat`` & ``semi_input.dat`` files)

-  Models (ie: isochrone + extinction +distance modulus + mass
   distribution + binarity) are now evaluated *each time the GA selects
   them as a solution*, thus a new mass distribution is generated
   (`#186 <https://github.com/asteca/asteca/issues/186>`__). This has a
   performance cost, but provides higher accuracy in the best model
   assignment process since a single model can now be evaluated with a
   slightly different mass distribution several times (only with GA,
   *Brute Force* method will only process a model once).
-  Added an *exit switch* to the decontamination algorithm. It will stop
   iterations if the MPs converged to 0.1% tolerance values for all the
   stars in the cluster region (compared to the previous iteration).
   This speeds up the function considerably
   (`#185 <https://github.com/asteca/asteca/issues/185>`__).
-  The upper mass value in the IMF can now be `modified via the input
   parameters
   file <https://github.com/asteca/asteca/commit/4b1a897d69cf85b1c0263d738cf2132d9924eb9c>`__.
-  Code can now read ``params_input_XX.dat`` files when `using lazy
   parallelization <https://github.com/asteca/asteca/commit/f2508355d8136c2d5a6216093e6f9eda02bd99c1>`__.
-  Number of field regions `can now be set
   individually <https://github.com/asteca/ASteCA/commit/dc4c9223b0ec0a02904e30025eec50dfdc13637d>`__
   via the ``semi_input.dat`` file.
-  `Added ‘bb’ binning
   method <https://github.com/asteca/ASteCA/commit/d35c5611708d249e730bef77b0ee14226cce14de>`__
   based on `Bonnato & Bica
   (2007) <http://adsabs.harvard.edu/abs/2007MNRAS.377.1301B>`__. Sets
   bin widths of 0.25 and 0.5 for colors and magnitudes respectively.
-  Fixed bug in ``manual`` mode when `displaying
   errors <https://github.com/asteca/asteca/commit/2e4b1d8f8a084e78bc56d52df494a796a6909de6>`__.
-  Fixed bug when narrow frames were plotted
   (`#168 <https://github.com/asteca/asteca/issues/168>`__).
-  Moved text box outside of synthetic cluster’s plot to improve its
   visibility (`#205 <https://github.com/asteca/asteca/issues/205>`__).
-  Closed `#13 <https://github.com/asteca/asteca/issues/13>`__. Saha’s W
   likelihood needs the number of model points to be fixed, which
   prevents it from being used when the mass varies. There’s nothing to
   be gained by adding this function.
-  Caveat dragged from version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__ is
   `resolved <https://github.com/asteca/ASteCA/commit/ff3b240ec3d1b2339ce51cf262e71810a33b6517>`__.


`[v0.1.8] <https://github.com/asteca/asteca/releases/tag/v0.1.8>`__ - 2015-04-09
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(**Warning**: this release breaks compatibility with the previous
version of the ``params_input.dat`` file)

-  Added ``local`` and ``mp_05`` methods to the selection of which stars
   to use in the best fit cluster parameter assignation process
   (`#180 <https://github.com/asteca/asteca/issues/180>`__,
   `#183 <https://github.com/asteca/asteca/issues/183>`__).
-  Added an *automatic update checker* function that notifies the user
   if an updated version of ``ASteCA`` is available for download
   (`#179 <https://github.com/asteca/asteca/issues/179>`__).
-  Added grid lines over the photometric diagrams of the observed and
   synthetic cluster, showing the binning made by the method selected in
   each case (`#131 <https://github.com/asteca/asteca/issues/131>`__).
-  Best fit synthetic cluster found is now saved to file
   (`#154 <https://github.com/asteca/asteca/issues/154>`__).
-  Correctly obtain approximate number of members (``n_memb``) and
   contamination index (``CI``) when the cluster radius extends beyond
   the RDP, thus making the field star density value (``field_dens``)
   unreliable (`#111 <https://github.com/asteca/asteca/issues/111>`__).
-  Added ``f10`` flag to alert when the ``memb_par`` value is greater
   than +-0.33, which means that there are twice as many estimated true
   members in either method
   (`#175 <https://github.com/asteca/asteca/issues/175>`__).
-  Improved ``top_tiers`` plotting and saved file
   (`#184 <https://github.com/asteca/asteca/issues/184>`__).

**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.7] <https://github.com/asteca/asteca/releases/tag/v0.1.7>`__ - 2015-03-26
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(**Warning**: this release breaks compatibility with the previous
version of the ``params_input.dat`` file)

-  Re-write ``lowexp`` `error rejection
   method <https://github.com/asteca/asteca/commit/6b2857aefa2878ee5aba245a7fbf9cc1f423820b>`__,
   now uses *prediction bands* instead of *confidence intervals*.
-  Force ``matplotlib``\ ’s backend to make the code `work in
   servers <https://github.com/asteca/asteca/commit/197af6439baabd3e9db4039775aba721d84047a2>`__.
-  Fixed ``eyefit`` method for `error
   rejection <https://github.com/asteca/asteca/commit/d92be0c8e398739fba562d59ba35b11eeac9a9a0>`__.
   It changed after fixing
   `#169 <https://github.com/asteca/asteca/issues/169>`__.
-  Added `SDSS
   CMDs <https://github.com/asteca/asteca/commit/2324a70f402ddbe9fdde203c3745f93b6d6dc545>`__
   ``g vs (u-g)`` & ``g vs (g-r)``, at the request of Tjibaria Pijloo
   (Department of Astrophysics, Radboud University Nijmegen).
-  Fixed bug in binarity generation for the CMDs of the form
   ``X vs (X-Y)``
   (`#181 <https://github.com/asteca/asteca/issues/181>`__).
-  Smarter selection of stars to be used by the best fit function,
   improvements in several plots
   (`#171 <https://github.com/asteca/asteca/issues/171>`__,
   `#172 <https://github.com/asteca/asteca/issues/172>`__).
-  Best fit function can now accept a *minimum magnitude* value instead
   of just a *minimum probability* value
   (`#115 <https://github.com/asteca/asteca/issues/115>`__).
-  Added a ``memb_par`` parameter to compare the number of approximate
   cluster members obtained via the structural analysis and via the
   decontamination algorithm
   (`#162 <https://github.com/asteca/asteca/issues/162>`__).
-  Code is now able to correctly read the names of files with `more than
   one dot in it’s
   name <https://github.com/asteca/asteca/commit/c0358ed9526b835bfeeddf75804002ad51c69610>`__.
-  Fixed bad `alphabetical
   ordering <https://github.com/asteca/asteca/commit/b6ca2a2df8b7e614dc9beb38e99400e3b69208bf>`__
   of input cluster files.
-  Better limits for photometric diagram
   (`#173 <https://github.com/asteca/asteca/issues/173>`__).
-  Fixed ``DeprecationWarning``
   `issue <https://github.com/asteca/asteca/commit/97d77f1d7f36adf6af6398a2f4a5b944598fda8f>`__.
-  Invert x axis when `RA cords are
   used <https://github.com/asteca/asteca/commit/e99da37a398c446d71c59c43f4547434d0c9f7e7>`__
   (improved
   `here <https://github.com/asteca/asteca/commit/aeb7d7d097eb40289d2bb4c83adf433567bb28d0>`__).
-  Several fixes and improvements made to plotted diagrams
   (`5c7dc7f <https://github.com/asteca/asteca/commit/5c7dc7f9f348bf2bedb3eb86daf7decbbf83df33>`__;
   `1642349 <https://github.com/asteca/asteca/commit/16423496d22bb843294189fd121a0ed8a0c6e783>`__;
   `b57028c <https://github.com/asteca/asteca/commit/b57028c93259afbf3cbebc905c482349fcb6ef7a>`__;
   `240178a <https://github.com/asteca/asteca/commit/240178a3c797910d6a807a41a8dd6c2f94d82cfb>`__;
   `9ec0ab8 <https://github.com/asteca/asteca/commit/9ec0ab8c3d966e0dbe19c6b5cff65e1cb381c939>`__;
   `fef14c4 <https://github.com/asteca/asteca/commit/fef14c476b88bc9f82bcd39e96cee222a0628cdd>`__;
   `db0df2a <https://github.com/asteca/asteca/commit/db0df2adc8d9821ab5122ba6b6482557627a779e>`__;
   `575ebe7 <https://github.com/asteca/asteca/commit/575ebe7de64c1c4da04eb7c18dfab4b8bd1b2751>`__;
   `#177 <https://github.com/asteca/asteca/issues/177>`__;
   `#178 <https://github.com/asteca/asteca/issues/178>`__).


**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.61] <https://github.com/asteca/asteca/releases/tag/v0.1.61>`__ - 2015-03-04
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Added `“lazy
   parallelization” <https://github.com/asteca/asteca/commit/b536c84c2ad085bbe8ff10a0b6535618ae1ba09a>`__
   ability. Now the user can run as many instances of the code as needed
   simply by creating extra ``asteca_xx.py`` and ``input_xx`` folders
   where ``xx`` are integers of the form: 01, 02,…, 99.
-  `Reposition <https://github.com/asteca/asteca/commit/e7dec4b75a62ff397ee62cb322345f6b17b74ff6>`__
   several text boxes in output images, newer versions of ``matplotlib``
   moved them from the previous position.
-  Fix `bad
   import <https://github.com/asteca/asteca/commit/9bed2166e9cc36faa7077c79c436c50e40801820>`__
   of ``rpy2`` package, positioned incorrectly in two functions.
-  Fix ``DeprecationWarning`` showing when ``exp_function`` was used
   (`#169 <https://github.com/asteca/asteca/issues/169>`__).


**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.5] <https://github.com/asteca/asteca/releases/tag/v0.1.5>`__ - 2015-03-03
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(**Warning**: this release breaks compatibility with the previous
version of the ``params_input.dat`` file)

-  Improved radius assignment algorithm
   (`#146 <https://github.com/asteca/asteca/issues/146>`__).
-  Detect cropped cluster region and use correct area when generating
   field regions
   (`#139 <https://github.com/asteca/asteca/issues/139>`__,
   `#157 <https://github.com/asteca/asteca/issues/157>`__).
-  Fixed bug that crashed the code when top tiers synthetic clusters
   with no stars were plotted
   (`#147 <https://github.com/asteca/asteca/issues/147>`__). Added
   minimum total mass of 10Mo.
-  Fixed bug where KDE p-values for field vs field comparison were
   artificially increased by comparing a field region with itself
   (`#138 <https://github.com/asteca/asteca/issues/138>`__).
-  Obtain KDE p-value even if one field region is defined
   (`#114 <https://github.com/asteca/asteca/issues/114>`__).
-  Fixed small bug that prevented integrated magnitude curves from being
   plotted (`#145 <https://github.com/asteca/asteca/issues/145>`__).
-  Fixed several smaller bugs and issues
   (`#110 <https://github.com/asteca/asteca/issues/110>`__,
   `#150 <https://github.com/asteca/asteca/issues/150>`__,
   `#140 <https://github.com/asteca/asteca/issues/140>`__,
   `#142 <https://github.com/asteca/asteca/issues/142>`__,
   `#141 <https://github.com/asteca/asteca/issues/141>`__,
   `#149 <https://github.com/asteca/asteca/issues/149>`__,
   `#95 <https://github.com/asteca/asteca/issues/95>`__,
   `#148 <https://github.com/asteca/asteca/issues/148>`__,
   `#136 <https://github.com/asteca/asteca/issues/136>`__,
   `#163 <https://github.com/asteca/asteca/issues/163>`__,
   `#143 <https://github.com/asteca/asteca/issues/143>`__).


**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.4] <https://github.com/asteca/asteca/releases/tag/v0.1.4>`__ - 2014-12-18
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Improved plotting of crowded fields
   (`#62 <https://github.com/asteca/asteca/issues/62>`__).
-  Function to generate image is now more stable
   (`#112 <https://github.com/asteca/asteca/issues/112>`__). Re-arranged
   plots in output image.
-  Add *Top tiers* models output
   (`#130 <https://github.com/asteca/asteca/issues/130>`__).
-  Fixed small bug in KDE p-values function
   (`#134 <https://github.com/asteca/asteca/issues/134>`__).
-  Minor re-arrangement with semi-input data.


**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.3] <https://github.com/asteca/asteca/releases/tag/v0.1.3>`__ - 2014-12-10
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Accept arrays of non-equispaced parameter values instead of only
   equispaced ranges
   (`#121 <https://github.com/asteca/asteca/issues/121>`__).
-  Added support for log-normal `Chabrier
   (2001) <http://adsabs.harvard.edu/abs/2001ApJ...554.1274C>`__ IMF.
-  More precise encoding/decoding in genetic algorithm.
-  Functions separated into sections
   (`#125 <https://github.com/asteca/asteca/issues/125>`__).
-  Input parameters set as global variables
   (`#132 <https://github.com/asteca/asteca/issues/132>`__).


**Caveats**

-  Same as version
   `0.1.2 <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__.


`[v0.1.2] <https://github.com/asteca/asteca/releases/tag/v0.1.2>`__ - 2014-12-01
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Likelihood method now supports `Dolphin
   (2002) <http://adsabs.harvard.edu/abs/2002MNRAS.332...91D>`__
   *Poisson likelihood ratio* function.
-  Closed `#120 <https://github.com/asteca/asteca/issues/120>`__,
   `#101 <https://github.com/asteca/asteca/issues/101>`__,
   `#129 <https://github.com/asteca/asteca/issues/129>`__,
   `#124 <https://github.com/asteca/asteca/issues/124>`__,
   `#102 <https://github.com/asteca/asteca/issues/102>`__.
-  Minor `position
   fix <https://github.com/asteca/asteca/commit/00538bda879009bae0a4e7565b124c8939c75d0f>`__
   for synthetic cluster text box in output plot.
-  Brute force algorithm now returns `correct
   errors <https://github.com/asteca/asteca/commit/afe30cbdff561a90986a638c55a4b7247fd0bc53>`__.
-  Some fixes for when unique values in the input parameter ranges are
   used
   (`[1] <https://github.com/asteca/asteca/commit/7cc383d799f2af5c1f1f8a6dcfc80e639461f02d>`__,
   `[2] <https://github.com/asteca/asteca/commit/c6505025d4c3b6147a2913fad648dc18c125376b>`__).
-  Replaced deprecated `compiler
   package <https://github.com/asteca/asteca/commit/f9e8c5edba5f5ca8cc33ec1afb4d137f7167e8df>`__
   used to flatten list.


**Caveats**

-  Still not sure why *tolstoy* likelihood is biased towards high masses
   :confused:


`[v0.1.1] <https://github.com/asteca/asteca/releases/tag/v0.1.1>`__ - 2014-11-07
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*More stable release.*

-  Closed `#113 <https://github.com/asteca/asteca/issues/113>`__,
   `#116 <https://github.com/asteca/asteca/issues/116>`__.
-  Minor
   `change <https://github.com/asteca/asteca/commit/3cffb4faa0c1dc6956aae2217c73afb4f392e53d>`__
   to error function.
-  Closed *Known issues* from previous version.


**Caveats**

-  Same as previous version.


`[v0.1.0] <https://github.com/asteca/asteca/releases/tag/v0.1.0>`__ - 2014-10-08
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*First semi-stable buggy release*

-  Closed `#72 <https://github.com/asteca/asteca/issues/72>`__,
   `#99 <https://github.com/asteca/asteca/issues/99>`__,
   `#37 <https://github.com/asteca/asteca/issues/37>`__.
-  Changed the way the IMF was
   `sampled <https://github.com/Gabriel-p/asteca/commit/0671e74c52fbecde6bcbb1afb1c2624875156e57>`__,
   now it should be faster and more precise.
-  Some speed improvements (moved things around mainly).
-  Binary fraction is now a free parameter.

**Known issues**

-  **Serious bug**: if the DA is set to run but the *Best fit method*
   isn’t, the final plot can’t be produced since the ``syn_cl_err``
   function isn’t used
   (`fixed <https://github.com/Gabriel-p/asteca/commit/3e806bd0af5d7fcd7c8f2940716df880f4c1b67d>`__
   in next release).
-  Forgotten ``print`` prints out mass values every time the E/I
   operator is applied
   (`fixed <https://github.com/Gabriel-p/asteca/commit/8b313ef60fddccc41fd6fb7b9746f75f3e867d39>`__
   in next release).
-  If the number of points (``n_left``) in the radius finding function
   is smaller than 4, a very small radius is likely to be selected.
   `Fixed <https://github.com/Gabriel-p/asteca/commit/c247fd7fa4cca4d6bb341263434a4a43a4778efd>`__
   in next release.


**Caveats**

-  The total initial mass can be set as a free parameter but the
   likelihood function will select always synthetic clusters of high
   mass. Thus it is advised to leave this parameter fixed to 1000 solar
   masses.
-  The binary fraction found is not stored in the output data file.
-  Some density map plots for mass and binary fraction are missing.



`[v4.0.0-beta] <https://github.com/asteca/asteca/releases/tag/v4.0.0-beta>`__ - 2014-09-23
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Closed `#85 <https://github.com/asteca/asteca/issues/85>`__,
   `#70 <https://github.com/asteca/asteca/issues/70>`__,
   `#43 <https://github.com/asteca/asteca/issues/43>`__,
   `#86 <https://github.com/asteca/asteca/issues/86>`__.
-  Metallicity and age now take steps in the GA.
-  Add
   `checker <https://github.com/Gabriel-p/asteca/blob/master/functions/checker.py>`__
   function to make sure certain parameters are set correctly before
   running.
-  Number of points in ``get_radius`` increased 20% –> 25% of `the
   RDP <https://github.com/Gabriel-p/asteca/commit/a2e9b8f16111d5adafe66fed1eb64ed8bc03997b>`__.



`[v3.0.0-beta] <https://github.com/asteca/asteca/releases/tag/v3.0.0-beta>`__ - 2014-09-16
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Closed: `#89 <https://github.com/asteca/asteca/issues/89>`__,
   `#77 <https://github.com/asteca/asteca/issues/77>`__,
   `#80 <https://github.com/asteca/asteca/issues/80>`__.
-  The ``params_input.dat`` and ``semi_input.dat`` files are now located
   at the top level next to ``asteca.py``.
-  Cluster’s photometric files are not longer required to be stored
   inside a sub-folder to be picked-up by the code.



`[v2.0.1-beta] <https://github.com/asteca/asteca/releases/tag/v2.0.1-beta>`__ - 2014-09-15
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Correct version number.



`[v2.0.0-beta] <https://github.com/asteca/asteca/releases/tag/v2.0.0-beta>`__ - 2014-09-11
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-  Closed issues: `#15 <https://github.com/asteca/asteca/issues/15>`__,
   `#73 <https://github.com/asteca/asteca/issues/73>`__,
   `#53 <https://github.com/asteca/asteca/issues/53>`__,
   `#24 <https://github.com/asteca/asteca/issues/24>`__,
   `#75 <https://github.com/asteca/asteca/issues/75>`__,
   `#79 <https://github.com/asteca/asteca/issues/79>`__,
   `#81 <https://github.com/asteca/asteca/issues/81>`__,
   `#59 <https://github.com/asteca/asteca/issues/59>`__,
   `#83 <https://github.com/asteca/asteca/issues/83>`__,
   `#78 <https://github.com/asteca/asteca/issues/78>`__,
   `#69 <https://github.com/asteca/asteca/issues/69>`__,
   `#74 <https://github.com/asteca/asteca/issues/74>`__.
-  Changed name of package (OCAAT –> ASteCA).
-  Added separate function to handle the spatial 2D histogram.
-  Changes to ``get_center`` function (now hopefully simpler)
-  Added UBVI support for *V vs (U-V)*.
-  Added 2MASS CMD support for *J vs (J-H)*, *H vs (J-H)* and *K vs
   (H-K)*.
-  Improve field star regions integrated magnitudes curve averaging.
-  Simplify process of adding a new CMD.
-  Added details on how the integrated magnitude calculation is done in
   the manual.
-  Lots of minor edits/corrections.



`[v1.0.0-beta] <https://github.com/asteca/asteca/releases/tag/v1.0.0-beta>`__ - 2014-08-24
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*First beta release*

Version used (with some small changes) in the `original
article <http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html>`__.
