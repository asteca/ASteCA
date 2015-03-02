# Change Log

## [[0.1.5]](https://github.com/asteca/asteca/releases/tag/v0.1.5) - 2015-03-02

### Changed

* Improved radius assignment algorithm (#146).
* Detect cropped cluster region and use correct area when generating field regions (#139, #157).
* Fixed bug that crashed the code when top tiers synthetic clusters with no stars were plotted (#147). Added minimum total mass of 10Mo.
* Fixed bug where KDE p-values for field vs field comparision was artificially increased by comparing a field region with itself (#138).
* Obtain KDE p-value even if one field region is defined (#114).
* Fixed small bug that prevented integrated magnitude curves from being plotted (#145).
* Fixed several smaller bugs and issues (#110, #150, #140, #142, #141, #149, #95, #148, #136).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[0.1.4]](https://github.com/asteca/asteca/releases/tag/v0.1.4) - 2014-12-18

### Changed

* Improved plotting of crowded fields (#62).
* Function to generate image is now more stable (#112). Re-arranged plots in output image.
* Add _Top tiers_ models output (#130).
* Fixed small bug in KDE p-values function (#134).
* Minor re-arrangement with semi-input data.

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[0.1.3]](https://github.com/asteca/asteca/releases/tag/v0.1.3) - 2014-12-10

### Changed:

* Accept arrays of non-equispaced parameter values instead of only equispaced ranges (#121).
* Added support for lognormal [Chabrier (2001)](http://adsabs.harvard.edu/abs/2001ApJ...554.1274C) IMF.
* More precise encoding/decoding in genetic algorithm.
* Functions separated into sections (#125).
* Input parameters set as global variables (#132).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[0.1.2]](https://github.com/asteca/asteca/releases/tag/v0.1.2) - 2014-12-01

### Changed

* Likelihood method now supports [Dolphin (2002)](http://adsabs.harvard.edu/abs/2002MNRAS.332...91D) _Poisson likelihood ratio_ function.
* Closed #120, #101, #129, #124, #102.
* Minor [position fix](https://github.com/asteca/asteca/commit/00538bda879009bae0a4e7565b124c8939c75d0f) for synthetic cluster text box in output plot.
* Brute force algorithm now returns [correct errors](https://github.com/asteca/asteca/commit/afe30cbdff561a90986a638c55a4b7247fd0bc53).
* Some fixes for when unique values in the input parameter ranges are used [[1]](https://github.com/asteca/asteca/commit/7cc383d799f2af5c1f1f8a6dcfc80e639461f02d) [[2]](https://github.com/asteca/asteca/commit/c6505025d4c3b6147a2913fad648dc18c125376b).
* Replaced deprecated [compiler package](https://github.com/asteca/asteca/commit/f9e8c5edba5f5ca8cc33ec1afb4d137f7167e8df) used to flatten list.

### Caveats

 * Still not sure why _tolstoy_ likelihood is biased towards high masses :confused:

## [[0.1.1]](https://github.com/asteca/asteca/releases/tag/v0.1.1) - 2014-11-07

_More stable release._

### Changed

* Closed #113, #116.
* Minor [change](https://github.com/asteca/asteca/commit/3cffb4faa0c1dc6956aae2217c73afb4f392e53d) to error function.
* Closed _Known issues_ from previous version.

### Caveats

 * Same as previous version.

## [[0.1.0]](https://github.com/asteca/asteca/releases/tag/v0.1.0) - 2014-10-08

_First <s>semi-stable</s> buggy release_

### Changed

* Closed #72, #99, #37.
* Changed the way the IMF was [sampled](https://github.com/Gabriel-p/asteca/commit/0671e74c52fbecde6bcbb1afb1c2624875156e57), now it should be faster and more precise.
* Some speed improvements (moved things around mainly).
* Binary fraction is now a free parameter.

### Known issues

 * **Serious bug**: if the DA is set to run but the _Best fit method_ isn't, the final plot can't be produced since the `syn_cl_err` function isn't used ([fixed](https://github.com/Gabriel-p/asteca/commit/3e806bd0af5d7fcd7c8f2940716df880f4c1b67d) in next release).
 * Forgotten `print` prints out mass values every time the E/I operator is applied ([fixed](https://github.com/Gabriel-p/asteca/commit/8b313ef60fddccc41fd6fb7b9746f75f3e867d39) in next release).
 * If the number of points (`n_left`) in the radius finding function is smaller than 4, a very small radius is likely
to be selected. [Fixed](https://github.com/Gabriel-p/asteca/commit/c247fd7fa4cca4d6bb341263434a4a43a4778efd) in next release.

### Caveats

 * The total initial mass can be set as a free parameter but the likelihood function will select always synthetic clusters of high mass. Thus it is advised to leave this parameter fixed to 1000 solar masses.
 * The binary fraction found is not stored in the output data file.
 * Some density map plots for mass and binary fraction are missing.

## [[v4.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v4.0.0-beta) - 2014-09-23

### Changed

* Closed #85, #70, #43, #86.
- Metallicity and age now take steps in the GA.
- Add [checker](https://github.com/Gabriel-p/asteca/blob/master/functions/checker.py) function to make sure certain parameters are set correctly before running.
* Number of points in `get_radius` increased 20% --> 25% of the RDP (https://github.com/Gabriel-p/asteca/commit/a2e9b8f16111d5adafe66fed1eb64ed8bc03997b)

## [[v3.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v3.0.0-beta) - 2014-09-16

### Changed

* Closed: #89, #77, #80.
* The `params_input.dat` and `semi_input.dat` files are now located at the top level next to `asteca.py`.
* Cluster's photometric files are not longer required to be stored inside a sub-folder to be picked-up by the code.

## [[v2.0.1-beta]](https://github.com/asteca/asteca/releases/tag/v2.0.1-beta) - 2014-09-15

### Changed

* Correct version number.

## [[v2.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v2.0.0-beta) - 2014-09-11

### Changed

* Closed issues: #15, #73, #53, #24,  #75, #79, #81, #59, #83, #78, #69, #74.
* Changed name of package (OCAAT --> ASteCA).
* Added separate function to handle the spatial 2D histogram.
* Changes to `get_center` function (now hopefully simpler)
* Added UBVI support for _V vs (U-V)_.
* Added 2MASS CMD support for _J vs (J-H)_, _H vs (J-H)_ and _K vs (H-K)_.
* Improve field star regions integrated magnitudes curve averaging.
* Simplify process of adding a new CMD.
* Added details on how the integrated magnitude calculation is done in the manual.
* Lots of minor edits/corrections.

## [[v1.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v1.0.0-beta) - 2014-08-24

_First beta release_

Version used (with some small changes) in the [original article](http://arxiv.org/abs/1412.2366).

