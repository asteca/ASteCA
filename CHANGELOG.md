# Change Log

## [[v0.1.8]](https://github.com/asteca/asteca/releases/tag/v0.1.8) - 2015-04-09

### Changed

(**Warning**: this release breaks compatibility with the previous version of the `params_input.dat` file)

* Added `local`  and `mp_05` methods to the selection of which stars to use in the best fit cluster parameter assignation process (#180, #183).
* Added an _automatic update checker_ function that notifies the user if an updated version of `ASteCA` is available for download (#179).
* Added grid lines over the photometric diagrams of the observed and synthetic cluster, showing the binning made by the method selected in each case (#131).
* Best fit synthetic cluster found is now saved to file (#154).
* Correctly obtain approximate number of members (`n_memb`) and contamination index (`CI`) when the cluster radius extends beyond the RDP, thus making the field star density value (`field_dens`) unreliable (#111).
* Added `f10` flag to alert when the `memb_par` value is greater than +-0.33, which means that there are twice as many estimated true members in either method (#175).
* Improved `top_tiers` plotting and saved file (#184).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.7]](https://github.com/asteca/asteca/releases/tag/v0.1.7) - 2015-03-26

### Changed

(**Warning**: this release breaks compatibility with the previous version of the `params_input.dat` file)

* Re-write `lowexp` [error rejection method](https://github.com/asteca/asteca/commit/6b2857aefa2878ee5aba245a7fbf9cc1f423820b), now uses _prediction bands_ instead of _confidence intervals_.
* Force `matplotlib`'s backend to make the code [work in servers](https://github.com/asteca/asteca/commit/197af6439baabd3e9db4039775aba721d84047a2) .
* Fixed `eyefit` method for [error rejection](https://github.com/asteca/asteca/commit/d92be0c8e398739fba562d59ba35b11eeac9a9a0). It changed after fixing [#169](https://github.com/asteca/asteca/issues/169).
* Added [SDSS CMDs](https://github.com/asteca/asteca/commit/2324a70f402ddbe9fdde203c3745f93b6d6dc545) `g vs (u-g)` & `g vs (g-r)`, at the request of Tjibaria Pijloo (Department of Astrophysics, Radboud University Nijmegen).
* Fixed bug in binarity generation for the CMDs of the form `X vs (X-Y)` ([#181](https://github.com/asteca/asteca/issues/181)).
* Smarter selection of stars to be used by the best fit function, improvements in several plots ([#171](https://github.com/asteca/asteca/issues/171), [#172](https://github.com/asteca/asteca/issues/172)).
* Best fit function can now accept a _minimum magnitude_ value instead of just a _minimum probability_ value ([#115](https://github.com/asteca/asteca/issues/115)).
* Added a `memb_par` parameter to compare the number of approximate cluster members obtained via the structural analysis and via the decontamination algorithm ([#162](https://github.com/asteca/asteca/issues/162)).
* Code is now able to correctly read the names of files with [more than one dot in it's name](https://github.com/asteca/asteca/commit/c0358ed9526b835bfeeddf75804002ad51c69610).
* Fixed bad [alphabetical ordering](https://github.com/asteca/asteca/commit/b6ca2a2df8b7e614dc9beb38e99400e3b69208bf) of input cluster files.
* Better limits for photometric diagram ([#173](https://github.com/asteca/asteca/issues/173)).
* Fixed `DeprecationWarning` [issue](https://github.com/asteca/asteca/commit/97d77f1d7f36adf6af6398a2f4a5b944598fda8f).
* Invert x axis when [RA cords are used](https://github.com/asteca/asteca/commit/e99da37a398c446d71c59c43f4547434d0c9f7e7) (improved [here](https://github.com/asteca/asteca/commit/aeb7d7d097eb40289d2bb4c83adf433567bb28d0)).
* Several fixes and improvements made to plotted diagrams ([5c7dc7f](https://github.com/asteca/asteca/commit/5c7dc7f9f348bf2bedb3eb86daf7decbbf83df33); [1642349](https://github.com/asteca/asteca/commit/16423496d22bb843294189fd121a0ed8a0c6e783); [b57028c](https://github.com/asteca/asteca/commit/b57028c93259afbf3cbebc905c482349fcb6ef7a); [240178a](https://github.com/asteca/asteca/commit/240178a3c797910d6a807a41a8dd6c2f94d82cfb); [9ec0ab8](https://github.com/asteca/asteca/commit/9ec0ab8c3d966e0dbe19c6b5cff65e1cb381c939); [fef14c4](https://github.com/asteca/asteca/commit/fef14c476b88bc9f82bcd39e96cee222a0628cdd); [db0df2a](https://github.com/asteca/asteca/commit/db0df2adc8d9821ab5122ba6b6482557627a779e); [575ebe7](https://github.com/asteca/asteca/commit/575ebe7de64c1c4da04eb7c18dfab4b8bd1b2751); [#177](https://github.com/asteca/asteca/issues/177); [#178](https://github.com/asteca/asteca/issues/178)).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.61]](https://github.com/asteca/asteca/releases/tag/v0.1.61) - 2015-03-04

### Changed

* Added [_"lazy parallelization"_](https://github.com/asteca/asteca/commit/b536c84c2ad085bbe8ff10a0b6535618ae1ba09a) ability. Now the user can run as many instances of the code as needed simply by creating extra `asteca_xx.py` and `input_xx` folders where `xx` are integers of the form: 01, 02,..., 99.
* [Reposition](https://github.com/asteca/asteca/commit/e7dec4b75a62ff397ee62cb322345f6b17b74ff6) several text boxes in output images, newer versions of `matplotlib` moved them from the previous position.
* Fix [bad import](https://github.com/asteca/asteca/commit/9bed2166e9cc36faa7077c79c436c50e40801820) of `rpy2` package, positioned incorrectly in two functions.
* Fix `DeprecationWarning` showing when `exp_function` was used ([#169](https://github.com/asteca/asteca/issues/169)).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.5]](https://github.com/asteca/asteca/releases/tag/v0.1.5) - 2015-03-03

### Changed

(**Warning**: this release breaks compatibility with the previous version of the `params_input.dat` file)

* Improved radius assignment algorithm ([#146](https://github.com/asteca/asteca/issues/146)).
* Detect cropped cluster region and use correct area when generating field regions ([#139](https://github.com/asteca/asteca/issues/139), [#157](https://github.com/asteca/asteca/issues/157)).
* Fixed bug that crashed the code when top tiers synthetic clusters with no stars were plotted ([#147](https://github.com/asteca/asteca/issues/147)). Added minimum total mass of 10Mo.
* Fixed bug where KDE p-values for field vs field comparision were artificially increased by comparing a field region with itself ([#138](https://github.com/asteca/asteca/issues/138)).
* Obtain KDE p-value even if one field region is defined ([#114](https://github.com/asteca/asteca/issues/114)).
* Fixed small bug that prevented integrated magnitude curves from being plotted ([#145](https://github.com/asteca/asteca/issues/145)).
* Fixed several smaller bugs and issues ([#110](https://github.com/asteca/asteca/issues/110), [#150](https://github.com/asteca/asteca/issues/150), [#140](https://github.com/asteca/asteca/issues/140), [#142](https://github.com/asteca/asteca/issues/142), [#141](https://github.com/asteca/asteca/issues/141), [#149](https://github.com/asteca/asteca/issues/149), [#95](https://github.com/asteca/asteca/issues/95), [#148](https://github.com/asteca/asteca/issues/148), [#136](https://github.com/asteca/asteca/issues/136), [#163](https://github.com/asteca/asteca/issues/163), [#143](https://github.com/asteca/asteca/issues/143)).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.4]](https://github.com/asteca/asteca/releases/tag/v0.1.4) - 2014-12-18

### Changed

* Improved plotting of crowded fields ([#62](https://github.com/asteca/asteca/issues/62)).
* Function to generate image is now more stable ([#112](https://github.com/asteca/asteca/issues/112)). Re-arranged plots in output image.
* Add _Top tiers_ models output ([#130](https://github.com/asteca/asteca/issues/130)).
* Fixed small bug in KDE p-values function ([#134](https://github.com/asteca/asteca/issues/134)).
* Minor re-arrangement with semi-input data.

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.3]](https://github.com/asteca/asteca/releases/tag/v0.1.3) - 2014-12-10

### Changed

* Accept arrays of non-equispaced parameter values instead of only equispaced ranges ([#121](https://github.com/asteca/asteca/issues/121)).
* Added support for lognormal [Chabrier (2001)](http://adsabs.harvard.edu/abs/2001ApJ...554.1274C) IMF.
* More precise encoding/decoding in genetic algorithm.
* Functions separated into sections ([#125](https://github.com/asteca/asteca/issues/125)).
* Input parameters set as global variables ([#132](https://github.com/asteca/asteca/issues/132)).

### Caveats

* Same as version [0.1.2](https://github.com/asteca/asteca/releases/tag/v0.1.2).

## [[v0.1.2]](https://github.com/asteca/asteca/releases/tag/v0.1.2) - 2014-12-01

### Changed

* Likelihood method now supports [Dolphin (2002)](http://adsabs.harvard.edu/abs/2002MNRAS.332...91D) _Poisson likelihood ratio_ function.
* Closed [#120](https://github.com/asteca/asteca/issues/120), [#101](https://github.com/asteca/asteca/issues/101), [#129](https://github.com/asteca/asteca/issues/129), [#124](https://github.com/asteca/asteca/issues/124), [#102](https://github.com/asteca/asteca/issues/102).
* Minor [position fix](https://github.com/asteca/asteca/commit/00538bda879009bae0a4e7565b124c8939c75d0f) for synthetic cluster text box in output plot.
* Brute force algorithm now returns [correct errors](https://github.com/asteca/asteca/commit/afe30cbdff561a90986a638c55a4b7247fd0bc53).
* Some fixes for when unique values in the input parameter ranges are used [[1]](https://github.com/asteca/asteca/commit/7cc383d799f2af5c1f1f8a6dcfc80e639461f02d) [[2]](https://github.com/asteca/asteca/commit/c6505025d4c3b6147a2913fad648dc18c125376b).
* Replaced deprecated [compiler package](https://github.com/asteca/asteca/commit/f9e8c5edba5f5ca8cc33ec1afb4d137f7167e8df) used to flatten list.

### Caveats

 * Still not sure why _tolstoy_ likelihood is biased towards high masses :confused:

## [[v0.1.1]](https://github.com/asteca/asteca/releases/tag/v0.1.1) - 2014-11-07

_More stable release._

### Changed

* Closed [#113](https://github.com/asteca/asteca/issues/113), [#116](https://github.com/asteca/asteca/issues/116).
* Minor [change](https://github.com/asteca/asteca/commit/3cffb4faa0c1dc6956aae2217c73afb4f392e53d) to error function.
* Closed _Known issues_ from previous version.

### Caveats

 * Same as previous version.

## [[v0.1.0]](https://github.com/asteca/asteca/releases/tag/v0.1.0) - 2014-10-08

_First <s>semi-stable</s> buggy release_

### Changed

* Closed [#72](https://github.com/asteca/asteca/issues/72), [#99](https://github.com/asteca/asteca/issues/99), [#37](https://github.com/asteca/asteca/issues/37).
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

* Closed [#85](https://github.com/asteca/asteca/issues/85), [#70](https://github.com/asteca/asteca/issues/70), [#43](https://github.com/asteca/asteca/issues/43), [#86](https://github.com/asteca/asteca/issues/86).
- Metallicity and age now take steps in the GA.
- Add [checker](https://github.com/Gabriel-p/asteca/blob/master/functions/checker.py) function to make sure certain parameters are set correctly before running.
* Number of points in `get_radius` increased 20% --> 25% of [the RDP](https://github.com/Gabriel-p/asteca/commit/a2e9b8f16111d5adafe66fed1eb64ed8bc03997b).

## [[v3.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v3.0.0-beta) - 2014-09-16

### Changed

* Closed: [#89](https://github.com/asteca/asteca/issues/89), [#77](https://github.com/asteca/asteca/issues/77), [#80](https://github.com/asteca/asteca/issues/80).
* The `params_input.dat` and `semi_input.dat` files are now located at the top level next to `asteca.py`.
* Cluster's photometric files are not longer required to be stored inside a sub-folder to be picked-up by the code.

## [[v2.0.1-beta]](https://github.com/asteca/asteca/releases/tag/v2.0.1-beta) - 2014-09-15

### Changed

* Correct version number.

## [[v2.0.0-beta]](https://github.com/asteca/asteca/releases/tag/v2.0.0-beta) - 2014-09-11

### Changed

* Closed issues: [#15](https://github.com/asteca/asteca/issues/15), [#73](https://github.com/asteca/asteca/issues/73), [#53](https://github.com/asteca/asteca/issues/53), [#24](https://github.com/asteca/asteca/issues/24),  [#75](https://github.com/asteca/asteca/issues/75), [#79](https://github.com/asteca/asteca/issues/79), [#81](https://github.com/asteca/asteca/issues/81), [#59](https://github.com/asteca/asteca/issues/59), [#83](https://github.com/asteca/asteca/issues/83), [#78](https://github.com/asteca/asteca/issues/78), [#69](https://github.com/asteca/asteca/issues/69), [#74](https://github.com/asteca/asteca/issues/74).
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

Version used (with some small changes) in the [original article](http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html).
