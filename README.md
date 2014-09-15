# ASteCA

The code is still under development. Version 1.0.0 will be released as soon as the following list of issues is addressed:

* [Issues left until version 1.0.0](https://github.com/Gabriel-p/asteca/milestones/v1.0.0)

A list of **Beta** releases can be accessed [here](https://github.com/Gabriel-p/asteca/releases).

See the code's [site][1] for more information. <-- *Under development*.

## Contents

- [Installation](#installation)
 - [Python](#python)
 - [Dependencies](#dependencies)
 - [Extra packages](#extra-packages)
 - [Python distributions](#python-distributions)
 - [Theoretical isochrones](#theoretical-isochrones)
- [To do](#to-do)
- [Versioning](#versioning)

<!-- end toc -->

## Installation

The packaged zip or tarball for the latest release can be downloaded from:

* [.zip](https://github.com/Gabriel-p/asteca/releases)
* [.tar.gz](https://github.com/Gabriel-p/asteca/releases)

You can alternatively checkout the entire project via `git`. Simply locate
yourself in the folder you want the code to be downloaded and run:

    git clone https://github.com/Gabriel-p/asteca.git

This will create a new sub-folder named _/asteca_ with all the code
stored inside.

### Python

The code has been tested with the July 2014 release of [python](www.python.org):

* [Python - 2.7.8](https://www.python.org/download/releases/2.7.8/)

It should also work on older versions, if you encounter any problems please contact me or [open a new issue](https://github.com/Gabriel-p/asteca/issues/new).

### Dependencies

The packages listed below are required to run ASteCA.

* [Numpy - 1.8.2][4] -- `pip install numpy`
* [Matplotlib - 1.3.1][6] -- `pip install matplotlib`
* [SciPy  - 0.14.0][5] -- `pip install scipy`

These versions are the ones I used, the code could work with older versions of the
packages but I can't guarantee it.

### Python distributions

An alternative to installing packages separately is to download a Python distribution which comes with many packages already installed:

* [Anaconda](https://store.continuum.io/cshop/anaconda/)
* [Canopy](https://www.enthought.com/products/canopy/)

The best and easiest way to install and manage several versions of python and its packages without affecting your system is [pyenv](https://github.com/yyuu/pyenv).

### Extra packages
If you want to use the [function](https://github.com/Gabriel-p/asteca/blob/master/functions/get_p_value.py) that obtains the cluster probability of being a true cluster, the following software is needed:

* [R - 3.0.1][3] -- `sudo apt-get install r-base`

After it is installed, open with sudo privileges (`sudo R`) and install, _in order_, the packages:

* rgl -- `install.packages("rgl")`
* mvtnorm -- `install.packages("mvtnorm")`
* misc3d -- `install.packages("misc3d")`
* ks -- `install.packages("ks")`

The package that allows `python` and `R` to coomunicate is also needed:

* [rpy2  -2.4.3](http://rpy.sourceforge.net/) -- `pip install rpy2`

These extra packages are not mandatory and ASteCA will still run without them, just not that particular function.

### Theoretical isochrones

ASteCA needs at least one set of theoretical isochrones stored in a `/asteca/isochrones` folder to be able to apply the function that estimates the clusters' parameters.
[Girardi isochrones][7] are currently the default but any set can in practice be used (with minor changes made to the code)

The isochrones can be downloaded manually or the package [ezPadova-2][8] can be used to automatically fetch them from the site.

## To do

* List of open [bugs][9].
* List of planned [enhancements][10].
* High priority [issues][11].
* All open [issues][12].

***

## Versioning

Version numbering follows the [Semantic Versioning](http://semver.org/) guidelines. Releases will be numbered with the following format:

`<major>.<minor>.<patch>-<build>`

Constructed with the following guidelines:

* A new *major* release indicates a large change where backwards compatibility is broken.
* A new *minor* release indicates a normal change that maintains backwards compatibility.
* A new *patch* release indicates a bugfix or small change which does not affect compatibility.
* A new *build* release indicates this is a pre-release of the version.


[1]: http://gabriel-p.github.io/asteca/
[3]: http://www.r-project.org/
[4]: http://www.numpy.org/
[5]: http://www.scipy.org/
[6]: http://matplotlib.org/
[7]: http://stev.oapd.inaf.it/cgi-bin/cmd
[8]: https://github.com/Gabriel-p/ezpadova
[9]: https://github.com/Gabriel-p/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[10]: https://github.com/Gabriel-p/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[11]: https://github.com/Gabriel-p/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Ap%3Ahigh
[12]: https://github.com/Gabriel-p/asteca/issues
