# OCAAT

This code is still under heavy development.

See the code's [site][1] for more information.

## Software

The code has been tested with the latest 2.7.x `python` version:

* [Python - 2.7.8][2]

It should also work on older versions. Iif you encounter any problems with older versions, please contact me or open a new issue.

The best and easiest way to install several versions of `python` and its packages without affecting your system is [pyenv](https://github.com/yyuu/pyenv).

### Python dependencies

The packages listed below are required to run the code.

* [Numpy - 1.8.2][4] -- `pip install numpy`
* [Matplotlib - 1.3.1][6] -- `pip install matplotlib`
* [SciPy  - 0.14.0][5] -- `pip install scipy`

The versions are the ones I used, could work with older versions but I can't guarantee it.

#### Extra packages
If you want to use the function to obtain the cluster probability of being a true cluster, the following packages are needed:

* [R - 3.0.1][3] -- `sudo apt-get install r-base`
* [rpy2  -2.4.3](http://rpy.sourceforge.net/) -- `pip install rpy2`

These packages are not mandatory and OCAAT will still run without them, just not that particular function.

### Theoretical isochrones

OCAAT needs at least one set of theoretical isochrones stored in a `/isochrones/`
folder to be able to apply the function that estimates the clusters' parameters. [Girardi isochrones][7] are currently the default but any set can in practice be 
used (with minor changes made to the code)

The isochrones can be downloaded manually or the package [ezPadova-2][8] can be
used to automatically fetch them from the site.

## To do

* List of open [bugs][9].
* List of planned [enhancements][10].
* High priority [issues][11].
* All open [issues][12].


[1]: http://gabriel-p.github.io/ocaat/
[2]: www.python.org
[3]: http://www.r-project.org/
[4]: http://www.numpy.org/
[5]: http://www.scipy.org/
[6]: http://matplotlib.org/
[7]: http://stev.oapd.inaf.it/cgi-bin/cmd
[8]: https://github.com/Gabriel-p/ezpadova
[9]: https://github.com/Gabriel-p/ocaat/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[10]: https://github.com/Gabriel-p/ocaat/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[11]: https://github.com/Gabriel-p/ocaat/issues?q=is%3Aopen+is%3Aissue+label%3Aprior%3Ahigh
[12]: https://github.com/Gabriel-p/ocaat/issues
