# OCAAT

This code is still under heavy development.

See the code's [site][1] for more information.

## Software

Both `python` and `R` are needed to run OCAAT. I've tested the code
using the versions shown below.

* [Python v2.7.3][2]
* [R v3.0.2][3]

### Python dependencies

The packages listed are required to run the code. The versions are the ones I used,
could work with older versions but I can't guarantee it.

* [Numpy 1.8.0][4] -- `sudo pip install numpy`
* [SciPy 0.12.0][5] -- `sudo pip install scipy`
* [Matplotlib 1.2.1][6] -- `sudo pip install matplotlib`

### Theoretical isochrones

OCAAT needs at least one set of theoretical isochrones stored in a `/isochrones/`
folder to be able to apply the function that estimates the clusters' parameters. [Girardi isochrones][8] are currently the default but any set can in practice be 
used (with minor changes made to the code)

The isochrones can be downloaded manually or the package [ezPadova-2][9] can be
used to automatically fetch them from the site.

## To do

List of planned and still open [enhancements][7].


[1]: http://gabriel-p.github.io/ocaat/
[2]: www.python.org
[3]: http://www.r-project.org/
[4]: http://www.numpy.org/
[5]: http://www.scipy.org/
[6]: http://matplotlib.org/
[7]: https://github.com/Gabriel-p/ocaat/search?o=asc&p=1&q=label%3Aenhancement+state%3Aopen&s=created&type=Issues
[8]: http://stev.oapd.inaf.it/cgi-bin/cmd
[9]: https://github.com/Gabriel-p/ezpadova