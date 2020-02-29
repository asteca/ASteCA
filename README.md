[![AA](https://img.shields.io/badge/A%26A-576--A6,%202015-yellowgreen.svg)][12]
[![License](https://img.shields.io/badge/license-GPLv3-red.svg)][13]
________________________________________________________________________________


# ASteCA [Automated Stellar Cluster Analysis]

The [ASteCA][1] package is designed to fully automatize the usual tests
applied on [star clusters][2] in order to determine their characteristics:
center, radius, stars membership probabilities and associated
intrinsic/extrinsic parameters: metallicity, age, reddening, distance, total
mass, binarity fraction, etc.

**IMPORTANT**: until the release of v1.0.0 the package will be *under heavy
development*. Keep this in mind if you want to use it in your research.


## Releases

The latest release can always be accessed [here][5].

See the [CHANGELOG][6] file for a list of previous releases, the changes
made in each one, and possible compatibility breaks with older versions.


## Requirements

ASteCA requires the following packages:

* [`astropy`][14], [`numpy`][15], [`scipy`][16], [`matplotlib`][17], [`ptemcee`][18], [`emcee`][19], [`corner.py`][20]

The `emcee` package is not essential but some functions will not work without it. The packages `ptemcee` and `corner.py` are included within the code for simplicity, and don't need to be installed.



## Referencing

The accompanying article describing the code in detail can be accessed
[via A&A][7], and referenced using the following BibTeX entry:

````
@article{Perren_2015,
    author = {{Perren, G. I.} and {V\'azquez, R. A.} and {Piatti, A. E.}},
    title = {ASteCA: Automated Stellar Cluster Analysis},
    DOI= "10.1051/0004-6361/201424946",
    url= "http://dx.doi.org/10.1051/0004-6361/201424946",
    journal = {A\&A},
    year = 2015,
    volume = 576,
    pages = "A6",
    month = "04",
}
````

This article will rapidly become outdated as new versions are released. For an
up to date description of **ASteCA** please refer to the online
[documentation][4].


## To do

* List of open [bugs][8].
* List of planned [enhancements][9].
* High priority [issues][10].
* All open [issues][11].

________________________________________________________________________________
[1]: http://asteca.github.io
[2]: https://en.wikipedia.org/wiki/Star_cluster
[4]: http://asteca.rtfd.org/
[5]: https://github.com/asteca/asteca/releases/latest
[6]: https://github.com/asteca/ASteCA/blob/master/CHANGELOG.md
[7]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[8]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[9]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[10]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Ap%3Ahigh
[11]: https://github.com/asteca/asteca/issues
[12]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[13]: http://www.gnu.org/licenses/gpl-3.0.en.html
[14]: https://www.astropy.org/
[15]: https://numpy.org/
[16]: https://www.scipy.org/
[17]: https://matplotlib.org/
[18]: https://github.com/willvousden/ptemcee
[19]: https://emcee.readthedocs.io/
[20]: https://corner.readthedocs.io/
