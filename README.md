[![AA](https://img.shields.io/badge/A%26A-576--A6%2C%202015-yellowgreen.svg)](http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html) [![License](http://img.shields.io/badge/license-GPLv3-red.svg)](http://www.gnu.org/licenses/gpl-3.0.en.html) [![Stories in Ready](https://badge.waffle.io/asteca/asteca.svg?label=ready&title=Ready)](http://waffle.io/asteca/asteca) [![Stories in Progress](https://badge.waffle.io/asteca/asteca.svg?label=in_prog&title=In%20Progress)](http://waffle.io/asteca/asteca)

# ASteCA [Automated Stellar Cluster Analysis]

The [ASteCA][3] package is designed to fully automatize the usual tests applied on [star clusters][12] in order to determine their characteristics: center, radius, stars membership probabilities and associated intrinsic/extrinsic parameters: metallicity, age, reddening, distance, total mass, binarity fraction, etc.</p>

Read the code's [documentation][9]. <-- *Still under development*.

The accompanying article describing the code in detail can be accessed [via A&A][10],
and referenced using the following BibTeX entry:

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


## Releases

The code is still under (heavy) development. The latest release can be accessed [here][1].

See the [CHANGELOG.md][11] file for a list of releases, the changes made in each and possible compatibility breaks with the previous version.

## To do

* List of open [bugs][4].
* List of planned [enhancements][5].
* High priority [issues][6].
* All open [issues][7].

***

## Versioning

Version numbering follows the [Semantic Versioning][8] guidelines. Releases will be numbered
with the following format:

`<major>.<minor>.<patch>-<build>`

Constructed with the following guidelines:

* A new *major* release indicates a large change where backwards compatibility is broken.
* A new *minor* release indicates a normal change that maintains backwards compatibility.
* A new *patch* release indicates a bugfix or small change which does not affect compatibility.
* A new *build* release indicates this is a pre-release of the version.

***

If you distribute a copy or make a fork of the project, you have to credit this project as source.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

Copyright (c) 2014 Gabriel Perren - Released under the GPL v3 license.

[1]: https://github.com/asteca/asteca/releases/latest
[3]: http://asteca.github.io
[4]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[5]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[6]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Ap%3Ahigh
[7]: https://github.com/asteca/asteca/issues
[8]: http://semver.org/
[9]: http://asteca.rtfd.org/
[10]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[11]: https://github.com/asteca/asteca/blob/master/CHANGELOG.md
[12]: https://en.wikipedia.org/wiki/Star_cluster
