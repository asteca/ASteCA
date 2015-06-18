[![AA](https://goo.gl/YT7as5)][20]
[![License](http://goo.gl/38AkxE)][21]
[![Stories in Ready](https://goo.gl/Jg3mg0)][22]
[![Stories in Progress](https://goo.gl/XeAMah)][22]
________________________________________________________________________________


# ASteCA [Automated Stellar Cluster Analysis]

The [ASteCA][1] package is designed to fully automatize the usual tests
applied on [star clusters][2] in order to determine their characteristics:
center, radius, stars membership probabilities and associated
intrinsic/extrinsic parameters: metallicity, age, reddening, distance, total
mass, binarity fraction, etc.</p>

Read the code's [documentation][3] for details on how to use **ASteCA**.
Notice that the docs are still under development and many sections will be
incomplete or just empty, until I have the time to finish them.

The accompanying article describing the code in detail can be accessed
[via A&A][4], and referenced using the following BibTeX entry:

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

The code is still *under heavy development*. The latest release can be
accessed [here][5].

See the [CHANGELOG.md][6] file for a list of releases, the changes made in
each and possible compatibility breaks with the previous version.

## To do

* List of open [bugs][7].
* List of planned [enhancements][8].
* High priority [issues][9].
* All open [issues][10].

________________________________________________________________________________
[1]: http://asteca.github.io
[2]: https://en.wikipedia.org/wiki/Star_cluster
[3]: http://asteca.rtfd.org/
[4]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[5]: https://github.com/asteca/asteca/releases/latest
[6]: https://github.com/asteca/asteca/blob/active/CHANGELOG.md
[7]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[8]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[9]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Ap%3Ahigh
[10]: https://github.com/asteca/asteca/issues


#### License

If you distribute a copy or make a fork of the project, you have to credit
this project as source.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

Copyright (c) 2014 Gabriel Perren - Released under the GPL v3 license.

#### Versioning

Version numbering follows the [Semantic Versioning][19] guidelines. Releases
will be numbered
with the following format:

`<major>.<minor>.<patch>-<build>`

Constructed with the following guidelines:

* A new *major* release indicates a large change where backwards compatibility
  is broken.
* A new *minor* release indicates a normal change that maintains backwards
  compatibility.
* A new *patch* release indicates a bugfix or small change which does not
  affect compatibility.
* A new *build* release indicates this is a pre-release of the version.

[19]: http://semver.org/
[20]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[21]: http://www.gnu.org/licenses/gpl-3.0.en.html
[22]: http://waffle.io/asteca/asteca
