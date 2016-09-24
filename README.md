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
mass, binarity fraction, etc.

**IMPORTANT**: until the release of v1.0.0 the package will be *under heavy
development*. Keep this in mind if you want to use it on your research.


## Running

After downloading the package, and assuming all necessary dependencies are
installed (see the code's [requirements][3]), simply open a terminal, move
into the unpacked **ASteCA** folder and run:

````
$ python asteca.py
````

Read the code's [documentation][4] for more details on how to use **ASteCA**.

Notice that, as the code itself, the docs are *still under development* and
many sections will be incomplete or just empty. The docs will be improved as
new versions are released.


## Releases

The latest release can always be accessed [here][5].

See the [CHANGELOG][6] file for a list of previous releases, the changes
made in each one, and possible compatibility breaks with older versions.


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
up to date description of **ASteCA** please refer to the online documentation.


## To do

* List of open [bugs][8].
* List of planned [enhancements][9].
* High priority [issues][10].
* All open [issues][11].

________________________________________________________________________________
[1]: http://asteca.github.io
[2]: https://en.wikipedia.org/wiki/Star_cluster
[3]: http://asteca.readthedocs.org/en/latest/requirements.html
[4]: http://asteca.rtfd.org/
[5]: https://github.com/asteca/asteca/releases/latest
[6]: https://github.com/asteca/ASteCA/blob/master/CHANGELOG.md
[7]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[8]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Abug
[9]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement
[10]: https://github.com/asteca/asteca/issues?q=is%3Aopen+is%3Aissue+label%3Ap%3Ahigh
[11]: https://github.com/asteca/asteca/issues


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

________________________________________________________________________________
[19]: http://semver.org/
[20]: http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html
[21]: http://www.gnu.org/licenses/gpl-3.0.en.html
[22]: http://waffle.io/asteca/asteca
