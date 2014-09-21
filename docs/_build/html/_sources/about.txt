About
=====

ASteCA (Automated Stellar Cluster Analysis) is an open source code
developed entirely in python designed to fully automatize the usual tests
applied to stellar clusters (SC) in order to determine their characteristics:
center, radius, stars' membership probabilities and associated
intrinsic/extrinsic parameters: metallicity, age, extinction, distance
and eventually total mass.

The code is designed modularly thus allowing the user to select which
tests to run and which to skip, all managed through a simple input data
file. Thanks to this modularity functions can even be swapped, added or
removed so, for example, more than one decontamination algorithm can be
made available for the user to choose.

In its current version (v1.0.0) ASteCA can be applied to a number of color
magnitude diagrams (CMD) in the systems of Johnson, Washington and 2MASS
but it can be easily extended to work on any number of CMDs. Eventually
the code can be generalized to work with an arbitrary nubmer of colors
instead of being limited to a single color, like it is when using a CMD.

ASteCA is intended to run without the need for user intervention but both
a semi-automatic and a manual mode are made available in case user input
is necessary or required.

The sections below present a brief description of the processes applied by
ASteCA to an input photometric data file (ie: an observed cluster) in the
order in which they run. More information about each of these processes
can be found in the code's user manual.
