.. ASteCA documentation master file, created by
   sphinx-quickstart on Sat Feb 17 12:19:51 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/asteca_icon.png
   :scale: 50%
   :align: center

ASteCA
======

**ASteCA** is pure Python open source tool designed to perform the analysis applied to stellar clusters in order to determine their intrinsic/extrinsic fundamental parameters: extinction, distance, metallicity, age, binarity, mass, etc..


.. warning::
   Version 0.5.0 is a complete re-write of **ASteCA**. In previous versions it was
   a large Python script and from this version onward it is a proper Python package.
   This version also removed many of the features that were previously
   available, to concentrate on the fundamental parameters estimation.
   Feel free to `contact me`_ if you have questions about using this code in your
   research, or `open a new issue`_ in the code's repository.

.. _contact me: mailto:gabrielperren@gmail.com
.. _open a new issue: https://github.com/Gabriel-p/asteca/issues/new


.. toctree::
   :maxdepth: 2
   :caption: User guide

   contents/installation
   contents/faq

.. toctree::
   :maxdepth: 3
   :caption: Basic configuration

   basic/cluster_load
   basic/isochrones_load
   basic/synthetic_clusters

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/quickstart

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api/cluster
   api/isochrones
   api/synthcluster
   api/likelihood



License & Attribution
---------------------

Copyright 2015-2024 Gabriel I Perren.

**ASteCA** is free software made available under the MIT License. For details
see the `LICENSE`_.

The accompanying article describing the code in detail can be accessed
via `A&A <http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html>`_,
and referenced using the following BibTeX entry:

.. code-block:: Bibtex

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


.. _LICENSE: https://github.com/asteca/ASteCA/blob/master/LICENSE.txt


Changelog
---------

.. include:: CHANGELOG.rst

