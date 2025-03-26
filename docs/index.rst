.. ASteCA documentation master file, created by
   sphinx-quickstart on Sat Feb 17 12:19:51 2024.

.. image:: _static/asteca_icon.webp
   :scale: 10%
   :align: center
   :class: no-scaled-link

ASteCA
======

**ASteCA** is an open-source tool developed in Python for the analysis of stellar clusters. It is designed to determine their structural parameters, membership probabilities, and both intrinsic and extrinsic fundamental properties, including extinction, distance, metallicity, age, binarity, and mass.


.. important::
   Version |ProjectVersion| released on the 24th of January, 2025. See :ref:`changelog`
   for a detailed list of the changes implemented.

Feel free to `contact me`_ if you have questions about using this code in your
research, and please `open a new issue`_ in the code's repository if you find
something either wrong, confusing or missing.


.. _contact me: mailto:gabrielperren@gmail.com
.. _open a new issue: https://github.com/Gabriel-p/asteca/issues/new


License & Attribution
=====================

Copyright 2015-2024 Gabriel I Perren.

**ASteCA** is free software made available under the MIT License. For details
see the `LICENSE`_.

If you make use of **ASteCA** in your research, please cite its `original
article <http://www.aanda.org/articles/aa/abs/2015/04/aa24946-14/aa24946-14.html>`_
using the following BibTeX entry:

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


TOC
===

.. toctree::
   :maxdepth: 1
   :caption: User guide

   contents/installation
   contents/changelog
   contents/api_index


.. toctree::
   :maxdepth: 1
   :caption: Basic configuration

   basic/cluster_load
   basic/membership
   basic/isochrones_load
   basic/synthetic_clusters

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/membership
   tutorials/paramfit



