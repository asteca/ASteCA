Adding CMD support
==================

In order to add support for a new CMD the following functions should be
modified:

#. ``get-in-params``, add names of filters and positions needed to form
   the color.

#. ``get-isoch-params``, add the columns in the theoretical isochrone
   where the code should look for the filters used.

#. ``move-isochrone``, add the extinction equations that define how this
   effect along with the distance modulus affect the magnitude and color
   in the CMD.
