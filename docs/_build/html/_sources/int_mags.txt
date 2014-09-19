Integrated magnitudes
=====================

To obtain the integrated magnitude for the cluster region *cleaned* from
field stars contamination, the following equations are used.

First, assume:

.. math::
   I_{cl+fl} = I_{cl} + I_{fl}
   :label: integ-cl-fl

then:

.. math:: V^*_{cl} - V^*_{cl+fl} = -2.5 \log(I_{cl}/I_{cl+fl})

using Eq. :eq:`integ-cl-fl` we have:

.. math::
   V^*_{cl} - V^*_{cl+fl} = -2.5 \log(1 - I_{fl}/I_{cl+fl})
   :label: integ-cl-fl2

and given:

.. math::
   V^*_{fl} - V^*_{cl+fl} = -2.5 \log(I_{fl}/I_{cl+fl}) \Rightarrow \frac{I_{fl}}{I_{cl+fl}} = 10^{(V^*_{fl} - V^*_{cl+fl})/-2.5}
   :label: integ-cl-fl3

we can combine now Eqs. :eq:`integ-cl-fl2` and :eq:`integ-cl-fl3` to obtain:

.. math:: V^*_{cl} = -2.5 \log(1 - 10^{(V^*_{fl} - V^*_{cl+fl})/-2.5}) + V^*_{cl+fl}