
Coarse Documentation
====================

**coarse** allows to calculate life cycle inventories for vehicles of different:

* powertrain technologies,
* year of operation
* and sizes.

It is based on the model developed in `Uncertain environmental footprint of current and future battery electric
vehicles by Cox, et al (2018) <https://pubs.acs.org/doi/10.1021/acs.est.8b00261>`_.

Objective
---------
The objective is to produce life cycle inventories in a transparent way, to be further used in prospective
life cycle assessment of mobility technologies.

How to install this package?
----------------------------

In a terminal::

    pip install git+(Git repository address)

How to use it?
--------------

From a Python environment::

   from coarse import *
   cip = CarInputParameters()
   cip.stochastic(10)
   dcts, array = fill_xarray_from_input_parameters(cip)
   cm = CarModel(array)
   cm.set_all()




.. toctree::
   :maxdepth: 4
   :caption: Table of contents
