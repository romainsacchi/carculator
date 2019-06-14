
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

    pip install git+https://github.com/romainsacchi/coarse.git

will install the package and the required dependencies.


How to use it?
--------------

From a Python environment::

   from coarse import *
   cip = CarInputParameters()
   cip.stochastic(10)

This would effectively load all vehicle profiles as well as generate 10 random values for each parameters, that can be further used for error propagation analyses.

Alternatively, if one wishes to work with static parameters values instead::

    cip.static()

Then::

   dcts, array = fill_xarray_from_input_parameters(cip)

will generate a 3-dimensional array `array` to store the generated parameters values, with the following dimensions:

0. Vehicle size, e.g. "small", "medium". str.
1. Powertrain, e.g. "ICE-p", "BEV". str.
2. Year. int.
3. Samples.

Finally::

   cm = CarModel(array)
   cm.set_all()

generate a CarModel object and calculate the energy consumption and components mass for all vehicle profiles.


.. automodule:: coarse.car_input_parameters
    :members:



