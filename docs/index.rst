
Coarse Documentation
====================

**coarse** allows to calculate life cycle inventories for 112 different car profiles, according to selected:

* powertrain technologies (8): petrol engine, diesel engine, electric motor, hybrid, etc.,
* year of operation (2): 2017, 2040
* and sizes (7): Mini, Large, etc.

It is based on the model developed in `Uncertain environmental footprint of current and future battery electric
vehicles by Cox, et al (2018) <https://pubs.acs.org/doi/10.1021/acs.est.8b00261>`_.

Objective
---------

The objective is to produce life cycle inventories in a transparent and comprehensive way, to be further used in prospective
life cycle assessment of transportation technologies.

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

Hence, the tank-to-wheel energy requirement per km driven per powertrain technology for a SUV in 2017 can be obtained::

    TtW_energy = cm.array.sel(size='SUV', year=2017, parameter='TtW energy').values * 1/3600 * 100
    powertrain = cm.array.powertrain

    import numpy as np
    plt.bar(powertrain,[a[0] for a in TtW_energy])
    plt.ylabel('kWh/100 km')
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/blob/master/docs/fig_kwh_100km.png
  :width: 400
  :alt: Alternative text

.. automodule:: coarse.car_input_parameters
    :members:



