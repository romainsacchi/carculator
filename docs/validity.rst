Validity tests
==============

Driving cycle, velocity and acceleration
----------------------------------------

There are eleven possible driving cycles to select from:

    * WLTC
    * WLTC 3.1
    * WLTC 3.2
    * WLTC 3.3
    * WLTC 3.4
    * CADC Urban
    * CADC Road
    * CADC Motorway
    * CADC Motorway 130
    * CADC
    * NEDC

Velocity, driving distance, driving time, and acceleration are calculated.

Manually, such parameters can be obtained the following way::
    
    import pandas as pd
    import numpy as np
    # Retrieve the driving cycle WLTC 3 from the UNECE
    driving_cycle = pd.read_excel('http://www.unece.org/fileadmin/DAM/trans/doc/2012/wp29grpe/WLTP-DHC-12-07e.xls',
                sheet_name='WLTC_class_3', skiprows=6, usecols=[2,4,5])

    # Calculate velocity (km/h -> m/s)
    velocity = driving_cycle['km/h'].values * 1000 / 3600

    # Retrieve driving distance (-> km)
    driving_distance = velocity.sum() * 1000

    # Retrieve driving time (-> s)
    driving_time = len(driving_cycle.values)

    # Retrieve acceleration by calculating the delta of velocity per time interval of 2 seconds
    acceleration = np.zeros_like(velocity)
    acceleration[1:-1] = (velocity[2:] - velocity[:-2])/2

Using `carculator`, these parameters can be obtained the following way::

    from coarse.energy_consumption import EnergyConsumptionModel
    ecm = EnergyConsumptionModel('WLTC')

    # Access the driving distance
    ecm.velocity.sum() * 1000

    # Access the driving time
    len(ecm.velocity)

    # Access the acceleration
    ecm.acceleration
    
Both approaches should return identical results::

    print(np.array_equal(velocity, ecm.velocity))
    print(driving_distance == ecm.velocity.sum()*1000)
    print(driving_time == len(ecm.velocity))
    print(np.array_equal(acceleration, ecm.acceleration))

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/coarse.png
    :width: 900
    :alt: Alternative text

Modules
-------

Composed of four modules:

    * Driving cycle module
    * Mass module
    * Auxiliary energy module
    * Motive energy module
    
Driving cycle module
--------------------

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/driving_cycle.png
    :width: 400
    :alt: Alternative text
    
Mass module
-----------

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/mass_module.png
    :width: 900
    :alt: Alternative text
    
Auxiliary energy module
-----------------------

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/aux_energy.png
    :width: 900
    :alt: Alternative text
    
Motive energy module
--------------------

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/motive_energy.png
    :width: 900
    :alt: Alternative text
