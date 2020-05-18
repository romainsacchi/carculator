Modeling and assumptions
========================

The modeling of passenger vehicles in the past, present and future is complex and relies on many assumptions.
With **carculator**, we wish to be transparent about those: assumptions and modeling approaches should ideally be easily
critiqued and modified.

We there try to give a comprehensive list of assumptions and modeling choices on this page, and describe how, as a user, you
can change those.

Parameters' names are indicated ``verbatim`` and are to be used in **carculator**.

Vehicle sizing
**************
**carculator** models vehicles along four dimensions:

* their powertrain (e.g., gasoline-run internal combustion engine, battery electric vehicle, etc.),
* their size (e.g., mini, medium, large, etc.),
* their year of production (2000, 2010, 2017 and 2040)
* and a parameter dimension (i.e., input and calculated parameters).

When **carculator** sizes the vehicles for the different powertrains, sizes and years, it starts with the
input parameter's value for the `glider base mass`, which is essentially an initial guess for the mass of the vehicle's
glider without anything on it.

Then it adds the following components and their associated mass:

* ``fuel mass``: mass of the fuel in the fuel tank (only applicable to vehicles using liquid or gaseous fuels),
* ``fuel tank mass``: mass of the fuel tank (empty),
* ``charger mass``: mass of the onboard battery charger (for battery electric and plugin hybrid vehicles only),
* ``converter mass``: mass of the onboard electricity AC/DC converter (for battery electric and plugin hybrid vehicles only),
* ``inverter mass``: mass of the onboard electricity DC/AC converter (for battery electric and plugin hybrid vehicles only),
* ``power distribution unit mass``: mass of the onboard power distribution unit (for battery electric and plugin hybrid vehicles only),
* ``combustion engine mass``: mass of the internal combustion engine (if applicable),
* ``electric engine mass``: mass of the electric motor (if applicable),
* ``powertrain mass``: mass of the powertrain excluding the mass of the engine (e.g., transmission, drive shafts, differentials, etc.),
* ``fuel cell stack mass``: mass of the fuel cell stack (only for fuel cell electric vehicles),
* ``fuel cell ancillary BoP mass``: mass of the ancillary part of the Balance of Plant of the fuel cell stack (only for fuel cell electric vehicles),
* ``fuel cell essential BoP mass``: mass of the essential part of the Balance of Plant of the fuel cell stack (only for fuel cell electric vehicles),
* ``battery cell mass``: mass of the battery cells. Two types of batteries are distinguished: power and energy batteries,
* ``battery BoP mass``: mass of the Balance of Plant of the battery.


Adding the mass of the glider to the mass of these components constitutes a first attempt at guessing the ``curb mass`` of
the vehicle, that is its mass in working order, but without passengers and cargo.
The ``driving mass`` of the vehicle is then obtained by summing the ``curb mass`` to the mass of the passengers
(``average passengers`` x ``average passenger mass``) and cargo transported (``cargo mass``).

A second step consists into calculating the mass of the combustion and electric engine, based on the following relations:

    power demand (``power``) [kW] = ``power-to-mass ratio`` [kW/kg] x ``curb mass`` [kg]

    electrical power demand (``electric power``) [kW] = power demand (``power``) [kW] x (1 - ``combustion power share`` [%])

    ``electric engine mass`` [kW] = (``electric power`` [kW] x ``electric mass per power`` [kg/kW]) + ``electric fixed mass`` [kg]

    combustion power demand (``combustion power``) [kW] = ``power`` [kW] x ``combustion power share`` [%]

    ``combustion engine mass`` [kW] = (``combustion power`` [kW] x ``combustion mass per power`` [kg/kW]) + ``combustion fixed mass`` [kg]


As well as for the mass of the powertrain:

    ``powertrain mass`` [kg] = (``power`` [kW] x ``powertrain mass per power`` [kg/kW]) + ``powertrain fixed mass`` [kg]

With the mass of these new components recalculated (``electric engine mass``, ``combustion engine mass`` and ``powertrain mass``),
the curb mass of the vehicle is calculated once again. This iterative process stops when the curb mass of the vehicle
stabilizes (i.e., when recalculating the mass of the engine and powertrain does not lead to a change in the new curb
mass of more than one percent).

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/mass_module.png
    :width: 900
    :alt: Structure of mass module

Four initial input parameters are therefore of importance:
* `glider base mass`:the initial mass of the glider
* `power to mass ratio`: the power-to-mass ratio
* `combustion power share`: how much of the power is provided by an internal combustion engine
* `combustion mass per power`: the mass of the combustion engine per unit of power

For electric vehicles (i.e., BEV and FCEV), ``combustion power share`` = 0.
For internal combustion engine vehicles (i.e., ICEV-p, ICEV-d and ICEV-g),
``combustion power share`` = 1 in the early years (until 2020). However, starting 2020 on, this value drops progressively
to 0.85 by 2050, as we assumed a mild-hybridization of the powertrain to a level similar to that of non-plugin hybrids nowadays (i.e., HEV-p and HEV-d).
While it is uncertain whether ICEVs will exist in the future, it was assumed that a way for them to comply with future
emission standards was to be assisted by an electric engine. This mild-hybridization allows to reduce the size of the combustion engine and recover energy during braking.

For non-plugin hybrids, ``combustion power share`` is usually set at around 0.75.

For plugin hybrid vehicles, things are modeled differently: a purely electric vehicle is modeled, as well as a purely
combustion-based vehicle. Later on, when the range of the purely-electric vehicle is calculated, a ``electric utility ratio``
is obtained, which is used to fusion both vehicles. This ratio, which is dependent on the range, is usually between 0.6 and 0.7.
This means that plugin hybrid vehicles are made of between 60 and 70% of a purely electric vehicle and 30 to 40% of a purely combustion-based vehicle.

If I know already the ``curb mass`` of a vehicle, can I override its value?
---------------------------------------------------------------------------
Yes. After having created the CarModel() object and executed the :meth:`.set_all` method, you can override the
calculated ``curb mass`` value. Here is an example for a diesel car of medium size in 2020::

    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    cm.array.loc[dict(parameter="curb mass",
                  powertrain="ICEV-d",
                  year=2020,
                  size="Medium")] = 1600

How to prevent the mild-hybridization of ICEVs?
-----------------------------------------------

With **carculator online**:

In the Parameters section, search for `combustion power share` and add the parameter for the vehicles you wish to modify.

With **carculator**:

You can simply override the default value by "1" in ``array`` before passing it to CarModel()::

    dict_param = {('Powertrain',  ('ICEV-d', 'ICEV-p', 'ICEV-g'), 'all', 'combustion power share', 'none'): {
                                                                                        (2000, 'loc'): 1,
                                                                                        (2010, 'loc'): 1,
                                                                                        (2017, 'loc'): 1,
                                                                                        (2040, 'loc'): 1}
                                                                                        }
    modify_xarray_from_custom_parameters(dict_param, array)

You can also just override the default value of a specific powertrain of a specific size, for a specific year::

    dict_param = {('Powertrain',  'ICEV-d', 'Medium', 'combustion power share', 'none'): {
                                                                                        (2017, 'loc'): 1
                                                                                        }
    modify_xarray_from_custom_parameters(dict_param, array)

The ``curb mass`` values obtained for the vehicles in 2000, 2010 and 2017 are calibrated against a passenger cars database
`Car2DB <https://car2db.com/>`_. The calibration of the ``curb mass`` for vehicles for the year 2000 is done against vehicles in
the Car2DB database with a production year in the range of 1998-2002, against 2008-2012 and 2015-2018 for vehicles for the years
2010 and 2017, respectively.
The value of the input parameter ``glider base mass`` was adjusted to fit the distribution shown in the plots below.

Calibration of vehicles' curb mass for the year 2000

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/curb_mass_calibration_2000.png
    :width: 900
    :alt: Calibration for year 2000 vehicles

Calibration of vehicles' curb mass for the year 2010

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/curb_mass_calibration_2010.png
    :width: 900
    :alt: Calibration for year 2010 vehicles

Calibration of vehicles' curb mass for the year 2017

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/mass_comparison.png
    :width: 900
    :alt: Calibration for year 2017 vehicles

For the year 2040, the value for input parameters ``glider base mass``, ``combustion mass per power``, ``power to mass ratio`` are
adjusted according to the following studies:

* Hirschberg (Editor) S, Bauer C, Cox B, Heck T, Hofer J, Schenler W, et al. Opportunities and challenges for electric mobility: an interdisciplinary assessment of passenger vehicles Final report of the THELMA project in co-operation with the Swiss Competence Center for Energy Research "Efficient technologies and systems for mobil. 2016.
* Del Duce, Andrea; Gauch, Marcel; Althaus, Hans-Jörg: "Electric passenger car transport and passenger car life cycle inventories in ecoinvent version 3", International Journal of Life Cycle Assessment, Vol. 21, pp. 1314-1326, (2016)
* E. A. Grunditz and T. Thiringer, "Performance Analysis of Current BEVs Based on a Comprehensive Review of Specifications," in IEEE Transactions on Transportation Electrification, vol. 2, no. 3, pp. 270-289, Sept. 2016, doi: 10.1109/TTE.2016.2571783.

What happens what I inter-/extrapolate to other years?
------------------------------------------------------

If the default years of 2000, 2010, 2017 and 2040 are of no interest, it is possible to inter-/extrapolate the vehicle
models to any year between 2000 and 2050. When such inter-/extrapolation is done, all the *physical* input parameters' values
are inter-/extrapolated **linearly**.

With **carculator online**:
Simply drag the desired years from the left frame to the right frame.

With **carculator**:
After creating ``array``, which is a `DataArray` object from the library ``xarray``, it is possible to use the `.interp()`
method, like so::

     array = array.interp(year=np.arange(2015, 2051, 5),  kwargs={'fill_value': 'extrapolate'})

Here, the years under study are from 2015 to 2050 by step of 5 years.

This is slightly different for cost input parameters' values, which are usually following a decay-like cost curve, to account
for a learning rate.
Hence, parameters such as ``fuel tank cost per kg``, ``fuel cell cost per kW``, ``energy battery cost per kWh``, ``power battery cost per kW``,
or ``combustion powertrain cost per kW`` would be of shape: a*exp(b) + c. Coefficients *a*, *b* and *c* are defined to fit the literature and projections.

Projection of energy battery cost per kWh for BEV and FCEV.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/cost_energy_battery_projection.png
    :width: 900
    :alt: Projection of energy battery cost per kWh


Tank-to-wheel energy consumption
********************************
Once the vehicle and its powertrain has been sized, it is possible to calculate the motive energy required along
a specific drivign cycle to overcome the following forces:

* rolling resistance
* aerodynamic resistance
* air resistance
* road gradient resistance (if provided)

on top of the *kinetic energy* needed to move the vehicle.

To that amount of energy is subtracted the *energy recuperated* during braking, if the vehicle is equipped with
an electric motor (to the extent of the power of the motor, discounted with an `recuperation efficiency` currently set at 72%).

To calculate the tank-to-wheel energy, the following parameters are needed:

* the ``driving mass`` of the vehicle
* its ``rolling resistance coefficient``
* its ``aerodynamic drag coefficient``
* its ``frontal area``
* its tank-to-wheel efficiency (``TtW efficiency``)
* its ``recuperation efficiency``
* and the power of its electric motor, if any (``electric power``)

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/motive_energy.png
    :width: 900
    :alt: Calculation of the motive energy

Here is plotted the second-by-second power requirement for a large-sized battery electric vehicle, along the WLTC driving cycle:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/kw_bev_wltc.png
    :width: 900
    :alt: Calculation of the motive energy


In parallel, the ``TtW efficiency`` is calculated as the product of the following inefficiency parameters:

* ``battery discharge efficiency``
* ``fuel cell system efficiency``
* ``drivetrain efficiency``
* ``engine efficiency``

It represents the loss of energy between the energy storage and the wheels.

The power required for each second of the driving cycle is therefore summed up, and divided by the ``TtW efficiency``,
to obtain the amount of kilojoules needed in the tank (or battery) per km.

Finally, the `auxillary` energy, that is the energy needed to operate onboard equipment, is also calculated.
The sum of the `motive` and the `auxillary` energy gives the tank-to-wheel energy (``TtW energy``) of the vehicle.

Parameters such as ``battery discharge efficiency``, ``fuel cell system efficiency``, ``drivetrain efficiency``,
``engine efficiency`` and therefore, indirectly, ``TtW efficiency``, have been calibrated to obtain ``TtW energy``
figures that fit what is observed in reality.

For 2010 and 2017 vehicles, the tank-to-wheel energy use (``TtW energy``) and underlying parameters have been calibrated
against the database from the `Monitoring of CO2 emissions from passenger cars <https://www.eea.europa.eu/data-and-maps/data/co2-cars-emission-16>`_
program from the European Environment Agency. This database lists energy and emission measurement for each new passenger
car registered in the European Union, based on the NEDC and WLTC driving cycles.

Tank-to-wheel energy calibration for 2010 vehicles

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/EU_energy_comparison_2010.png
    :width: 900
    :alt: Tank-to-wheel energy calibration for 2010 vehicles


Tank-to-wheel energy calibration for 2017 vehicles

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/EU_energy_comparison.png
    :width: 900
    :alt: Tank-to-wheel energy calibration for 2017 vehicles

For the year 2000, such energy and emission measurement data was not available. Hence, we relied on the `International
Council on Clean Transportation data <https://theicct.org/chart-library-passenger-vehicle-fuel-economy>`_ that provides
historical time series on the measured fuel efficiency of diesel and petrol engines based on the WLTC driving cycle,
including its evolution between 2000 and 2010 (-20%). Therefore, the underlying parameters of ``TtW efficiency`` have
been adjusted to produce ``TtW energy`` figures about 20% more important than those observed in 2010.

Here is a comparison of the ``TtW energy`` based on the WLTC driving cycle for 2000, 2010 and 2017 vehicles:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/EU_energy_comparison_2000.png
    :width: 900
    :alt: Tank-to-wheel energy calibration for 2000 vehicles

Knowing the tank-to-wheel energy requirement allows to calculate the range (in km) of a vehicle on a full tank since:

    ``range`` [km] = (``fuel mass`` [kg] x ``LHV fuel MJ per kg`` [Mj/kg] x 1000) / ``TtW energy``

In the case of battery electric vehicles and hybrid vehicles, things are similar:

    ``range`` [km] = (``electric energy stored`` [kWh] x ``battery DoD`` [%] x 3.6 x 1000) / ``TtW energy``



How can I override the tank-to-wheel efficiency?
------------------------------------------------

With **carculator online**:

In the Parameters section, search for any or all of ``battery discharge efficiency``, ``fuel cell system efficiency``, ``drivetrain efficiency``,
 ``engine efficiency`` parameters and add them for the vehicles you wish to modify. The ``TtW efficiency`` is the
 product of those. Currently, it is not possible to modify directly the parameter ``TtW efficiency``, as it would be recalculated.
In order to do so, you need to use instead the Python library *carculator* (see next section).


With **carculator**:

Yes. After having created the CarModel() object and executed the :meth:`.set_all` method, you can override the
calculated ``TtW efficiency`` value and recalculate the TtW energy with the :meth:`.calculate_ttw_energy` method.
Here is an example for a diesel car of medium size in 2020, for which we want to set the TtW efficiency at 30% (instead of 24%)::

    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    cm.array.loc[dict(parameter="TtW efficiency",
                  powertrain="ICEV-d",
                  year=2020,
                  size="Medium")] = 0.3
    cm.calculate_ttw_energy()

If I know already the fuel consumption of a vehicle, can I override it?
-------------------------------------------------------------------

With **carculator online**:

Currently, it is not possible to modify directly the parameter ``TtW energy``, as it would be recalculated.
In order to do so, you need to use instead the Python library *carculator* (see next section).


With **carculator**:

Yes. After having created the CarModel() object and executed the :meth:`.set_all` method, you can override the
calculated ``TtW energy`` value (in kilojoules). Here is an example for a diesel car of medium size in 2020::

    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    cm.array.loc[dict(parameter="TtW energy",
                  powertrain="ICEV-d",
                  year=2020,
                  size="Medium")] = 2800



Fuel-related direct emissions
*****************************

Hot pollutants emissions
************************

Components origin
*****************

Background inventory
********************

