Modeling and assumptions
========================

The modeling of passenger vehicles in the past, present and future is complex and relies on many assumptions.
With **carculator** and **carculator online**, we wish to be transparent about those: assumptions and modeling approaches should ideally be easily
critiqued and modified.

We try here to give a comprehensive list of assumptions and modeling choices, and describe how, as a user, you
can change those.

Parameters' names are indicated ``verbatim`` and are to be used in **carculator**. The can also be accessed and modified
via its online graphical user interface **carculator online**,  via the search bar in the *Car Parameters* section.

Vehicle sizing
**************
**carculator** models vehicles along four dimensions:

* their powertrain (e.g., gasoline-run internal combustion engine, battery electric vehicle, etc.),
* their size (e.g., mini, medium, large, etc.),
* their year of production (2000, 2010, 2017 and 2040)
* and a parameter dimension (i.e., input and calculated parameters).

When **carculator** sizes the vehicles for the different powertrains, sizes and years, it starts with the
input parameter's value for the ``glider base mass``, which is essentially an initial guess for the mass of the vehicle's
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

* ``glider base mass``:the initial mass of the glider
* ``power to mass ratio``: the power-to-mass ratio
* ``combustion power share``: how much of the power is provided by an internal combustion engine
* ``combustion mass per power``: the mass of the combustion engine per unit of power

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

With **carculator online**:

Currently, it is not possible to modify directly the calculated parameter ``curb mass``, as it would be recalculated.
In order to do so, you need to use instead the Python library **carculator** (see next section). You can however
modify any of the input parameters ``glider base mass``, ``power to mass ratio``, ``combustion power share``
and ``combustion mass per power`` used to calculate ``curb mass``.
To do so, type their name in the search field of the Parameters section.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/power_to_mass_change.png
    :width: 900
    :alt: Change parameters affecting the curb mass


With **carculator**:

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

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/combustion_power_share.png
    :width: 900
    :alt: Change combustion power share parameter

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

How can I modify the battery capacity of a battery electric car?
-----------------------------------------------------------------------

Two parameters are of importance, ``energy battery mass`` [kg] and ``battery cell energy density`` [kWh/kg], so that:

``battery cell mass`` [kg] = ``energy battery mass`` [kg] × ``battery cell mass share`` [%]

``energy stored`` [kWh] = ``battery cell energy density`` [kWh/kg] x ``battery cell mass`` [kg]

Hence, by modifying either of them (or both), you can affect the capacity of the battery for a given size class.

With **carculator online**:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/battery_capacity_change.png
    :width: 900
    :alt: Change battery capacity


With **carculator**:

You can simply override the default values in ``array`` before passing it to CarModel()::

    dict_param = {('Energy Storage',  'BEV', 'Medium', 'energy battery mass', 'none'): {
                                                                                        (2000, 'loc'): 100,
                                                                                        (2010, 'loc'): 150,
                                                                                        (2017, 'loc'): 180,
                                                                                        (2040, 'loc'): 200}
                                                                                        },
                 ('Energy Storage',  'BEV', 'Medium', 'battery cell energy density', 'none'): {
                                                                                        (2000, 'loc'): 0.05,
                                                                                        (2010, 'loc'): 0.1,
                                                                                        (2017, 'loc'): 0.2,
                                                                                        (2040, 'loc'): 0.3}
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

What happens when I inter-/extrapolate to other years?
------------------------------------------------------

If the default years of 2000, 2010, 2017 and 2040 are of no interest, it is possible to inter-/extrapolate the vehicle
models to any year between 2000 and 2050. When such inter-/extrapolation is done, all the *physical* input parameters' values
are inter-/extrapolated **linearly**.

With **carculator online**:

In the Scope section, simply drag the desired years from the left frame to the right frame.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/select_vehicle_tutorial.gif
    :width: 900
    :alt: Select years

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

The `tank-to-wheel` energy consumption is the sum of:

* the `motive energy` needed to move the vehicle over 1 km
* the `auxilliary` energy needed to operate onboard equipment as well as to provide heating and cooling over 1 km

Motive energy
-------------

Once the vehicle and its powertrain has been sized, it is possible to calculate the `motive energy` required along
a specific driving cycle to overcome the following forces:

* rolling resistance
* aerodynamic resistance
* air resistance
* road gradient resistance (if provided)

on top of the *kinetic energy* needed to move the vehicle.

To calculate the motive energy, the following parameters are needed:

* the ``driving mass`` of the vehicle
* its ``rolling resistance coefficient``
* its ``aerodynamic drag coefficient``
* its ``frontal area``
* its tank-to-wheel efficiency (``TtW efficiency``)
* its ``recuperation efficiency``
* and the power of its electric motor, if any (``electric power``)

To that amount of energy is subtracted the *energy recuperated* during braking, if the vehicle is equipped with
an electric motor (to the extent of the power of the motor, discounted with a ``recuperation efficiency``).

* ``recuperation efficiency`` [%] = ``drivetrain efficiency`` [%] x ``battery charge efficiency`` [%]


.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/motive_energy.png
    :width: 900
    :alt: Calculation of the motive energy

Also, ``distance``, ``velocity`` and ``acceleration`` are derived from the driving cycle.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/driving_cycle.png
    :width: 400
    :alt: Driving cycle


In parallel, the ``TtW efficiency`` (the loss of energy between the energy storage and the wheels) is calculated as the product of the following efficiency parameters:

* ``battery discharge efficiency``
* ``fuel cell system efficiency``
* ``drivetrain efficiency``
* ``engine efficiency``

The `motive energy` is calculated as the sum of:

* rolling resistance [kg.m.s^-2] = ``driving mass`` [kg] x ``rolling resistance coefficient`` [%] x 9.81 [m/s^2]
* air resistance [kg.m.s^-2] = ``velocity`` ^2 [m^2/s^2] x (``frontal area`` [m^2] x ``aerodynamic drag coefficient`` [%] x air density [kg/m^3] / 2)
* road gradient resistance [kg.m.s^-2] = ``driving mass`` [kg] x 9.81 [m/s^2] x sin(gradient)
* kinetic force [kg.m.s^-2] = ``acceleration`` [m/s^2] x ``driving mass`` [kg]

This gives:

* force required [kg.m.s^-2] = rolling resistance + air resistance + road gradient resistance + kinetic force

Then, the gross power required is calculated as:

* power [W or kg.m^2.s^-3] = force required [kg.m.s^-2] x velocity [m/s]

The recuperated power, via electro-braking is calculated as the decelerating power (when power is negative) comprised
between 0 and the electric engine power *-1, times the recuperation efficiency:

* recuperated power [W] = power [W] * recuperation efficiency [%], when power between (-1 x electric engine power [W]) and 0

Finally, to obtain the `motive energy` the gross power minus the recuperated power (which is negative!) are summed along the driving cycle duration:

 * `motive energy` [joules] = sum ((power [W or joules/s] + recuperated power [W or joules/s]) / distance [m] / ``TtW efficiency`` [%] / 1000 [j/kj])

The `motive energy` is divided by the ``TtW efficiency`` to obtain the amount of kilojoules needed in the tank (or battery) to move the vehicle over 1 km.

Here is plotted the second-by-second gross power requirement for a large-sized battery electric vehicle, along the WLTC driving cycle:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/kw_bev_wltc.png
    :width: 900
    :alt: Calculation of the motive energy

How can I add a road gradient?
------------------------------

By default, the vehicles are compared based on a driving cycle on a flat road.



Auxilliary energy
----------------

The `auxilliary` energy, that is the energy needed to operate onboard equipment and heating and cooling systems, is also calculated
as the sum of the power demand over time.

This power demand entails:

* the average power demand for heating
* the average power demand for cooling
* the average power demand for onboard electronics

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/aux_energy.png
    :width: 900
    :alt: Auxilliary energy

This power demand is modeled as:

* ``auxilliary power demand`` [W] = ``auxilliary power base demand`` [W] + (``heating thermal demand`` [W] x ``heating energy consumption`` [0-1]) + (``cooling thermal demand`` [W] x ``cooling energy consumption`` [0-1])

``auxilliary power demand`` is summed over the driving time defined by the driving cycle and divided by the ``engine efficiency``.


The power demand for heating varies between 200 Watts and 350 Watts depending on the car size.
The power demand for cooling varies between 200 Watts and 350 Watts depending on the car size.

Note that, unlike battery electric vehicles, internal combustion engine vehicles satisfy the power demand in heating
without the additional use of energy, because ``heating energy consumption`` = 0.


Tank-to-wheel energy
--------------------

The sum of the `motive` and the `auxilliary` energy gives the tank-to-wheel energy (``TtW energy``) of the vehicle.

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

The following lower heating values (LHV) for the liquid and gaseous fuels, in Mj/kg, are used:

* conventional gasoline: 42.4
* conventioanl diesel: 42.8
* compressed natural gas: 55.5
* hydrogen: 120

Those can be changed by modifying the value of the ``LHV fuel MJ per kg`` in ``array`` before passing it to ``CarModel``.
For example, we can decrease the LHV of diesel::

    dict_param = {('Powertrain',  'ICEV-d', 'all', 'LHV fuel MJ per kg', 'none'): {
                                                                                        (2000, 'loc'): 44,
                                                                                        (2010, 'loc'): 44,
                                                                                        (2017, 'loc'): 44,
                                                                                        (2040, 'loc'): 44
                                                                                        }
    modify_xarray_from_custom_parameters(dict_param, array)


How can I override the tank-to-wheel efficiency?
------------------------------------------------

With **carculator online**:

You cannot directly override ``TtW efficieny``.
However, you can adjust any of the four parameters affecting ``TtW efficiency`` in the Tank-to-wheel efficiency section.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/ttw_efficiency_change.png
    :width: 900
    :alt: Tank-to-wheel efficiency adjustment

With **carculator**:

After having created the CarModel() object and executed the :meth:`.set_all` method, you can override the
calculated ``TtW efficiency`` value and recalculate ``TtW energy`` with the :meth:`.calculate_ttw_energy` method.
Here is an example for a diesel car of medium size in 2020, for which we want to set the TtW efficiency at 30% (instead of 24%)::

    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    cm.array.loc[dict(parameter="TtW efficiency",
                  powertrain="ICEV-d",
                  year=2020,
                  size="Medium")] = 0.3
    cm.calculate_ttw_energy()

You can also adjust any of the input parameters that affect ``TtW efficiency``, namely ``battery discharge efficiency``
 (for battery electric cars only), ``fuel cell stack efficiency`` (for fuel cell cars only), ``engine efficiency`` and
 ``drivetrain efficiency``.

If I know already the fuel consumption of a vehicle, can I override it?
-----------------------------------------------------------------------

With **carculator online**:

Currently, it is not possible to modify directly the parameter ``TtW energy``, as it would be recalculated.
In order to do so, you need to use instead the Python library *carculator* (see next section):

With **carculator**:

Yes. After having created the CarModel() object and executed the :meth:`.set_all` method, you can override the
calculated ``TtW energy`` value (in kilojoules). Here is an example for a diesel car of medium size in 2020::

    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    cm.array.loc[dict(parameter="TtW energy",
                  powertrain="ICEV-d",
                  year=2020,
                  size="Medium")] = 2800


Fuel blends
***********

The user can define fuel blends. The following fuel types are available, along with their lower heating value (in MJ/kg)
and CO2 emission factor (kg CO2/kg fuel.

Hydrogen technologies (LHV: 120 MJ/kg)
......................................

* 'electrolysis'
* 'smr - natural gas'
* 'smr - natural gas with CCS'
* 'smr - biogas'
* 'smr - biogas with CCS'
* 'coal gasification'
* 'wood gasification'
* 'wood gasification with CCS'

Natural gas technologies
------------------------

* 'cng' (55.5 MJ/kg, 2.65 kg CO2/kg CNG)
* 'biogas' (55.5 MJ/kg, 2.65 kg CO2/kg CNG)
* 'syngas' (55.5 MJ/kg, 2.65 kg CO2/kg CNG)

Diesel technologies
-------------------

* 'diesel' (42.8 MJ/kg, 3.14 kg CO2/kg)
* 'biodiesel - algae' (31.7 MJ/kg, 2.85 kg CO2/kg)
* 'biodiesel - cooking oil' (31.7 MJ/kg, 2.85 kg CO2/kg)
* 'synthetic diesel' (43.3 MJ/kg, 3.16 kg CO2/kg)

Petrol technologies
-------------------

* 'petrol' (42.4 MJ/kg, 3.18 kg CO2/kg)
* 'bioethanol - wheat straw' (26.8 MJ/kg, 1.91 kg CO2/kg)
* 'bioethanol - maize starch' (26.8 MJ/kg, 1.91 kg CO2/kg)
* 'bioethanol - sugarbeet' (26.8 MJ/kg, 1.91 kg CO2/kg)
* 'bioethanol - forest residues' (26.8 MJ/kg, 1.91 kg CO2/kg)
* 'synthetic gasoline' (42.4 MJ/kg, 3.18 kg CO2/kg)

Once the fuel blend is defined, the range is calculated once again, now considering the new energy amount stored in the tank.
Therefore, a car solely running on bio-ethanol will have a reduced range, increasing the fuel consumption and emissions
related to the growing of crops and supply of fuel. The tailpipe CO2 emissions may not necessarily increase as biofuels
have generally lower CO2 emission factors.

It is important to note that CO2 emissions of biogenic origin from biofuels are characterized with a similar Global Warming Potential factor
as those for conventional fossil fuels. However, CO2 uptake is considered during biomass growth.

Fuel-related direct emissions
*****************************

Carbon dioxide emissions from fuel combustion are calculated based on the fuel blend defined by the user (see above).

    carbon dioxide emission [kg/km] = CO2_fuel x ``fuel mass`` [kg] x share_fuel / ``range`` [km]

This is calculated for every fuel type found in the blend (primary and secondary fuel).

Other emissions based on fuel combustion are considered, from Spielmann et al., Transport Services Data v.2 (2007).
However those only apply when conventional diesel or conventional gasoline is burnt:

* Cadmium
* Chromium and Chromium VI
* Copper
* Nickel
* Selenium
* Zinc


Hot pollutants emissions
************************

**carculator** quantifies the emissions of the following substances:

* Hydrocarbons
* Carbon monoxide
* Nitrogen oxides
* Particulate matters
* Methane
* NMVOC
* Lead
* Sulfur dioxide
* Dinitrogen oxide
* Ammonia
* Benzene

It does so by correlating the emission of a substance at a given speed and the speed given for each second of the driving cycle.

The emission of substances function of the speed level is sourced from the
`Handbook Emission Factors for Road Transport <https://www.hbefa.net/e/index.html>`_ for vehicles of various emission
standards (from Euro-0 to Euro-6d).

Here is such correlation plotted for gasoline-run vehicles with a Euro-6d emission standard:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/hbefa_petrol_euro6d.png
    :width: 900
    :alt: Substance emission versus speed, petrol, Euro-6d

Given the years selected, the corresponding emission factors are chosen:

* before 1993: Euro-0
* between 1993 and 1997: Euro-1
* between 1998 and 2000: Euro-2
* between 2001 and 2005: Euro-3
* between 2006 and 2010: Euro-4
* between 2011 and 2014: Euro-5
* above 2015: Euro-6

Emissions are summed over the duration of the driving cycle. Furthermore, some driving cycles have distinct parts
corresponding to different driving environments: urban, suburban, highway, etc. These driving environments are used
to further split emissions and be more precise on the fate of the substances and the exposure of the population.

Noise emissions
***************

Given the driving cycle, where speed [km/h] is given along time [s], noise levels (in dB) are calculated for each of the
8 octaves (or frequency ranges) to obtain `propulsion` and `rolling noise` levels, based on the
`CNOSSOS model <https://ec.europa.eu/jrc/en/publication/reference-reports/common-noise-assessment-methods-europe-cnossos-eu>`_.

For electric engines, `special coefficients apply <https://hal.archives-ouvertes.fr/hal-01355872/document>`_.

Also, electric cars are added a warning signal of 56 dB at speed levels lower than 20 km/h.
Hybrid cars are assumed to use an electric engine up to a speed level of 30 km/h, beyond which the combustion engine is used.
The sum of the propulsion and rolling noise levels is converted to noise power (in joules) and divided by the distance
driven to obtain the noise power par km driven (joules/km), for each octave.

Noise emissions are further compartmented into urban, sub-urban and rural geographical environments based on speed
intervals given by the driving cycle.
The study from  `Cucurachi et al. 2014 <https://www.ncbi.nlm.nih.gov/pubmed/24035845>`_ is used to characterize noise
emissions against midpoint and endpoint indicators, expressed in Person-Pascal-second and DALYs, respectively.

Overall, propulsion noise emissions dominate in urban environments, thereby justifying the use of electric cars in that
regard. In sub-urban and rural environments, rolling noise emissions dominate above a speed level around 50 km/h.

It is important to note that although **carculator** differentiates noise coefficients by powertrain
(internal combustion engine, electric and hybrid), it is not possible to differentiate them by size class.
Therefore, the noise produced by a `small` vehicle will be similar to that produced by a `large` vehicle.


Vehicle inventory
*****************
This section presents the vehicle inventory once its size, mass, energy consumption and emissions are known.

.. csv-table:: Vehicle inventory
    :file: table_1.csv
    :widths: 10 10 30 30 10 10
    :header-rows: 1


Fuel pathways
*************

Different fuel pathways can be selected for a given powertrain type.
The table below lists them.

.. csv-table:: Fuel pathways
    :file: table_2.csv
    :widths: 15 25 30 30
    :header-rows: 1

Electricity mixes for battery charging and hydrogen production
**************************************************************

**carculator** has national electricity mixes for more than 80 countries, gathered from the following sources:

* European Union State members and the UK: `EU Reference Scenario 2016 <https://ec.europa.eu/energy/en/data-analysis/energy-modelling/eu-reference-scenario-2016>`_
* Switzerland: STEM model - Panos E, Kober T, Wokaun A. Long term evaluation of electric storage technologies vs alternative flexibility options for the Swiss energy system. Appl Energy 2019;252:113470
* African countries: `TEMBA <http://www.osemosys.org/temba.html>`_ model
* Other countries: `IEA World Energy outlook 2017 <https://www.iea.org/reports/world-energy-outlook-2017>`_

Unless a specific electricity mix is indicated by the user, such national mixes are used when modeling the energy chain
for battery and fuel cell electric vehicles (BEV, FCEV), for battery charging and the production of hydrogen via electrolysis, respectively.

Knowing the production year of the vehicle, considered to be its first year of use, as well as its annual mileage,
the number of years of use is calculated. Hence, **the electricity mix used is the kilometer-distributed mix over the
years of use of the vehicle**.

If the annual mileage of the vehicle is evenly distributed throughout its lifetime, the electricity mix used therefore
equals the average of the year-by-year national mixes comprised between Year 0 and Year 0 + the number of years of use.


Background inventory
********************

Besides datasets adapted from the literature, vehicle inventories also rely on a number of datasets provided by the database `ecoinvent cutoff 3.6 <https://www.ecoinvent.org>`_
such as "market for glider, passenger car", "market for diesel", etc.

However **carculator** does not directly use the database as is: the database and its datasets are modified according to
projections provided by the Integrated Assessment Model `REMIND <https://www.pik-potsdam.de/research/transformation-pathways/models/remind/remind>`_.

REMIND provides projections for different regions in the world until 2150, following different energy scenarios,
described `here <https://github.com/romainsacchi/rmnd-lca/blob/master/rmnd_lca/data/remind_output_files/description.md>`_.

Projection outputs include the expected change over time in efficiency for power plants, steel making, cement production, etc.

Using the Python library `rmnd_lca <https://github.com/romainsacchi/rmnd-lca/tree/master/rmnd_lca>`_, we produce a number
of ecoinvent databases with the inclusion of REMIND projections, so that future improvements in electricity production, among others,
propagate into the datasets involved in the vehicles' inventories.

**carculator** comes with pre-calculated impact values for ecoinvent datasets from the following databases:

* 2005 - ecoinvent-REMIND, SSP2-Base
* 2010 - ecoinvent-REMIND, SSP2-Base
* 2020 - ecoinvent-REMIND, SSP2-Base
* 2030 - ecoinvent-REMIND, SSP2-Base
* 2040 - ecoinvent-REMIND, SSP2-Base
* 2050 - ecoinvent-REMIND, SSP2-Base
* 2005 - ecoinvent-REMIND, SSP2-PkBudg1100
* 2010 - ecoinvent-REMIND, SSP2-PkBudg1100
* 2020 - ecoinvent-REMIND, SSP2-PkBudg1100
* 2030 - ecoinvent-REMIND, SSP2-PkBudg1100
* 2040 - ecoinvent-REMIND, SSP2-PkBudg1100
* 2050 - ecoinvent-REMIND, SSP2-PkBudg1100

Depending on the year of analysis and the energy scenario demanded, **carculator** picks the corresponding datasets.
If year of analysis in between the available years is demanded, a linear interpolation is used.

With **carculator online**, the results provided only use the "SSP2-Base" energy scenario of REMIND, projecting a global
atmospheric temperature increase by 3.5 degrees Celsius by 2100.
