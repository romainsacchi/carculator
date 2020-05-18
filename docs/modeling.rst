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

A second step consists into calculating the mass of the combustion and electric engine, based on the following:

.. math::

    power demand (``power``) [kW] = ``power-to-mass ratio`` [kW/kg] x ``curb mass`` [kg]
    electrical power demand (``electric power``) [kW] = power demand (``power``) [kW] x (1 - ``combustion power share`` [%])
    ``electric engine mass`` [kW] = (``electric power`` [kW] x ``electric mass per power`` [kg/kW]) + ``electric fixed mass`` [kg]
    combustion power demand (``combustion power``) [kW] = ``power`` [kW] x ``combustion power share`` [%]
    ``combustion engine mass`` [kW] = (``combustion power`` [kW] x ``combustion mass per power`` [kg/kW]) + ``combustion fixed mass`` [kg]

As well as for the mass of the powertrain:

.. math::

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

With ***carculator online***:
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

Fuel-related direct emissions
*****************************

Hot pollutants emissions
************************

Components origin
*****************

Background inventory
********************

**carculator** is a parameterized model that allows to generate and characterize life cycle inventories for different
vehicle configurations, according to selected:

* powertrain technologies (9): petrol engine, diesel engine, electric motor, hybrid, plugin-hybrid, etc.,
* year of operation (2): 2000, 2010, 2017, 2040 (with the possibility to interpolate in between, and up to 2050)
* and sizes (7): Mini, Large, etc.

The methodology used to develop `carculator` is explained in:
carculator: an open-source tool for prospective environmental and economic life cycle assessment of vehicles. When, Where and How can battery-electric vehicles help reduce greenhouse gas emissions?
Romain Sacchi, Christian Bauer, Brian Cox, Christopher Mutel
Environmental Modelling and Software (2020, submitted)

At the moment, the tool has a focus on passenger cars.

It is initially based on the model developed in `Uncertain environmental footprint of current and future battery electric
vehicles by Cox, et al (2018) <https://pubs.acs.org/doi/10.1021/acs.est.8b00261>`_.

More specifically, **carculator** generates `Brightway2 <https://brightwaylca.org/>`_ inventories, but also directly provides characterized
results against several midpoint indicators from the impact assessment method ReCiPe as well as life cycle cost indicators.

**carculator** is a special in the way that it uses time- and energy-scenario-differentiated background inventories for the future,
resulting from the coupling between the `ecoinvent 3.6 database <https://ecoinvent.org>`_ and the scenario outputs of PIK's
integrated assessment model `REMIND <https://www.pik-potsdam.de/research/transformation-pathways/models/remind/remind>`_.
This allows to perform prospective study while consider future expected changes in regard to the production of electricity,
cement, steel, heat, etc.

Objective
---------

The objective is to produce life cycle inventories for vehicles in a transparent, comprehensive and quick manner,
to be further used in prospective LCA of transportation technologies.

Why?
----

Many life cycle assessment (LCA) models of passenger cars exist. Yet, because LCA of vehicles, particularly for electric battery vehicles,
are sensitive to assumptions made in regards to electricity mix used for charging, lifetime of the battery, etc., it has led
to mixed conclusions being published in the scientific literature. Because the underlying calculations are kept undocumented,
it is not always possible to explain the disparity in the results given by these models, which can contribute to adding confusion among the public.

Because **carculator** is kept **as open as possible**, the methods and assumptions behind the generation of results are
easily identifiable and adjustable.
Also, there is an effort to keep the different modules (classes) separated, so that improving certain areas of the model is relatively
easy and does not require changing extensive parts of the code. In that regard, contributions are welcome.

Finally, beside being more flexible and transparent, **carculator** provides interesting features, such as:

* a stochastic mode, that allows fast Monte Carlo analyses, to include uncertainty at the vehicle level
* possibility to override any or all of the 200+ default input car parameters (e.g., number of passengers, drag coefficient) but also calculated parameters (e.g., driving mass).
* hot pollutants emissions as a function of the driving cycle, using `HBEFA <https://www.hbefa.net/e/index.html>`_ 4.1 data, further divided between rural, suburban and urban areas
* noise emissions, based on `CNOSSOS-EU <https://ec.europa.eu/jrc/en/publication/reference-reports/common-noise-assessment-methods-europe-cnossos-eu>`_ models for noise emissions and `Noise footprint from personal land‐based mobility by Cucurachi, et al (2019) <https://onlinelibrary.wiley.com/doi/full/10.1111/jiec.12837>`_ for inventory modelling and mid- and endpoint characterization of noise emissions, function of driving cycle and further divided between rural, suburban and urban areas
* export of inventories as an Excel file, to be used with Brightway2 or Simapro (in progress), including uncertainty information. This requires the user to have `ecoinvent 3.6 cutoff` installed on the LCA software the car inventories are exported to.
* export inventories directly into Brightway2, as a LCIImporter object to be registered. Additionally, when run in stochastic mode, it is possible to export arrays of pre-sampled values using the `presamples <https://pypi.org/project/presamples/>`_ library to be used together with the Monte Carlo function of Brightway2.
* development of an online graphical user interface: `carculator online <https://carculator.psi.ch>`_

How to install this package?
----------------------------

**carculator** is a Python package, and is primarily to be used from within a Python 3.x environment.
Because **carculator** is still at an early development stage, it is a good idea to install it in a separate environment,
such as a conda environment::

    conda create -n <name of the environment> python=3.7

Once your environment created, you should activate it::

    conda activate <name of the environment>

And install the **carculator** library in your new environment via Conda::

    pip install carculator

This will install the package and the required dependencies.

How to use it?
--------------

Static vs. Stochastic mode
**************************

Note: many examples are given in this `notebook <https://github.com/romainsacchi/carculator/blob/master/examples/Examples.ipynb>`_ that you can run directly on your computer..

The inventories can be calculated using the most likely value of the given input parameters ("static" mode), but also using
randomly-generated values based on a probability distribution for those ("stochastic" mode).

For example, the drivetrain efficiency of SUVs in 2017, regardless of the powertrain, is given the most likely value (i.e., the mode) of 0.38,
but with a triangular probability distribution with a minimum and maximum of 0.3 and 0.4, respectively.

Creating car models in static mode will use the most likely value of the given parameters to dimension the cars, etc., such as:

.. code-block:: python

   from carculator import *
   cip = CarInputParameters()
   cip.static()
   dcts, array = fill_xarray_from_input_parameters(cip)
   cm = CarModel(array)
   cm.set_all()


Alternatively, if one wishes to work with probability distributions as parameter values instead:

.. code-block:: python

    from carculator import *
    cip = CarInputParameters()
    cip.stochastic(800)
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array)
    cm.set_all()


This effectively creates 800 iterations of the same car models, picking pseudo-random value for the given parameters,
within the probability distributions defined. This allows to assess later the effect of uncertainty propagation on
characterized results.

In both case, a CarModel object is returned, with a 4-dimensional array `array` to store the generated parameters values, with the following dimensions:

0. Vehicle sizes (called "size"):
    * Mini
    * Small
    * Lower medium
    * Medium
    * Large
    * SUV
    * Van

1. Powertrains:
    * ICEV-p, ICEV-d, ICEV-g: vehicles with internal combustion engines running on gasoline, diesel and compressed gas, respectively.
    * HEV-p, HEV-d: vehicles with internal combustion engines running on gasoline and diesel, assisted with an electric engine.
    * PHEV-p, PHEV-d: vehicles with internal combustion engines running on gasoline and diesel, assisted with a plugin electric engine.
    * BEV: battery electric vehicles.
    * FCEV: fuel cell electric vehicles.

2. Year. Anything between 2000 and 2050.

3. Iteration number (length = 1 if static(), otherwise length = number of iterations).


:meth:`cm.set_all()` generates a CarModel object and calculates the energy consumption, components mass, as well as
exhaust and non-exhaust emissions for all vehicle profiles.

Custom values for given parameters
**********************************

You can pass your own values for the given parameters, effectively overriding the default values.

For example, you may think that the *base mass of the glider* for large diesel and petrol cars is 1600 kg in 2017
and 1,500 kg in 2040, and not 1,500 kg as defined by the default values. It is easy to change this value.
You need to create first a dictionary and define your new values as well as a probability distribution if needed :

.. code-block:: python

    dic_param = {
    ('Glider', ['ICEV-d', 'ICEV-p'], 'Large', 'glider base mass', 'triangular'): {(2017, 'loc'): 1600.0,
                                                                 (2017, 'minimum'): 1500.0,
                                                                 (2017, 'maximum'): 2000.0,
                                                                 (2040, 'loc'): 1500.0,
                                                                 (2040, 'minimum'): 1300.0,
                                                                 (2040, 'maximum'): 1700.0}}

Then, you simply pass this dictionary to `modify_xarray_from_custom_parameters(<dic_param or filepath>, array)`, like so:

.. code-block:: python

    cip = CarInputParameters()
    cip.static()
    dcts, array = fill_xarray_from_input_parameters(cip)
    modify_xarray_from_custom_parameters(dic_param, array)
    cm = CarModel(array, cycle='WLTC')
    cm.set_all()

Alternatively, instead of a Python dictionary, you can pass a file path pointing to an Excel spreadsheet that contains
the values to change, following `this template <https://github.com/romainsacchi/carculator/raw/master/docs/template_workbook.xlsx>`_.

The following probability distributions are accepted:
* "triangular"
* "lognormal"
* "normal"
* "uniform"
* "none"

Inter and extrapolation of parameters
*************************************

**carculator** creates by default car models for the year 2000, 2010, 2017 and 2040.
It is possible to inter and extrapolate all the parameters to other years simply by writing:

.. code-block:: python

    array = array.interp(year=[2018, 2022, 2035, 2040, 2045, 2050],  kwargs={'fill_value': 'extrapolate'})

However, we do not recommend extrapolating for years before 2000 or beyond 2050.

Changing the driving cycle
**************************

**carculator** gives the user the possibility to choose between several driving cycles. Driving cycles are determinant in
many aspects of the car model: hot pollutant emissions, noise emissions, tank-to-wheel energy, etc. Hence, each driving
cycle leads to slightly different results. By default, if no driving cycle is specified, the WLTC driving cycle is used.
To specify a driving cycle, simply do:

.. code-block:: python

    cip = CarInputParameters()
    cip.static()
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array, cycle='WLTC 3.4')
    cm.set_all()

In this case, the driving cycle *WLTC 3.4* is chosen (this driving cycle is in fact a sub-part of the WLTC driving cycle,
mostly concerned with driving on the motorway at speeds above 80 km/h). Driving cycles currently available:

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

The user can also create custom driving cycles and pass it to the :class:`CarModel` class:

.. code-block:: python

    import numpy as np
    x = np.linspace(1, 1000)
    def f(x):
        return np.sin(x) + np.random.normal(scale=20, size=len(x)) + 70

    cycle = f(x)
    cm = CarModel(array, cycle=cycle)

Accessing calculated parameters of the car model
************************************************
Hence, the tank-to-wheel energy requirement per km driven per powertrain technology for a SUV in 2017 can be obtained
from the CarModel object:

.. code-block:: python

    TtW_energy = cm.array.sel(size='SUV', year=2017, parameter='TtW energy', value=0) * 1/3600 * 100

    plt.bar(TtW_energy.powertrain, TtW_energy)
    plt.ylabel('kWh/100 km')
    plt.show()

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/fig_kwh_100km.png
    :width: 400
    :alt: Alternative text

Note that if you call the :meth:`stochastic` method of the :class:`CarInputParameters`, you would have several values stored for a given calculated parameter
in the array. The number of values correspond to the number of iterations you passed to :meth:`stochastic`.

For example, if you ran the model in stochastic mode with 800 iterations as shown in the section above, instead of one
value for the tank-to-wheel energy, you would have a distribution of values:

.. code-block:: python

    l_powertrains = TtW_energy.powertrain
    [plt.hist(e, bins=50, alpha=.8, label=e.powertrain.values) for e in TtW_energy]
    plt.ylabel('kWh/100 km')
    plt.legend()

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/stochastic_example_ttw.png
    :width: 400
    :alt: Alternative text

Any other attributes of the CarModel class can be obtained in a similar way.
Hence, the following code lists all direct exhaust emissions included in the inventory of an petrol Van in 2017:

List of all the given and calculated parameters of the car model:

.. code-block:: python

    list_param = cm.array.coords['parameter'].values.tolist()

Return the parameters concerned with direct exhaust emissions (we remove noise emissions):

.. code-block:: python

    direct_emissions = [x for x in list_param if 'emission' in x and 'noise' not in x]

Finally, return their values and display the first 10 in a table:

.. code-block:: python

    cm.array.sel(parameter=direct_emissions, year=2017, size='Van', powertrain='BEV').to_dataframe(name='direct emissions')



Or we could be interested in visualizing the distribution of non-characterized noise emissions, in joules:

.. code-block:: python

    noise_emissions = [x for x in list_param if 'noise' in x]
    data = cm.array.sel(parameter=noise_emissions, year=2017, size='Van', powertrain='ICEV-p', value=0)\
        .to_dataframe(name='noise emissions')['noise emissions']
    data[data>0].plot(kind='bar')
    plt.ylabel('joules per km')

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/example_noise_emissions.png
    :width: 400
    :alt: Alternative text

Modify calculated parameters
****************************

As input parameters, calculated parameters can also be overridden. For example here, we override the `driving mass`
of large diesel vehicles for 2010 and 2017:

.. code-block:: python

    cm.array.loc['Large','ICEV-d', 'driving mass', [2010, 2017]] = [[2000],[2200]]

Characterization of inventories (static)
****************************************

**carculator** makes the characterization of inventories easy. You can characterize the inventories directly from
**carculator** against midpoint impact assessment methods.

For example, to obtain characterized results against the midpoint impact assessment method ReCiPe for all cars:

.. code-block:: python

    ic = InventoryCalculation(cm.array)
    results = ic.calculate_impacts()


Hence, to plot the carbon footprint for all medium cars in 2017:

.. code-block:: python

    results.sel(size='Medium', year=2017, impact_category='climate change', value=0).to_dataframe('impact').unstack(level=1)['impact'].plot(kind='bar',
                stacked=True)
    plt.ylabel('kg CO2-eq./vkm')
    plt.show()

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/example_carbon_footprint.png
    :width: 400
    :alt: Alternative text

Note that, for now, only the ReCiPe method is available for midpoint characterization. Also, once the instance of the :class:`CarModel`
class has been created, there is no need to re-create it in order to calculate additional environmental impacts (unless you wish to
change values of certain input or calculated parameters, the driving cycle or go from static to stochastic mode).

Characterization of inventories (stochastic)
********************************************

In the same manner, you can obtain distributions of results, instead of one-point values if you have run the model in
stochastic mode (with 500 iterations and the driving cycle WLTC).

.. code-block:: python

    cip = CarInputParameters()
    cip.stochastic(500)
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    scope = {
        'powertrain':['BEV', 'PHEV'],
    }
    ic = InventoryCalculation(cm.array, scope=scope)

    results = ic.calculate_impacts()

    data_MC = results.sel(impact_category='climate change').sum(axis=3).to_dataframe('climate change')
    plt.style.use('seaborn')
    data_MC.unstack(level=[0,1,2]).boxplot(showfliers=False, figsize=(20,5))
    plt.xticks(rotation=70)
    plt.ylabel('kg CO2-eq./vkm')

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/example_stochastic_BEV_PHEV.png
    :width: 400
    :alt: Alternative text

Many other examples are described in a Jupyter Notebook in the ``examples`` folder.

Export of inventories (static)
******************************

Inventories can be exported as:
    * a Python list of exchanges
    * a Brightway2 bw2io.importers.base_lci.LCIImporter object, ready to be imported in a Brigthway2 environment
    * an Excel file, to be imported in a Brigthway2 environment

.. code-block:: python

    ic = InventoryCalculation(cm.array)

    # export the inventories as a Python list
    mylist = ic.export_lci()
    # export the inventories as a Brightway2 object
    import_object = ic.export_lci_to_bw()
    # export the inventories as an Excel file (returns the file path of the created file)
    filepath = ic.export_lci_to_excel()

Export of inventories (stochastic)
**********************************

If you had run the model in stochastic mode, the export functions return in addition an array that contains pre-sampled values
for each parameter of each car, in order to perform Monte Carlo analyses in Brightway2.

.. code-block:: python

    ic = InventoryCalculation(cm.array)

    # export the inventories as a Python list
    mylist, presamples_arr = ic.export_lci()
    # export the inventories as a Brightway2 object
    import_object, presamples_arr = ic.export_lci_to_bw()
    # export the inventories as an Excel file (note that this method does not return the presamples array)
    filepath = ic.export_lci_to_excel()

Import of inventories (static)
******************************

The background inventory is originally a combination between ecoinvent 3.6 and outputs from PIK's REMIND model.
Outputs from PIK's REMIND are used to project expected progress in different sectors into ecoinvent. For example, the efficiency
of electricity-producing technologies as well as the electricity mixes in the future for the main world regions
are built upon REMIND outputs.
The library used to create hybrid versions of the ecoinvent database from PIK's REMIND is called `rmnd_lca <https://github.com/romainsacchi/rmnd-lca>`_.
This means that, as it is, the inventory cannot properly link to ecoinvent 3.6 unless some transformation is performed
before. These transformations are in fact performed by default when exporting the inventory. Hence, when doing:

.. code-block:: python

    ic.export_lci_to_excel()

the resulting inventory should properly link to ecoinvent 3.6. Should you wish to export an inventory to link with a
REMIND-modified version of ecoinvent, just export the inventory with the `ecoinvent_compatibility` argument
set to `False`.

.. code-block:: python

    ic.export_lci_to_excel(ecoinvent_compatibility=False)

In that case, the inventory will only link to a custom ecoinvent database produced by `rmnd_lca`.

But in any case, the following script should successfully import the inventory:

.. code-block:: python

    import brightway2 as bw
    bw.projects.set_current("test_carculator")
    import bw2io
    fp = r"C:\file_path_to_the_inventory\lci-test.xlsx"

    i = bw2io.ExcelImporter(fp)
    i.apply_strategies()

    if 'additional_biosphere' not in bw.databases:
        i.create_new_biosphere('additional_biosphere')
    i.match_database("name_of_the_ecoinvent_db", fields=('name', 'unit', 'location', 'reference product'))
    i.match_database("biosphere3", fields=('name', 'unit', 'categories'))
    i.match_database("additional_biosphere", fields=('name', 'unit', 'categories'))
    i.match_database(fields=('name', 'unit', 'location'))

    i.statistics()
    i.write_database()
