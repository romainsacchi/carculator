Introduction
============

``carculator`` is a parameterized model that allows to generate and characterize life cycle inventories for different
vehicle configurations, according to selected:

* powertrain technologies (8): petrol engine, diesel engine, electric motor, hybrid, etc.,
* year of operation (2): 2017, 2040 (with the possibility to interpolate in between)
* and sizes (7): Mini, Large, etc.

At the moment, the tool has a focus on passenger cars.

It is initially based on the model developed in `Uncertain environmental footprint of current and future battery electric
vehicles by Cox, et al (2018) <https://pubs.acs.org/doi/10.1021/acs.est.8b00261>`_.

More specifically, ``carculator`` generates `Brightway2 <https://brightwaylca.org/>`_ inventories, but also directly provides characterized
results against several midpoint indicators from the impact assessment method ReCiPe as well as life cycle cost indicators.

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

Because ``carculator`` is kept **as open as possible**, the methods and assumptions behind the generation of results are
easily identifiable and adjustable.
Also, there is an effort to keep the different modules (classes) separated, so that improving certain areas of the model is relatively
easy and does not require changing extensive parts of the code. In that regard, contributions are welcome.

Finally, beside being more flexible and transparent, ``carculator`` provides interesting features, such as:

* a stochastic mode, that allows fast Monte Carlo analyses, to include uncertainty at the vehicle level
* possibility to override any or all of the 200+ default input car parameters (e.g., number of passengers, drag coefficient)
but also calculated parameters (e.g., driving mass).
* hot pollutants emissions as a function of the driving cycle, using `HBEFA <https://www.hbefa.net/e/index.html>`_ 4.1 data, further divided between rural, suburban and urban areas
* noise emissions, based on `CNOSSOS-EU <https://ec.europa.eu/jrc/en/publication/reference-reports/common-noise-assessment-methods-europe-cnossos-eu>`_ models for noise emissions and `Noise footprint from personal land‐based mobility by Cucurachi, et al (2019) <https://onlinelibrary.wiley.com/doi/full/10.1111/jiec.12837>`_ for inventory modelling and mid- and endpoint characterization of noise emissions, function of driving cycle and further divided between rural, suburban and urban areas
* export of inventories as an Excel file, to be used with Brightway2 or Simapro (in progress), including uncertainty information. This requires the user to have `ecoinvent 3.6 cutoff` installed on the LCA software the car inventories are exported to.
* export inventories directly into Brightway2, as a LCIImporter object to be registered. Additionally, when run in stochastic mode, it is possible to export arrays of pre-sampled values using the `presamples <https://pypi.org/project/presamples/>`_ library to be used together with the Monte Carlo function of Brightway2.

How to install this package?
----------------------------

``carculator`` is a Python package, and is primarily to be used from within a Python 3.x environment.
Because ``carculator`` is still at an early development stage, it is a good idea to install it in a separate environment,
such as a conda environment::

    conda create -n <name of the environment> python=3.7

Once your environment created, you should activate it::

    conda activate <name of the environment>

And install the development version of the ``carculator`` package in your new environment via Conda::

    conda install -c conda-forge -c pascallesage -c cmutel -c romainsacchi/label/nightly carculator-dev

This will install the package and the required dependencies.

How to use it?
--------------

Static vs. Stochastic mode
**************************

The inventories can be calculated using the most likely value of the given input parameters, but also using
randomly-generated values based on a probability distribution for those.
For example, the drivetrain efficiency of SUVs, regardless of the powertrain, is given the most likely value of 0.38,
but with a triangular probability distribution with a minimum and maximum of 0.3 and 0.4, respectively.

Creating car models in static mode will use the most likely value of the given parameters to dimension the cars, etc., such as::

   from carculator import *
   cip = CarInputParameters()
   cip.static()
   dcts, array = fill_xarray_from_input_parameters(cip)
   cm = CarModel(array)
   cm.set_all()


Alternatively, if one wishes to work with probability distributions as parameter values instead::

    from carculator import *
    cip = CarInputParameters()
    cip.stochastic(800)
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array)
    cm.set_all()


This effectively creates 800 iterations of the same car models, picking pseudo-random value for the given parameters,
within the probability distributions defined. This allows to assess later the effect of uncertainty propagation on
characterized results.

In both case, a CarModel object is returned, with a 3-dimensional array `array` to store the generated parameters values, with the following dimensions:

0. Vehicle size, e.g. "small", "medium". str.
1. Powertrain, e.g. "ICE-p", "BEV". str.
2. Year. int.
3. Values. float


:meth:`cm.set_all()` generates a CarModel object and calculates the energy consumption, components mass, as well as
exhaust and non-exhaust emissions for all vehicle profiles.

Custom values for given parameters
**********************************

You can pass your own values for the given parameters, effectively overriding the default values.

For example, you may think that the *base mass of the glider* for large diesel and petrol cars is 1600 kg in 2017
and 1500 kg in 2040, and not 1,500 kg as defined by the default values. It is easy to change this value.
You need to create first a dictionary and define your new values as well as a probability distribution if needed ::

    dic_param = {
    ('Glider', ['ICEV-d', 'ICEV-p'], 'Large', 'glider base mass', 'triangular'): {(2017, 'loc'): 1600.0,
                                                                 (2017, 'minimum'): 1500.0,
                                                                 (2017, 'maximum'): 2000.0,
                                                                 (2040, 'loc'): 1500.0,
                                                                 (2040, 'minimum'): 1300.0,
                                                                 (2040, 'maximum'): 1700.0}}

Then, you simply pass this dictionary to `modify_xarray_from_custom_parameters(<dic_param or filepath>, array)`, like so::

    cip = CarInputParameters()
    cip.static()
    dcts, array = fill_xarray_from_input_parameters(cip)
    modify_xarray_from_custom_parameters(dic_param, array)
    cm = CarModel(array, cycle='WLTC')
    cm.set_all()

Alternatively, instead of a Python dictionary, you can pass a file path pointing to an Excel spreadsheet that contains
the values to change, following `this template <https://github.com/romainsacchi/coarse/raw/master/docs/template_workbook.xlsx>`_.

Inter and extrapolation of parameters
*************************************

It is posisble to inter and extrapolate all the parameters to other years simply by doing::

    array = array.interp(year=[2018, 2022, 2035, 2040, 2045, 2050],  kwargs={'fill_value': 'extrapolate'})


Changing the driving cycle
**************************

``carculator`` gives the user the possibility to choose between several driving cycles. Driving cycles are determinant in
many aspects of the car model: hot pollutant emissions, noise emissions, tank-to-wheel energy, etc. Hence, each driving
cycle leads to slightly different results. By default, if no driving cycle is specified, the WLTC driving cycle is used.
To specify a driving cycle, simply do::

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

The user can also create custom driving cycles and pass it to the :class:`CarModel` class::

    import numpy as np
    x = np.linspace(1, 1000)
    def f(x):
        return np.sin(x) + np.random.normal(scale=20, size=len(x)) + 70

    cycle = f(x)
    cm = CarModel(array, cycle=cycle)

Accessing calculated parameters of the car model
************************************************
Hence, the tank-to-wheel energy requirement per km driven per powertrain technology for a SUV in 2017 can be obtained
from the CarModel object::

    TtW_energy = cm.array.sel(size='SUV', year=2017, parameter='TtW energy', value=0) * 1/3600 * 100

    plt.bar(TtW_energy.powertrain, TtW_energy)
    plt.ylabel('kWh/100 km')
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/fig_kwh_100km.png
    :width: 400
    :alt: Alternative text

Note that if you call the :meth:`stochastic` method of the :class:`CarInputParameters`, you would have several values stored for a given calculated parameter
in the array. The number of values correspond to the number of iterations you passed to :meth:`stochastic`.

For example, if you ran the model in stochastic mode with 800 iterations as shown in the section above, instead of one
value for the tank-to-wheel energy, you would have a distribution of values::

    l_powertrains = TtW_energy.powertrain
    [plt.hist(e, bins=50, alpha=.8, label=e.powertrain.values) for e in TtW_energy]
    plt.ylabel('kWh/100 km')
    plt.legend()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/stochastic_example_ttw.png
    :width: 400
    :alt: Alternative text

Any other attributes of the CarModel class can be obtained in a similar way.
Hence, the following lists all direct exhaust emissions included in the inventory of an petrol Van in 2017:

List of all the given and calculated parameters of the car model::

    list_param = cm.array.coords['parameter'].values.tolist()

Return the parameters concerned with direct exhaust emissions (we remove noise emissions)::

    direct_emissions = [x for x in list_param if 'emission' in x and 'noise' not in x]

Finally, return their values and display the first 10 in a table::

    cm.array.sel(parameter=direct_emissions, year=2017, size='Van', powertrain='BEV').to_dataframe(name='direct emissions')

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_direct_emissions.png
    :width: 400
    :alt: Alternative text


Or we could be interested in visualizing the distribution of non-characterized noise emissions, in joules::

    noise_emissions = [x for x in list_param if 'noise' in x]
    data = cm.array.sel(parameter=noise_emissions, year=2017, size='Van', powertrain='ICEV-p', value=0)\
        .to_dataframe(name='noise emissions')['noise emissions']
    data[data>0].plot(kind='bar')
    plt.ylabel('joules per km')

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_noise_emissions.png
    :width: 400
    :alt: Alternative text

Modify calculated parameters
****************************

As input parameters, calculated parameters cna also be overridden. For exmaple here, we override the `driving mass`
of large diesel vehicles in, for all years::

    cm.array.loc['Large','ICEV-d','driving mass',:] = [[2000],[2200]]



Characterization of inventories (static)
****************************************

``carculator`` makes the characterization of inventories easy. You can characterize the inventories directly from
``carculator`` against midpoint, endpoint and single score impact assessment methods.

For example, to obtain characterized results against the midpoint impact assessment method ReCiPe for all cars::

    ic = InventoryCalculation(cm.array)
    results = ic.calculate_impacts()



Hence, to plot the carbon footprint for all medium cars in 2017::

    results.sel(size='Medium', year=2017, impact_category='climate change', value=0).to_dataframe('impact').unstack(level=1)['impact'].plot(kind='bar',
                stacked=True)
    plt.ylabel('kg CO2-eq./vkm')
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_carbon_footprint.png
    :width: 400
    :alt: Alternative text

Note that, for now, only the ReCiPe method is available for midpoint characterization. Also, once the instance of the :class:`CarModel`
class has been created, there is no need to re-create it in order to calculate additional environmental impacts (unless you wish to
change values of certain input or calculated parameters, the driving cycle or go from static to stochastic mode).

Characterization of inventories (stochastic)
********************************************

In the same manner, you can obtain distributions of results, instead of one-point values if you have run the model in
stochastic mode (with 500 iterations and the driving cycle WLTC)::

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

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_stochastic_BEV_PHEV.png
    :width: 400
    :alt: Alternative text

Many other examples are described in a Jupyter Notebook in the ``examples`` folder.