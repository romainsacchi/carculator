Introduction
============

**Carculator** is a parametrized model that allows to generate and characterize life cycle inventories for 112 different
passenger car configurations, according to selected:

* powertrain technologies (8): petrol engine, diesel engine, electric motor, hybrid, etc.,
* year of operation (2): 2017, 2040
* and sizes (7): Mini, Large, etc.

It is initially based on the model developed in `Uncertain environmental footprint of current and future battery electric
vehicles by Cox, et al (2018) <https://pubs.acs.org/doi/10.1021/acs.est.8b00261>`_.

More specifically, **Carculator** generates Brightway2- and Simapro-compatible inventories, but also directly provides characterized
results against several mid- and endpoint indicators (e.g., ReCiPe, Ecological Scarcity) and life cycle cost indicators.

Why?
----
Many life cycle assessment models of passenger cars exist, often leading to varying conclusions being published.
Unfortunately, because the underlying calculations are kept undocumented, it is not possible to explain the disparity
in the results given by these models.

Because **Carculator** is kept as open as possible, the methods and assumptions behind the generation of results are easily identifiable.
Also, there is an effort to keep the different modules (classes) separated, so that improving certain areas of the model is relatively
easy and does not require changing extensive parts of the code. In that regard, contributions are welcome.

Finally, beside being more flexible and transparent, **Carculator** provides interesting features, such as:

* a stochastic mode, that allows fast Monte Carlo analyses (1,000 iterations in 22 seconds, on a computer equipped with an i5 CPU)
* possibility to override any or all of the 200+ default car parameters (e.g., number of passengers, drag coefficient)
* hot pollutants emissions function of driving cycle, using HBEFA 3.3 data, further divided between rural, suburban and urban areas
* noise emissions, based on `CNOSSOS-EU <https://ec.europa.eu/jrc/en/publication/reference-reports/common-noise-assessment-methods-europe-cnossos-eu>`_ models for noise emissions and `Noise footprint from personal land‐based mobility by Cucurachi, et al (2019) <https://onlinelibrary.wiley.com/doi/full/10.1111/jiec.12837>`_ for inventory modelling and mid- and endpoint characterization of noise emissions, function of driving cycle and further divided between rural, suburban and urban areas
* exports inventories as csv files to be used with Brightway2 or Simapro (in progress), including uncertainty information. This requires
the user to have `ecoinvent 3.5 cutoff` installed on the LCA software the car inventories are exported to.
* exports inventories directly into Brightway2, as a database object to be registered. Additionally, when run in stochastic mode,
it is possible to export arrays of pre-sampled values using the `presamples <https://pypi.org/project/presamples/>`_ library
to be used together with the Monte Carlo function of Brightway2.

Objective
---------

The objective is to produce life cycle inventories for passenger cars in a transparent and comprehensive way,
to be further used in prospective life cycle assessment of transportation technologies.

How to install this package?
----------------------------

**Carculator** is a Python package, and is primarily to be used from within a Python 3.x environment.
Because **Carculator** is still at an early development stage, it is a good idea to install it in a separate environment,
such as a conda environment::

    conda create -n <name of the environment> python=3.7

Once your environment created, you should activate it::

    conda activate <name of the environment>

And install the development version of the **Carculator** package in your new environment via Conda::

    conda install -c conda-forge -c pascallesage -c cmutel -c romainsacchi/label/nightly carculator-dev

This will install the package and the required dependencies.


How to use it?
--------------

From a Python environment::

   from carculator import *
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

generates a CarModel object and calculate the energy consumption, components mass, exhaust and non-exhaust emissions for all vehicle profiles.



Accessing attributes of the car model
*************************************
Hence, the tank-to-wheel energy requirement per km driven per powertrain technology for a SUV in 2017 can be obtained::

    import numpy as np
    TtW_energy = cm.array.sel(size='SUV', year=2017, parameter='TtW energy') * 1/3600 * 100

    plt.bar(TtW_energy.powertrain,np.squeeze(TtW_energy))
    plt.ylabel('kWh/100 km')
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/fig_kwh_100km.png
    :width: 400
    :alt: Alternative text

Any other attributes of the CarModel class can be obtained in a similar way. For example, all calculated parameters that start with
`_lci` are inputs or outputs to the inventory.
For example, the following lists all direct exhaust emissions included in the inventory of an petrol Van in 2017:

List of all the given and calculated parameters of the car model::

    list_param = cm.array.coords['parameter'].values.tolist()

Return the parameters concerned with direct exhaust emissions (we remove noise emissions)::

    direct_emissions = [x for x in list_param if x.startswith('_lci_direct_') and 'noise' not in x]

Finally, return their values and display the first 10 in a table::

    cm.array.sel(parameter=direct_emissions, year=2017, size='Van', powertrain='BEV').to_dataframe(name='direct emissions')

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_direct_emissions.png
    :width: 400
    :alt: Alternative text


Or we could be interested in visualizing teh distriution of uncharacterized noise emissions, in joules::

    noise_emissions = [x for x in list_param if 'noise' in x]
    data = cm.array.sel(parameter=noise_emissions, year=2017, size='Van', powertrain='ICEV-p', value=0)\
        .to_dataframe(name='noise emissions')['noise emissions']
    data[data>0].plot(kind='bar')
    plt.ylabel('joules per km')

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/example_noise_emissions.png
    :width: 400
    :alt: Alternative text






