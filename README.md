# Carculator

<p align="center">
  <img style="height:130px;" src="https://github.com/romainsacchi/coarse/raw/master/docs/mediumsmall.png">
</p>

<p align="center">
  <a href="https://travis-ci.org/romainsacchi/carculator" target="_blank"><img src="https://travis-ci.org/romainsacchi/carculator.svg?branch=master"></a>
  <a href="https://ci.appveyor.com/project/romainsacchi/carculator" target="_blank"><img src="https://ci.appveyor.com/api/projects/status/github/romainsacchi/carculator?svg=true"></a>
  <a href="https://coveralls.io/github/romainsacchi/coarse" target="_blank"><img src="https://coveralls.io/repos/github/romainsacchi/coarse/badge.svg"></a>
  <a href="https://coarse-lci.readthedocs.io/en/latest/" target="_blank"><img src="https://readthedocs.org/projects/coarse_lci/badge/?version=latest"></a>
</p>

A fully parameterized Python model developed by the [Technology Assessment group](https://www.psi.ch/en/ta) of the
[Paul Scherrer Institut](https://www.psi.ch/en) to calculate the life cycle material and energy requirements of passenger cars.

More specifically, *Carculator* allows to:
* produce life cycle assessment results that include conventional midpoint impact assessment indicators as well cost indicators
* calculate hot pollutant and noise emissions based on a specified driving cycle
* produce error propagation analyzes (i.e., Monte Carlo) while preserving relations between inputs and outputs
* control all the parameters sensitive to the foreground model (i.e., the vehicles) but also to the background model
(i.e., supply of fuel, battery chemistry, etc.)
* and easily export the vehicle models as inventories to be further imported in the `brightway2` LCA framework

*Carculator* integrates well with [Brightway](https://brightwaylca.org/) and [presamples](https://github.com/PascalLesage/brightway2-presamples).

Extended from the initial work described in [Uncertain environmental footprint of current and future battery electric vehicles by Cox, et al (2018)](https://pubs.acs.org/doi/abs/10.1021/acs.est.8b00261).


## Documentation

See [Documentation](https://coarse-lci.readthedocs.io/en/latest/index.html).

## Usage

Calculate the ``Tank to wheel`` energy requirement (in kWh/100 km) of current SUVs for the driving cycle WLTC 3.4
over 800 Monte Carlo iterations:

    from carculator import *
    cip = CarInputParameters()
    cip.stochastic(800)
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array, cycle='WLTC 3.4')
    cm.set_all()
    TtW_energy = cm.array.sel(size='SUV', year=2017, parameter='TtW energy', value=0) * 1/3600 * 100

    plt.bar(TtW_energy.powertrain, TtW_energy)
    plt.ylabel('kWh/100 km')
    plt.show()
    
![MC results](https://github.com/romainsacchi/coarse/raw/master/docs/stochastic_example_ttw.png)

Compare the carbon footprint of electric vehicles with that of rechargeable hybrid vehicles for different size categories today and in the future
over 500 Monte Carlo iterations:

    from carculator import *
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
    
![MC results](https://github.com/romainsacchi/coarse/raw/master/docs/example_stochastic_BEV_PHEV.png)

For more examples, see [examples](https://github.com/romainsacchi/carculator/blob/master/examples/Examples.ipynb).

## Web graphical user interface

Carculator has a graphical user interface for fast comparisons of vehicles.
See [carculator_online](http://carculator.psi.ch).

## Installation

*Carculator* is at an early stage of development and is subject to continuous change and improvement.
Three ways of installing *Carculator* are suggested.

### Installation of a stable version from Pypi

    pip install carculator

### Installation of a development version from from conda

    conda install -c conda-forge -c pascallesage -c cmutel -c romainsacchi/label/nightly carculator-dev
    

### Installation of a development version from from GitHub

    pip install git+https://github.com/romainsacchi/carculator.git


## Contributing

Details on how to contribute: see [Issues](https://github.com/romainsacchi/carculator/issues).

## Support

Do not hesitate to contact the development team at [carculator@psi.ch](carculator@psi.ch).