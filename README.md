# ``carculator``

<p align="center">
  <img style="height:130px;" src="https://github.com/romainsacchi/coarse/raw/master/docs/mediumsmall.png">
</p>

<p align="center">
  <a href="https://travis-ci.org/romainsacchi/carculator" target="_blank"><img src="https://travis-ci.org/romainsacchi/carculator.svg?branch=master"></a>
  <a href="https://ci.appveyor.com/project/romainsacchi/carculator" target="_blank"><img src="https://ci.appveyor.com/api/projects/status/github/romainsacchi/carculator?svg=true"></a>
  <a href="https://coveralls.io/github/romainsacchi/coarse" target="_blank"><img src="https://coveralls.io/repos/github/romainsacchi/coarse/badge.svg"></a>
  <a href="https://coarse-lci.readthedocs.io/en/latest/" target="_blank"><img src="https://readthedocs.org/projects/coarse_lci/badge/?version=latest"></a>
</p>

## What is ``carculator``?

A fully parameterized Python model developed by the [Technology Assessment group](https://www.psi.ch/en/ta) of the
[Paul Scherrer Institut](https://www.psi.ch/en) to perform life cycle assessments (LCA) of passenger cars.

## What is Life Cycle Assessment?

Life Cycle Assessment (LCA) is a systematic way of accounting for environmental impacts along the relevant phases of the life of a product or service.
Typically, the LCA of a passenger vehicle includes the raw material extraction, the manufacture of the vehicle, its distribution, use and maintenance, as well as its disposal.
The compiled inventories of material and energy required along the life cycle of the vehicle is characterized against some impact categories (e.g., climate change).
In the research field of mobility, LCA is widely used to investigate the superiority of a technology over another one.

## Why is ``carculator`` needed?

``carculator`` allows to:
* produce [life cycle assessment (LCA)](https://en.wikipedia.org/wiki/Life-cycle_assessment) results that include conventional midpoint impact assessment indicators as well cost indicators
* calculate hot pollutant and noise emissions based on a specified driving cycle
* produce error propagation analyzes (i.e., Monte Carlo) while preserving relations between inputs and outputs
* control all the parameters sensitive to the foreground model (i.e., the vehicles) but also to the background model
(i.e., supply of fuel, battery chemistry, etc.)
* and easily export the vehicle models as inventories to be further imported in the [Brightway2](https://brightwaylca.org/) LCA framework

``carculator`` integrates well with [Brightway](https://brightwaylca.org/) and [presamples](https://github.com/PascalLesage/brightway2-presamples).

Extended from the initial work described in [Uncertain environmental footprint of current and future battery electric vehicles by Cox, et al (2018)](https://pubs.acs.org/doi/abs/10.1021/acs.est.8b00261).


## Documentation

See [Documentation](https://coarse-lci.readthedocs.io/en/latest/index.html).

## Installation

``carculator`` is at an early stage of development and is subject to continuous change and improvement.
Three ways of installing ``carculator`` are suggested.

### Installation of a stable release (0.0.4) from Pypi

    pip install carculator

### Installation of a development version from conda

    conda install -c conda-forge -c cmutel -c romainsacchi/label/nightly carculator
    

### Installation of a development version from GitHub

    pip install git+https://github.com/romainsacchi/carculator.git


## Usage examples

Calculate the fuel efficiency (or ``Tank to wheel`` energy requirement) in km/L of petrol-equivalent of current SUVs for the driving cycle WLTC 3.4
over 800 Monte Carlo iterations:
```python
    from carculator import *
    import matplotlib.pyplot as plt
    cip = CarInputParameters()
    cip.stochastic(800)
    dcts, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array, cycle='WLTC 3.4')
    cm.set_all()
    TtW_energy = 1 / (cm.array.sel(size='SUV', year=2017, parameter='TtW energy') / 42000) # assuming 42 MJ/L petrol
    
    l_powertrains = TtW_energy.powertrain
    [plt.hist(e, bins=50, alpha=.8, label=e.powertrain.values) for e in TtW_energy]
    plt.xlabel('km/L petrol-equivalent')
    plt.ylabel('number of iterations')
    plt.legend()
```
   
![MC results](https://github.com/romainsacchi/coarse/raw/master/docs/stochastic_example_ttw.png)

Compare the carbon footprint of electric vehicles with that of rechargeable hybrid vehicles for different size categories today and in the future
over 500 Monte Carlo iterations:
```python
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
```
    
![MC results](https://github.com/romainsacchi/coarse/raw/master/docs/example_stochastic_BEV_PHEV.png)

For more examples, see [examples](https://github.com/romainsacchi/carculator/blob/master/examples/Examples.ipynb).

## Web graphical user interface

``carculator`` has a graphical user interface for fast comparisons of vehicles.
See [carculator_online](http://carculator.psi.ch).


## Contributing

Details on how to contribute: see [Issues](https://github.com/romainsacchi/carculator/issues).

## Support

Do not hesitate to contact the development team at [carculator@psi.ch](mailto:carculator@psi.ch).