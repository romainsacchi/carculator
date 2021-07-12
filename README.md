# ``carculator``

<p align="center">
  <img style="height:130px;" src="https://github.com/romainsacchi/coarse/raw/master/docs/mediumsmall.png">
</p>

<p align="center">
  <a href="https://badge.fury.io/py/carculator" target="_blank"><img src="https://badge.fury.io/py/carculator.svg"></a>
  <a href="https://github.com/romainsacchi/carculator" target="_blank"><img src="https://github.com/romainsacchi/carculator/actions/workflows/main.yml/badge.svg?branch=master"></a>
  <a href="https://ci.appveyor.com/project/romainsacchi/carculator" target="_blank"><img src="https://ci.appveyor.com/api/projects/status/github/romainsacchi/carculator?svg=true"></a>
  <a href="https://coveralls.io/github/romainsacchi/carculator" target="_blank"><img src="https://coveralls.io/repos/github/romainsacchi/carculator/badge.svg"></a>
  <a href="https://carculator.readthedocs.io/en/latest/" target="_blank"><img src="https://readthedocs.org/projects/carculator/badge/?version=latest"></a>
  <a href="https://doi.org/10.5281/zenodo.3778259"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3778259.svg" alt="DOI"></a>
</p>

Prospective environmental and economic life cycle assessment of vehicles made blazing fast.

A fully parameterized Python model developed by the [Technology Assessment group](https://www.psi.ch/en/ta) of the
[Paul Scherrer Institut](https://www.psi.ch/en) to perform life cycle assessments (LCA) of vehicles.

See [the documentation](https://carculator.readthedocs.io/en/latest/index.html) for more detail, validation, etc.

## Table of Contents

- [Background](#background)
  - [What is Life Cycle Assessment](#what-is-life-cycle-assessment)
  - [Why carculator](#why-carculator)
- [Install](#install)
- [Usage](#usage)
  - [As a Python library](#as-a-python-library)
  - [As a web app](#as-a-web-app)
- [Support](#support)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Background

### What is Life Cycle Assessment?

Life Cycle Assessment (LCA) is a systematic way of accounting for environmental impacts along the relevant phases of the life of a product or service.
Typically, the LCA of a passenger vehicle includes the raw material extraction, the manufacture of the vehicle, its distribution, use and maintenance, as well as its disposal.
The compiled inventories of material and energy required along the life cycle of the vehicle is characterized against some impact categories (e.g., climate change).

In the research field of mobility, LCA is widely used to investigate the superiority of a technology over another one.

### Why ``carculator``?

``carculator`` allows to:
* produce [life cycle assessment (LCA)](https://en.wikipedia.org/wiki/Life-cycle_assessment) results that include conventional midpoint impact assessment indicators as well cost indicators
*  ``carculator`` uses time- and energy scenario-differentiated background inventories for the future, based on outputs of Integrated Asessment Model [REMIND](https://www.pik-potsdam.de/research/transformation-pathways/models/remind/remind). 
* calculate hot pollutant and noise emissions based on a specified driving cycle
* produce error propagation analyzes (i.e., Monte Carlo) while preserving relations between inputs and outputs
* control all the parameters sensitive to the foreground model (i.e., the vehicles) but also to the background model
(i.e., supply of fuel, battery chemistry, etc.)
* and easily export the vehicle models as inventories to be further imported in the [Brightway2](https://brightwaylca.org/) LCA framework
  or the [SimaPro](https://www.simapro.com/) LCA software.

``carculator`` integrates well with the [Brightway](https://brightwaylca.org/) LCA framework.

``carculator`` was built based on work described in [Uncertain environmental footprint of current and future battery electric vehicles by Cox, et al (2018)](https://pubs.acs.org/doi/abs/10.1021/acs.est.8b00261).

## Install

``carculator`` is at an early stage of development and is subject to continuous change and improvement.
Three ways of installing ``carculator`` are suggested.

We recommend the installation on **Python 3.7 or above**.

### Installation of the latest version, using conda

    conda install -c romainsacchi carculator

### Installation of a stable release (1.3.1) from Pypi

    pip install carculator

## Usage

### As a Python library

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

## As a Web app

``carculator`` has a [graphical user interface](https://carculator.psi.ch) for fast comparisons of vehicles.

## Support

Do not hesitate to contact the development team at [carculator@psi.ch](mailto:carculator@psi.ch).

## Maintainers

* [Romain Sacchi](https://github.com/romainsacchi)
* [Chris Mutel](https://github.com/cmutel/)

## Contributing

See [contributing](https://github.com/romainsacchi/carculator/blob/master/CONTRIBUTING.md).

## License

[BSD-3-Clause](https://github.com/romainsacchi/carculator/blob/master/LICENSE). Copyright 2020 Paul Scherrer Institut.
