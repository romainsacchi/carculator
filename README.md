# Carculator

A fully parameterized Python model to calculate the life cycle material and energy requirements of passenger cars.

Integrates well with [Brightway](https://brightwaylca.org/) and [presamples](https://github.com/PascalLesage/brightway2-presamples).

Extracted from [Uncertain environmental footprint of current and future battery electric vehicles by Cox, et al (2018)](https://pubs.acs.org/doi/abs/10.1021/acs.est.8b00261).

[![Build Status](https://travis-ci.org/romainsacchi/carculator.svg?branch=master)](https://travis-ci.org/romainsacchi/carculator) [![Build status](https://ci.appveyor.com/api/projects/status/github/romainsacchi/coarse?svg=true)](https://ci.appveyor.com/project/romainsacchi/coarse) [![Coverage Status](https://coveralls.io/repos/github/romainsacchi/coarse/badge.svg)](https://coveralls.io/github/romainsacchi/coarse) [![Documentation](https://readthedocs.org/projects/coarse_lci/badge/?version=latest)](https://coarse-lci.readthedocs.io/en/latest/)

## Documentation

See [Documentation](https://coarse-lci.readthedocs.io/en/latest/index.html).

## Usage

As an example, comparing the carbon footprint electric vehicles to rechargeable hybrid vehicles for different sizes today and in the future
over 500 Monte Carlo iterations is as easy as:

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

For examples, see [examples](https://github.com/romainsacchi/carculator/blob/master/examples/Examples.ipynb).

## Web graphical user interface

See [carculator_online](http://carculator.psi.ch).

## Installation from conda

    conda install -c conda-forge -c pascallesage -c cmutel -c romainsacchi/label/nightly carculator-dev

## Installation from GitHub

    pip install git+https://github.com/romainsacchi/carculator.git

will install the package in editable mode and the required dependencies.

## Contributing

Details on how to contribute: see Issues.
