---
title: 'Carculator: prospective environmental and economic life cycle assessment of vehicles'
tags:
  - Python
  - environment
  - transport
  - mobility
  - carbon footprint
authors:
  - name: Romain Sacchi
    orcid: 0000-0003-1440-0905
    affiliation: 1
  - name: Christopher Mutel
    orcid: 0000-0002-7898-9862
    affiliation: 1
  - name: Christian Bauer
    affiliation: 1
  - name: Brian Cox
    orcid: 0000-0002-4531-6709
    affiliation: 1
affiliations:
 - name: Technology Assessment, Paul Scherrer Institut, Villigen, Switzerland
   index: 1
date: 18 December 2019
bibliography: paper.bib
---

# Summary

Life Cycle Assessment (LCA) is a tool that accounts for harmful emissions
occurring in all the relevant phases of the life cycle of a product or service [@ISO2006].
With the rising urgency for climate change mitigation, LCA has therefore emerged as a tool
to support decisions for the design of policies. It is often used to compare technologies on the basis of a 
common functional unit against some environmental indicators (e.g., emissions of greenhouse gases) [@Sala2016]. 

In the field of mobility, LCA studies have largely focused on comparing fossil fuel-powered vehicles to
battery electric or hydrogen-powered vehicles -- see for example [@Bauer2015].
While LCA is fit for such purpose, results typically show that assumptions made upstream
from the use phase can be important, as emphasized, among others, in [@Bauer2015;@Helmers2017].
Because the importance of such assumptions is not always transparent to the reader and because not all of the
information required to reproduce the results are readily available,
it is possible to find similar studies with contradicting conclusions. It was the case with a study by [@Buchal2019],
where several inconsistencies and unfounded assumptions in the fuel pathway of electric cars led
to conclusions in contrast with the rest of the scientific literature.

As a response to this situation, ``carculator`` was developed to perform LCA of different
vehicle technologies in a transparent, open-source, reproducible and efficient manner.
 
``carculator`` is a fully parameterized Python model that performs prospective
LCA of passenger vehicles. The foreground modeling is based on a physical model that sizes vehicles of different types
and dimensions and calculates the energy to move them over a given distance based on a driving cycle.
The background modeling relies on the extensive life cycle inventory database ecoinvent v.3 [@Wernet2016].

Initially developed by @Cox2018, the code has been refactored into a library to conduct
fast calculations and offer a convenient way to control parameters both in
the foreground (i.e., vehicles) and background (i.e., fuel pathways and energy storage) aspects
of the model. ``carculator`` also handles uncertainty in parameters and can perform error propagation analyses
relatively quickly.

Therefore, performing complex analyses becomes relatively easy.
In this example, an error propagation analysis between a diesel-powered vehicle, bio-ethanol-powered vehicle and a
battery electric vehicle in 2040 in terms of Ozone Depletion Potential (ODP) is performed over 1,000 iterations.
In addition, an electricity mix for the battery charge based solely on hydro-power is specified, along with a lithium
iron phosphate (LFP) battery type, manufactured in Norway.

```python
    from carculator import *

    background_configuration = {
    # in this case, 100% hydropower
    'custom electricity mix' : [[1,0,0,0,0,0,0,0,0,0]],            
    'petrol technology': 'bioethanol - wheat straw',
    'battery technology': 'LFP',
    'battery origin': 'NO'
    }
    scope = {
        'powertrain':['ICEV-d', 'BEV', 'ICEV-p'],
        'size':['Large'],
        'year':[2040]
    }
    cip = CarInputParameters()
    # 1000 iterations
    cip.stochastic(1000)
    _, array = fill_xarray_from_input_parameters(cip)
    cm = CarModel(array, cycle='WLTC')
    cm.set_all()
    ic = InventoryCalculation(
                cm.array,
                scope=scope,
                background_configuration=background_configuration
                )
    results = ic.calculate_impacts()
```

![ODP results per vehicle-kilometer, 1,000 iterations](https://github.com/romainsacchi/coarse/raw/master/docs/MC_example_article.png)

The analysis is done in 72 seconds on a modern laptop.

Additionally, users can export and share the inventories of the car models as well as the uncertainty data,
to reuse them in other LCA tools such as Brightway2 [@Mutel2017].

``carculator`` benefits from an [online documentation](https://coarse-lci.readthedocs.io/en/latest/) as well as a [graphical user interface](http://carculator.psi.ch).


# Acknowledgements

The authors would like to acknowledge the financial contribution of InnoSuisse via the project
Swiss Competence Center for Energy Research (SCCER) Efficient Technologies and Systems for Mobility.
This work was also financially supported by ACT ELEGANCY, Project No 271498, which has received funding
from DETEC (CH), BMWi (DE), RVO (NL), Gassnova (NO), BEIS (UK), Gassco, Equinor and Total, and is
cofunded by the European Commission under the Horizon 2020 programme, ACT Grant Agreement No 691712.

# References
