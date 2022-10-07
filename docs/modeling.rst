Sections
********

`1 Overview of carculator modules <#overview-of-carculator-modules>`__

`2 Vehicle modelling <#vehicle-modelling>`__

`2.1 Size classes <#size-classes>`__

`2.2 Manufacture year and emission standard <#manufacture-year-and-emission-standard>`__

`2.3 Size and mass-related parameters and modeling <#size-and-mass-related-parameters-and-modeling>`__

`2.4 Electric energy storage <#electric-energy-storage>`__

`2.5 Fuel cell stack <#fuel-cell-stack>`__

`2.6 Light weighting <#light-weighting>`__

`2.7 Sizing of onboard energy storage <#sizing-of-onboard-energy-storage>`__

`2.7.1 Sizing of battery <#sizing-of-battery>`__

`2.8 Electric utility factor <#electric-utility-factor>`__

`3 Inventory modelling <#inventory-modelling>`__

`3.1 Road demand <#road-demand>`__

`3.2 Fuel properties <#fuel-properties>`__

`3.3 Exhaust emissions <#exhaust-emissions>`__

`3.3.1 NMHC speciation <#nmhc-speciation>`__

`3.4 Non-exhaust emissions <#non-exhaust-emissions>`__

`3.4.1 Engine wear emissions <#engine-wear-emissions>`__

`3.4.2 Abrasion emissions <#abrasion-emissions>`__

`3.4.3 Refrigerant emissions <#refrigerant-emissions>`__

`3.5 Noise emissions <#noise-emissions>`__

`3.6 Electricity mix calculation <#electricity-mix-calculation>`__

`3.7 Inventories for fuel pathways <#inventories-for-fuel-pathways>`__

`3.8 Inventories for energy storage components <#inventories-for-energy-storage-components>`__

`4 Life cycle impact assessment <#life-cycle-impact-assessment>`__

`References <#references>`__

This document intends to describe the *carculator* model, assumptions
and inventories as exhaustively as possible.
*carculator* is an open-source Python library. Its code is publicly
available via its `Github
repository <https://github.com/romainsacchi/carculator>`__. There is
also `an examples
notebook <https://github.com/romainsacchi/carculator/blob/master/examples/Examples.ipynb>`__,
to guide new users into performing life cycle analyses. Finally, there
is also an online graphical user interface available at
https://carculator.psi.ch.

Overview of *carculator* modules
********************************

The main module *model.py* builds
the vehicles and delegates the calculation of motive and auxiliary
energy, noise, abrasion and exhaust emissions to satellite modules. Once
the vehicles are fully characterized, the set of calculated parameters
are passed to *inventory.py* which derives life cycle inventories and
calculates life cycle impact indicators. Eventually, these inventories
can be passed to *export.py* to be exported to various LCA software.


Vehicle modelling
*****************

The modelling of vehicles along powertrain types, time and size classes
is described in this section. It is also referred to as *foreground*
modelling.

Size classes
------------

Originally, *carculator* defines nine size classes, namely: Micro, Mini,
Small, Lower medium, Medium, Large, Medium SUV, Large SUV and Van,
according to the following criteria, adapted from the work of [1]_ and
shown in Table 1.

Table 1 Criteria for size classes

+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| EU segment  | EU segment definition          | carculator    | Minimum footprint [m2]                                | Maximum footprint [m2]  | Minimum curb mass [kg]  | Maximum curb mass [kg]  | Examples                                                               |
+=============+================================+===============+=======================================================+=========================+=========================+=========================+========================================================================+
| L7e         | Micro cars/Heavy quadricycles  | Micro         |                                                       |                         | 400                     | 600                     | Renault Twizzy, Microlino                                              |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| A           | Mini cars                      | Mini          |                                                       | 3.4                     |                         | 1050                    | Renault Twingo, Smart ForTwo, Toyota Aygo                              |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| B           | Small cars                     | Small         | 3.4                                                   | 3.8                     | 900                     | 1'100                   | Renault Clio, VW Polo, Toyota Yaris                                    |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| C           | Medium cars                    | Lower medium  | 3.8                                                   | 4.3                     | 1'250                   | 1'500                   | VW Golf, Ford Focus, Mercedes Class A                                  |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| C           |                                | Medium        | 3.9                                                   | 4.4                     | 1'500                   | 1'750                   | VW Passat, Audi A4, Mercedes Class C                                   |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| D           | Large/Executive                | Large         | 4.4                                                   |                         | 1'450                   | 2'000                   | Tesla Model 3, BMW 5 Series, Mercedes E series                         |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| J           | Sport Utility                  | Medium SUV    | 4.5                                                   |                         | 1'300                   | 2'000                   | Toyota RAV4, Peugeot 2008, Dacia Duster                                |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| J           | Sport Utility                  | Large SUV     | 6                                                     |                         | 2'000                   | 2'500                   | Audi Q7, BMW X7, Mercedes-Benz GLS, Toyota Landcruiser, Jaguar f-Pace  |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+
| M           | Multi-Purpose Vehicles         | Van           | Defined by body type rather than mass and footprint.  |                         |                         |                         | VW Transporter, Mercedes Sprinter, Ford Transit                        |
+-------------+--------------------------------+---------------+-------------------------------------------------------+-------------------------+-------------------------+-------------------------+------------------------------------------------------------------------+

Example of Micro car (Microlino), Mini car (Smart) and Small/Compact car (VW Polo)

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image1.png
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image2.jpeg
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image3.png
    :width: 30%

Example of Lower medium car (VW Golf), Medium car (Peugeot 408) and Large car (Tesla Model 3)


.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image4.png
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image5.jpeg
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image6.png
    :width: 30%

Example of Medium SUV car (Peugeot 2008), Large SUV car (Audi Q7) and Van (Fiat Ducato)

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image7.png
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image8.png
    :width: 30%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image9.png
    :width: 30%

**Important remark**: Micro cars are not considered passenger cars in
the Swiss and European legislation, but heavy quadricycles. We do
however assimilate them to passenger cars. They are only modelled with a
battery electric powertrain.

**Important remark**: Sport Utility Vehicles (SUV) are considered more
as a body type than a size class. These vehicles have distinct
aerodynamic properties, but their curb mass can be as light as that of a
VW Polo or a Renault Clio (i.e., the Dacia Duster or Peugeot 2008 have a
curb mass of 1'150 kg, against 1'100-1'300 kg for a VW Polo) or as heavy
as a Mercedes Class E (i.e., the Audi Q7 has a minimum curb mass of
2'000 kg, against 1'900 kg for a Mercedes Class E). To assess the
impacts of very large SUV, the "Large SUV" category has been added, to
represent SUV models with a very high curb mass (2'000 kg and above) and
footprint.

Manufacture year and emission standard
--------------------------------------

Several emission standards are considered. For simplicity, it is assumed
that the vehicle manufacture year corresponds to the registration year,
as shown in Table 2.

Table 2 Correspondence between manufacture year and emission standards
used in *carculator*

+----------------+----------------+----------------+----------------+
|                | **Start of     | **End of       | **Manufacture  |
|                | registration** | registration   | year in**      |
|                |                | (incl.)**      | *carculator*   |
+================+================+================+================+
| **EURO-3**     | 2001           | 2005           | **2003**       |
+----------------+----------------+----------------+----------------+
| **EURO-4**     | 2006           | 2010           | **2008**       |
+----------------+----------------+----------------+----------------+
| **EURO-5**     | 2011           | 2014           | **2013**       |
+----------------+----------------+----------------+----------------+
| **EURO-6 a/b** | 2015           | 2017           | **2016**       |
+----------------+----------------+----------------+----------------+
| **EURO-6 c**   | 2018           | 2018           | **2018**       |
+----------------+----------------+----------------+----------------+
| **EURO-6 d     | 2019           | 2020           | **2019**       |
| (temp)**       |                |                |                |
+----------------+----------------+----------------+----------------+
| **EURO-6 d**   | 2021           | -              | **2021**       |
+----------------+----------------+----------------+----------------+

Size and mass-related parameters and modeling
---------------------------------------------

The vehicle glider and its components (powertrain, energy storage, etc.)
are sized according to engine power, which itself is conditioned by the
curb mass of the vehicle. The curb mass of the vehicle is the sum of
the vehicle components (excluding the driver and possible cargo) as
represented in Figure 2.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image10.png
   :width: 100%

Figure 2 Vehicle mass calculation workflow

This is an iterative process that stops when the curb mass of the
vehicle converges, as illustrated in Figure 3.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image11.png
   :width: 6.26667in
   :height: 3.99167in

Figure 3 Representation of the convergence of the sizing of the
passenger car model

Curb mass of the vehicle
________________________

.. math::

    m_{curb} = sum(m_{glider}, m_{charger}, m_{conv},
            m_{inv}, m_{distr}, m_{comb}, m_{elec},\\
            m_{pwt}, m_{fcstack}, m_{fcauxbop}, m_{fcessbop},
            m_{battcell}, m_{battbop}, m_{fueltank}, m_{fuel})

With:

-  :math:`m_{curb}` being the vehicle curb mass, in kg

-  :math:`m_{fuel}` being the fuel mass, in kg

-  :math:`m_{charger}` being the electric onboard charge mass, in kg

-  :math:`m_{conv}` being the current converter, in kg

-  :math:`m_{inv}` being the current AC/DC inverter, in kg

-  :math:`m_{distr}` being the power distribution unit, in kg

-  :math:`m_{comb}` being the combustion engine mass, in kg

-  :math:`m_{elec}` being the electric motor mass, in kg

-  :math:`m_{pwt}` being the powertrain mass, in kg

-  :math:`m_{fcstack}` being the fuel cell stack mass, in kg

-  :math:`m_{fcauxbop}` being the fuel cell auxiliary components mass, in kg

-  :math:`m_{battcell}` being the battery cell mass, in kg

-  :math:`m_{battbop}` being the battery auxiliary components mass, in kg

-  :math:`m_{fcessbop}` being the fuel cell essential components mass, in kg

-  :math:`m_{fueltank}` being the fuel tank mass, in kg


For each iteration, the tank-to-wheel energy consumption (i.e., the
motive energy minus any recuperated braking energy, together with the
needed auxiliary energy to power onboard electronics) of the vehicle is
calculated (i.e., to size the energy storage components, calculate the
fuel consumption, etc.), as described later in this section.

Cargo and driving mass of the vehicle
-------------------------------------

The cargo mass of the vehicle is the sum of the cargo mass and the
passenger mass.

.. math::

    m_{cargo} = m_{cargo} + m_{passenger}

where :math:`m_{cargo}` is the cargo mass, in kg,
and :math:`m_{passenger}` is the passenger mass.

The driving mass of the vehicle is the sum of the curb mass and the cargo mass.

.. math::

    m_{driving} = m_{curb} + m_{cargo}

where :math:`m_{curb}` is the curb mass, in kg, and :math:`m_{cargo}` is the
cargo mass, in kg, and :math:`m_{driving}` is the driving mass, in kg.

Light-weight rates
------------------

Because the LCI dataset used to represent the glider of the vehicle is
not representative of todays' use of light-weighting materials, such as
aluminium (i.e., the dataset "glider for passenger cars" only contains
0.5% of its mass in aluminium) and advanced high strength steel (AHSS),
an amount of such light-weighting materials is introduced to substitute
conventional steel and thereby reduce the mass of the glider.

As further explained in Section 2.6, the mass of the glider is reduced
by replacing steel with a mix of aluminium and AHSS. Hence, the amounts
of light weighting materials introduced depend on the rate of glider
light weighting in 2020 relative to 2000 (approximately 11% for
combustion engine vehicles). The amount of aluminium introduced is
further cross-checked with the amounts indicated in [2]_ and listed in
Table 3, and comes in addition to the aluminium already contained in the
LCI datasets for the engine and transmission.

**Important remark:** the light-weighting rate is for most vehicles
approximately 11% in 2020 relative to 2000. However, battery-equipped vehicles
are an exception to this: Medium, Large and Large SUV vehicles have
significantly higher light weighting rates to partially compensate for
the additional mass of their batteries. In order to match the battery
capacity and the curb mass of their respective size class, their light
weighting rate is increased to 14, 28 and 30%, respectively. This trend
is also confirmed by [2]_, showing that battery electric vehicles have
85% more aluminium than combustion engine vehicles, partly going into
the battery management system, and partly going into the chassis to
compensate for the extra mass represented by the battery.

These light-weighting rates have been fine adjusted to match the curb
mass of a given size class, while preserving the battery capacity. For
example, in the case of the Large SUV, its curb mass should
approximately be 2'200 kg, with an 80 kWh battery weighting 660 kg
(e.g., Jaguar i-Pace). This is possible with a 30% light weighting rate,
introducing approximately 460 kg of aluminium in the chassis (which
matches roughly with the value given for an Audi e-Tron in Table 3) and
1'008 kg of AHSS, in lieu of 2'034 kg of regular steel.

Table 3 Amount of aluminium in European passenger cars. Source: [2]_

+-------------------------------------------------------------------------------+--------+--------------+----------+------+----------+--------+------+------------------+
| Used in source                                                                | Basic  | Sub-Compact  | Compact  |      | Midsize  | Large  |      | Audi e-Tron      |
+===============================================================================+========+==============+==========+======+==========+========+======+==================+
| Used in carculator                                                            | Small  |              |          |      | Medium   | Large  |      | Large SUV (BEV)  |
+-------------------------------------------------------------------------------+--------+--------------+----------+------+----------+--------+------+------------------+
| Average aluminium content per vehicle [kg]                                    | 77     | 98           | 152      |      | 266      | 442    |      | 804              |
+-------------------------------------------------------------------------------+--------+--------------+----------+------+----------+--------+------+------------------+
| Share of aluminium mass in components other than engine and transmission [%]  | 66%    |              |          |      |          |        |      |                  |
+-------------------------------------------------------------------------------+--------+--------------+----------+------+----------+--------+------+------------------+
| Aluminium to be added to the glider [%]                                       | 65     |              |          | 175  |          | 292    | 530  |                  |
+-------------------------------------------------------------------------------+--------+--------------+----------+------+----------+--------+------+------------------+

Curb mass calibration
---------------------

The final curb mass obtained for each vehicle is calibrated against the
European Commission's database for CO\ :sub:`2` emission tests for
passenger cars (hereafter called EC-CO2-PC) using the NEDC/WLTP driving
cycles [3]_. Each vehicle registered in the European Union is tested and
several of the vehicle attributes are registered (e.g., dimension, curb
mass, driving mass, CO\ :sub:`2` emissions, etc.). This has represented
about 15+ million vehicles per year for the past five years.

The figure below shows such calibration for the years 2010, 2013, 2016,
2018, 2019 and 2020 -- to be representative of EURO-4, -5, 6 a/b, 6-c
and 6-d-temp vehicles. No measurements are available for 2003 (EURO-3)
or 2021 (EURO-6-d). After cleaning the data from the EC-CO2-PC database,
it represents 27 million points to calibrate the curb mass of the
vehicles with. Green vertical bars represent the span of 50% of the curb
mass distribution, and the red dots are the curb mass values modeled by
*carculator*.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image12.png

Figure 4 Calibration of the curb mass of the passenger car model against
the EC-CO2-PC database. Red dots: values modeled by carculator. Green
box-and-whiskers: values distribution from the EC-CO2-PC database (box:
50% of the distribution, whiskers: 90% of the distribution). Micro cars
are not represented in the EC-CO2-PC database. Sample size for each size
class is given above each chart. M = Mini, S = Small, L-M = Lower
medium, M = Medium, L = Large, L-SUV = Large SUV. Source for vehicle
tank-to-wheel energy consumption measurements: [4]_.

Table 4 shows the mass distribution for gasoline and battery electric
passenger cars resulting from the calibration. Mass information on other
vehicles is available in the vehicles' specifications spreadsheet.

Table 4 Mass distribution for gasoline and battery electric passenger
cars *in 2021*

+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Gasoline              |        |         |        | Battery electric  |        |        |         |        |            |
+=======================+========+=========+========+===================+========+========+=========+========+============+
| in kilograms          | Small  | Medium  | Large  | Large SUV         | Micro  | Small  | Medium  | Large  | Large SUV  |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Glider base mass      | 998    | 1'170   | 1'550  | 1'900             | 350    | 998    | 1'170   | 1'550  | 1'900      |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Light weighting       | -110   | -129    | -171   | -209              | -35    | -140   | -164    | -434   | -570       |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Glider mass           | 888    | 1'041   | 1'380  | 1'691             | 315    | 858    | 1'006   | 1'116  | 1'330      |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Powertrain mass       | 96     | 106     | 132    | 140               | 42     | 67     | 77      | 94     | 100        |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Engine or motor mass  | 111    | 125     | 157    | 168               | 29     | 61     | 73      | 96     | 102        |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Energy storage mass   | 72     | 85      | 104    | 104               | 120    | 276    | 360     | 580    | 660        |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Electronics mass      | 3      | 4       | 5      | 7                 | 23     | 23     | 23      | 23     | 23         |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Curb mass             | 1'170  | 1'361   | 1'777  | 2'110             | 529    | 1'285  | 1'540   | 1'910  | 2'215      |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Passenger mass        | 120    | 120     | 120    | 120               | 120    | 120    | 120     | 120    | 120        |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Cargo mass            | 20     | 20      | 20     | 20                | 20     | 20     | 20      | 20     | 20         |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+
| Driving mass          | 1'310  | 1'501   | 1'917  | 2'250             | 669    | 1'425  | 1'680   | 2'050  | 2'355      |
+-----------------------+--------+---------+--------+-------------------+--------+--------+---------+--------+------------+

Energy consumption
-------------------

The energy consumption model of *carculator* calculates the energy
required at the wheels by considering different types of resistance.
Some of these resistances are related to the vehicle size class. For
example, the frontal area of the vehicle influences the aerodynamic
drag. Also, the kinetic energy to overcome the vehicle's inertia is
influenced by the mass of the vehicle (which partially correlates to
with the size class or body type), but also by the acceleration required
by the driving cycle. Other resistances, such as the climbing effort,
are instead determined by the driving cycle (but the vehicle mass also
plays a role here). Once the energy required at the wheels is known, the
model goes on to calculate the energy required at the tank level by
considering additional losses along the drive train (i.e., axles,
gearbox, and engine losses). The different types of resistance
considered are depicted in Figure 5, and the module calculation workflow
is presented in Figure 6.

Powertrains that are partially or fully electrified have the possibility
to recuperate a part of the energy spent for propulsion during
deceleration or braking. The round-trip battery energy loss (which is
the sum of the charge and discharge battery loss, described in Figure 5)
is subtracted from the recuperated energy. For hybrid vehicles (i.e.,
HEV-p, HEV-d), this allows to downsize the combustion engine and improve
the overall tank-to-wheel efficiency, as explained in [5]_.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image13.png
   :width: 6.26667in
   :height: 2.49167in

Figure 5 Representation of the different types of resistance considered.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image15.png
   :width: 100%

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image14.png
   :width: 20%


Figure 6 Motive energy calculation workflow


Finally, for each second of the driving cycle, the auxiliary power load
is considered. It comprises an auxiliary base power load (i.e., to
operate onboard electronics), as well as the power load from heating and
cooling. While electric vehicles provide energy from the battery to
supply heating and cooling (i.e., thereby decreasing the available
energy available for traction), combustion vehicles recover enough waste
engine heat to supply adequate heating. The values considered for the
auxiliary base power load and for the power load for heating and cooling
are presented in Table 5. These values are averaged over the whole year,
based on maximum demand and share of operation.

Table 5 Auxiliary power demand. Source : [6]_

+----------------------------------+---------------------------+--------+---------------------+---------+--------+---------------------------+--------+--------+---------+--------+------------+
| Auxiliary power base demand [W]  | Heating power demand [W]  |        |                     |         |        | Cooling power demand [W]  |        |        |         |        |            |
+==================================+===========================+========+=====================+=========+========+===========================+========+========+=========+========+============+
|                                  |                           | Micro  | Small               | Medium  | Large  | Large SUV                 | Micro  | Small  | Medium  | Large  | Large SUV  |
+----------------------------------+---------------------------+--------+---------------------+---------+--------+---------------------------+--------+--------+---------+--------+------------+
| ICEV, HEV, PHEV                  | 94                        |        | Provided by engine  |         |        |                           |        | 250    | 320     | 350    | 350        |
+----------------------------------+---------------------------+--------+---------------------+---------+--------+---------------------------+--------+--------+---------+--------+------------+
| BEV, FCEV                        | 75                        | 200    | 250                 | 320     | 350    | 350                       | 0      | 250    | 320     | 350    | 350        |
+----------------------------------+---------------------------+--------+---------------------+---------+--------+---------------------------+--------+--------+---------+--------+------------+

.. math::

    P_{aux} = P_{base} + (P_{heating} \times D_{heating}) + (P_{cooling} \times D_{cooling})

where :math:`P_{base}` is the auxiliary base power load [W],
:math:`P_{heating}` is the power load for heating [W],
:math:`P_{cooling}` is the power load for cooling [W],
:math:`D_{heating}` is the demand for heating [0-1] (=0 for non-electric vehicles),
and :math:`D_{cooling}` is the demand for cooling [0-1].

To convert it into an energy consumption :math:`F_{aux}` [kj/km],
the auxiliary power load is multiplied by the time of the driving cycle
and divided by the distance driven:

.. math::

    F_{aux} = \frac{P_{aux} \times T}{D}

where :math:`T` is the driving cycle time [seconds] and D is the distance [m].

**Important remark:** Micro cars are not equipped with an air
conditioning system. Hence, their cooling energy requirement is set to
zero.

A driving cycle is used to calculate the tank-to-wheel energy required
by the vehicle to drive over one kilometer. For example, the WLTC
driving cycle comprises a mix of urban, sub-urban and highway driving.
It is assumed representative of average Swiss and European driving
profile - although this would likely differ in the case of intensive
mountain driving.

Figure 7 exemplifies such calculation for a medium battery electric
passenger car manufactured in 2020, using the WLTC driving cycle.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image16.png
   :width: 4.46667in
   :height: 2.69742in

Figure 7 Cumulated tank-to-wheel energy consumption, along the WLTC
driving cycle, for a mid-size battery electric vehicle from 2020

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image17.png
   :width: 3.84722in
   :height: 1.94653in

Figure 8 Driving cycle and related parameters

So, the tank-to-wheel energy consumption :math:`F_{ttw}` is the sum of the motive energy and the
energy required to power auxiliary equipment. It is calculated as:

.. math::

    F_{ttw} = F_{motive} + F_{aux}

where :math:`F_{motive}` is the motive energy, and :math:`F_{aux}` is the
auxiliary energy.

There are no fuel consumption measurements available for fuel cell
vehicles. Values found in the literature and from manufacturers data are
used to approximate the engine and transmission efficiency and to
calibrate the final energy consumption.

For diesel and gasoline hybrid vehicles, which are ICE vehicles equipped
with a small electric motor to allow for energy recuperation and
reducing the engine size, the drivetrain and engine efficiency are based
on [5]_ [7]_. The amount of energy recuperated is determined by the driving
cycle as well as the round-trip efficiency between the wheels and the
engine and cannot be superior to the power output of the engine. Further
on, the share of recuperated energy over the total negative motive
energy (i.e., the braking or deceleration energy) is used as a
discounting factor for brake wear particle emissions.

Engine and transmission efficiencies for the different powertrains are
fine-tuned until it aligns reasonably well with the fuel consumption
values from the EC-CO2-PC database, as shown in Figure 9.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image18.png
   :width: 100%

Figure 9 Validation of the tank-to-wheel energy consumption against the
monitoring program for passenger vehicle emissions database. Sample size
for each size class is given above each chart. M = Mini, S = Small, L-M
= Lower medium, M = Medium, L = Large, L-SUV = Large SUV. Horizontal
lines within the green boxes represent the median value. The green boxes
represent 50% of the distribution (25\ :sup:`th`-75\ :sup:`th`
percentiles). The whiskers represent 90% of the distribution
(5\ :sup:`th`-95\ :sup:`th` percentiles). Outliers are not shown. Source
for vehicle tank-to-wheel energy consumption measurements: [4]_

Electric energy storage
-----------------------

Battery electric vehicles can use different battery chemistries (Li-ion
NMC, Li-ion LFP, Li-ion NCA and Li-LTO) depending on the manufacturer's
preference or the location of the battery supplier. Unless specified
otherwise, all battery types are produced in China, as several sources,
among which BloombergNEF [8]_, seem to indicate that more than 75% of the
world's cell capacity is manufactured there.

Accordingly, the electricity mix used for battery cells manufacture and
drying, as well as the provision of heat are assumed to be
representative of the country (i.e., the corresponding providers are
selected from the LCI background database).

The battery-related parameters considered in *carculator* for 2020 are
shown in Table 6. For LFP batteries, "blade battery" or "cell-to-pack"
battery configurations are considered, as introduced by CATL [9]_ and BYD
[10]_, two major LFP battery suppliers in Asia. This greatly increases
the cell-to-pack ratio and the gravimetric energy density at the pack
level.

Overall, the gravimetric energy density values at the cell and system
levels presented in Table 6 are considered conservative: some
manufacturers perform significantly better than the average, and these
values tend to change rapidly over time, as it is being the focus of
much R&D. Hence, by 2050, the gravimetric energy density of NMC and NCA
cells are expected to reach 0.5 kWh/kg, while that of LFP cells plateaus
at 0.15 kWh/kg (but benefits from a high cell-to-pack ratio)..

Table 6 Specifications for the different battery types

+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
|                                                                             | Lithium Nickel Manganese Cobalt Oxide (LiNiMnCoO2) — NMC  | Lithium Iron Phosphate(LiFePO4) — LFP  | Lithium Nickel Cobalt Aluminum Oxide (LiNiCoAlO2) — NCA  | Source                              |
+=============================================================================+===========================================================+========================================+==========================================================+=====================================+
| Cell energy density [kWh/kg]                                                | 0.2 (0.5 in 2050)                                         | 0.15                                   | 0.23 (0.5 in 2050)                                       | [11]_                                |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Cell-to-pack ratio                                                          | 0.6 (0.65 in 2050)                                        | 0.8 (0.9 in 2050)                      | 0.5 (0.55 in 2050)                                       | [12]_                                |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Pack-level gravimetric energy density [kWh/kg]                              | 0.12                                                      | 0.12                                   | 0.14                                                     | Calculated from the two rows above  |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Share of cell mass in battery system [%]                                    | 70 to 80% (depending on chemistry, see third row above)   |                                        |                                                          | [5,12]_                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Maximum state of charge [%]                                                 | 100%                                                      | 100%                                   | 100%                                                     | [11,13]_                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Minimum state of charge [%]                                                 | 20%                                                       | 20%                                    | 20%                                                      |                                     |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Cycle life to reach 20% initial capacity loss  (80%-20% SoC charge cycle)   | 2'000                                                     | 7'000+                                 | 1'000                                                    | [14]_                                |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Corrected cycle life                                                        | 3'000                                                     | 7'000                                  | 1'500                                                    | Assumption                          |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Charge efficiency                                                           | 85% in 2020, 86% in 2050                                  |                                        |                                                          | [5,15]_ for passenger cars.          |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+
| Discharge efficiency                                                        | 88% in 2020, 89% in 2050                                  |                                        |                                                          | [5,16]_                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------+----------------------------------------+----------------------------------------------------------+-------------------------------------+

Note that the NMC battery cell used by default corresponds to a so-called NMC 6-2-2 chemistry: it exhibits three times the mass amount of *Ni*
compared to *Mn*, and *Co*, while *Mn* and *Co* are present in equal amount. Development aims at reducing the content of Cobalt and increasing the
Nickel share. A selection of other chemistry types can be chosen from.


On account that:

-  the battery cycle life values were obtained in the context of an
   experiment [14]_,

-  with loss of 20% of the initial capacity, the battery may still
   provide enough energy to complete the intended route,

cycle life values for NMC and NCA battery chemistries are corrected by
+50%.

**Important assumption**: The environmental burden associated with the
manufacture of spare batteries is entirely allocated to the vehicle use.
The number of battery replacements is rounded up.

Table 7 gives an overview of the number of battery replacements assumed
for the different battery electric vehicles in *carculator*.

Table 7 Number of battery replacements assumed or calculated for each
vehicle type by default

============================= === === ===
\                             NMC LFP NCA
Passenger car, electric, 2020 0   0   0
Passenger car, electric, 2050 0   0   0
============================= === === ===

Users are encouraged to test the sensitivity of end-results on the
number of battery replacements.

The number of battery replacement is calculated as follows:

.. math::

    n_{batt_repl} = \frac{L_{veh}}{L_{batt}} - 1

where :math:`L_{veh}` is the lifetime of the vehicle [km],
and :math:`L_{batt}` is the lifetime of the battery [km].


Liquid and gaseous energy storage
---------------------------------

The oxidation energy stored in liquid fuel tanks is calculated as:

.. math::

    E_{oxid} = m_{fuel} \times E_{lhv}

where :math:`m_{fuel}` is the mass of the fuel [kg], and :math:`E_{lhv}` is the lower heating value of hte fuel [MJ/kg].

The fuel tank mass is calculated as:

.. math::

    m_{tank} = m_{fuel} \times m_{tank unit}

where :math:`m_{tank unit}` being the tank mass per unit of energy [kg/MJ].

Note that the tank mass per unit of energy is different for liquid fuels (gasoline,
diesel), and for gaseous fuels (compressed gas, hydrogen). Also, compressed gas tanks
store at 200 bar, while hydrogen tanks store at 700 bar.

Fuel cell stack
---------------

All fuel cell electric vehicles use a proton exchange membrane
(PEM)-based fuel cell system.

Table 8 lists the specifications of the fuel cell stack and system used
in *carculator* in 2020. The durability of the fuel cell stack,
expressed in hours, is used to determine the number of replacements
needed - the expected kilometric lifetime of the vehicle as well as the
average speed specified by the driving cycle gives the number of hours
of operation. The environmental burden associated with the manufacture
of spare fuel cell systems is entirely allocated to vehicle use as no
reuse channels seem to be implemented for fuel cell stacks at the
moment.

Table 8 Specifications for fuel cell stack systems

+-----------------------+----------+--------+-----------------------+
|                       | 2020     | 2050   | Source                |
+-----------------------+----------+--------+-----------------------+
| Power [kW]            | 65 - 140 | 65 140 | Calculated.           |
+-----------------------+----------+--------+-----------------------+
| Fuel cell stack       | 55-58%   | 60%    | [5]_                   |
| efficiency [%]        |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Fuel cell stack own   | 15%      | 12%    |                       |
| consumption [% of kW  |          |        |                       |
| output]               |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Fuel cell system      | 45-50%   | 53%    |                       |
| efficiency [%]        |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Power density [W/cm2  | 0.9      | 1      | For passenger cars,   |
| cell]                 |          |        | [17]_.                 |
+-----------------------+----------+--------+-----------------------+
| Specific mass [kg     | 0.51     |        |                       |
| cell/W]               |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Platinum loading      | 0.13     |        |                       |
| [mg/cm2]              |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Fuel cell stack       | 4 000    | 5 625  | [18]_ [19]_           |
| durability [hours to  |          |        |                       |
| reach 20% cell        |          |        |                       |
| voltage degradation]  |          |        |                       |
+-----------------------+----------+--------+-----------------------+
| Fuel cell stack       | 1        | 0      | Calculated.           |
| lifetime replacements |          |        |                       |
| [unit]                |          |        |                       |
+-----------------------+----------+--------+-----------------------+

The fuel cell system efficiency :math:`r_{fcsys}` is calculated as:

.. math::

    r_{fcsys} = \frac{r_{fcstack}}{r_{fcown}}

where :math:`r_{fcstack}` is the fuel cell stack efficiency [%], and :math:`r_{fcown}` is the rate of autoconsumption [%].
For reference, the rate of auto-consumption in 2020 for a fuel cell system is 15% (i.e., 15% of the power produced by
the fuel cell system is consumed by it).

The fuel cell system power :math:`P_{fcsys}` is calculated as:

.. math::

    P_{fcsys} = P_{veh} \times r_{fcshare} \times r_{fcown}

where :math:`P_{veh}` is the vehicle engine power and :math:`r_{fcshare}` is the fuel cell system power relative
to the vehicle engine power [%].

Finally, the fuel cell stack mass is calculated as:

.. math::

    m_{fcstack} = 0.51 [kg/kW] \times P_{fcsys} \times \frac{800 [mW/cm^2]}{A_{fc}}

where :math:`P_{fcsys}` is the fuel cell system power [kW],
:math:`A_{fc}` is the fuel cell fuel cell power area density [kW/cm2],
and :math:`m_{fcstack}` is the fuel cell stack mass [kg].

**Important remark:** although fuel cell electric vehicles have a small
battery to enable the recuperation of braking energy, etc., we model
it as a power battery, not a storage battery.
For example, the Toyota Mirai is equipped with a 1.6 kWh
nickel-based battery.

The battery power is calculated as:

.. math::

    P_{batt} = P_{fcsys} \times (1 - r_{fcsys})


where :math:`P_{fcsys}` is the fuel cell system power [kW],
and :math:`r_{fcshare}` is the fuel cell system power share [%].

The number of fuel cell replacements is based on the average distance driven
with a set of fuel cells given their lifetime expressed in hours of use.
The number is replacement is rounded *up* as we assume no allocation of burden
with a second life. It is hence is calculated as:

.. math::

    n_{fcrep} = \frac{L_{veh}}{V_{avg} \times L_{fc}} - 1

where :math:`L_{veh}` is the lifetime of the vehicle [km],
:math:`V_{avg}` is the average speed of the driving cycle selected [km/h],
and :math:`L_{fc}` is the fuel cell lifetime in hours [h].

Light-weighting
---------------

The automotive industry has been increasingly using light weighting
materials to replace steel in engine blocks, chassis, wheels rims and
powertrain components [2]_. However, vehicles light weighting has not led
to an overall curb mass reduction for passenger cars and trucks, as
additional safety equipment compensate for it. According to [20]_,
passenger cars in the EU in 2016 were on average 10% heavier than in
2000.

The dataset used to represent the chassis of passenger cars (i.e.,
"glider, for passenger car") does not reflect today's use of light
weighting materials, such as aluminium and advanced high strength steel
(AHSS).

A report from the Steel Recycling Institute [21]_ indicates that every
kilogram of steel in a car glider can be substituted by 0.75 kilogram of
AHSS or 0.68 kilogram of aluminium. Looking at the material composition
of different car models three years apart, [22]_ show that steel is in
fact increasingly replaced by a combination of both aluminium and AHSS.
However, they also show that the use of AHSS is generally preferred to
aluminium as its mass reduction-to-cost ratio is preferable.

Hence, it is considered that, for a given mass reduction to reach,
two-third of the mass reduction comes from using AHSS, and one third
comes from using aluminium. This means that one kilogram of mass
reduction is achieved by replacing 3.57 kilogram of steel by:

-  1.76 kilogram of AHSS

-  0.8 kilogram of aluminium

Additionally, additional efforts is made to ensure that the final
aluminium content in the chassis corresponds to what is actually found
in current passenger car models, according to [2]_.

While ecoinvent v.3.8 has a LCI dataset for the supply of aluminium, it
is not the case for AHSS. However, an LCA report from the World Steel
Institute [23]_ indicates that AHSS has a similar carbon footprint than
conventional primary low-alloyed steel from a basic oxygen furnace route
(i.e., 2.3 kg CO\ :sub:`2`-eq./kg). We therefore use conventional steel
to represent the use of AHSS.

The amount of light-weighting obtained from the use of light-weighting
materials is:

    :math:`\Delta m_{glider} = m_{glider} \times r_{lightweighting}`

where :math:`\Delta m_{steel}` is the mass reduction of the glider [kg],
and :math:`r_{lightweighting}` is the light-weighting ratio [%].

Sizing of onboard energy storage
--------------------------------

Sizing of battery
_________________

The sizing of batteries for battery electric vehicles is conditioned by
the battery mass, which is defined as an input parameter for each size
class. The battery masses given for the different size classes are
presented in Figure 10 using the battery chemistry NMC, and is based on
representative battery storage capacities available today on the market
- which are represented in relation to the curb mass. The data is
collected from the vehicle's registry of Touring Club Switzerland.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image19.png
   :width: 6.24306in
   :height: 3.25326in

Figure 10 Energy storage capacity for current battery electric cars,
shown in relation to curb mass. Red dots are the energy storage
capacities used for Small, Medium and Large battery electric vehicles in
*carculator*.

Seventy percent of the overall battery mass is assumed to be represented
by the battery cells in the case of NMC and NCA batteries. Given the
energy density of the battery cell considered, this yields the storage
capacity of the battery. A typical depth of discharge of 80% is used to
calculate the available storage capacity.

Table 9 Parameters for battery sizing for battery electric vehicles
using NMC battery chemistry

+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Unit                                                    | Micro     | Small                      | Medium            | Large                                                                         | Large SUV                   |                |
+=========================================================+===========+============================+===================+===============================================================================+=============================+================+
| Storage capacity (reference)                            | kWh       | 14                         | 35                | 45                                                                            | 70                          | 80             |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Commercial models with similar energy storage capacity  |           | Microlino, Renault Twizzy  | VW e-Up!, BMW i3  | Citroen e-C4, DS 3 E.Tense, Peugeot 2008, Peugeot 208, Opel Corsa-e, VW ID.3  | Audi e-Tron, Tesla Model 3  | Jaguar i-Pace  |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Battery mass (system)                                   | Kilogram  | 120                        | 291               | 375                                                                           | 583                         | 660            |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Battery cell mass                                       | %         | ~70%                       |                   |                                                                               |                             |                |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Battery cell mass                                       | Kilogram  | 72                         | 175               | 225                                                                           | 330                         | 400            |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Balance of Plant mass                                   | Kilogram  | 48                         | 116               | 150                                                                           | 233                         | 260            |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Energy density                                          | kWh/kg    | 0.2                        |                   |                                                                               |                             |                |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Storage capacity                                        | kWh       |                            | 35                | 45                                                                            | 70                          | 80             |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Depth of discharge                                      | %         | 80%                        |                   |                                                                               |                             |                |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+
| Storage capacity (available)                            | kWh       | 14                         | 28                | 36                                                                            | 56                          | 65             |
+---------------------------------------------------------+-----------+----------------------------+-------------------+-------------------------------------------------------------------------------+-----------------------------+----------------+

Hence, the battery cell mass is calculated as:

.. math::

    m_{cell} = m_{pack} \times s_{cell}

where :math:`m_{pack}` is the mass of the pack, and :math:`s_{cell}` is the cell-to-pack ratio.

And the electricity stored in the battery is calculated as:

.. math::

    E_{battery} = m_{cell} \times C_{cell}

where :math:`E_{battery}` being battery capacity [kWh], :math:`C_{cell}` is the cell energy density [kg/kWh], and :math:`m_{cell}` is the cell mass [kg].

By deduction, the balance of plant mass is:

.. math::

    m_{BoP} = m_{battery} - m_{cell}

where :math:`m_{battery}` is the mass of the battery [kg], and :math:`m_{cell}` is the cell mass [kg].

Finally, the range autonomy is calculated as:

.. math::

    R_{autonomy} = \frac{C_{battery} \times r_{discharge}}{F_{ttw}}

where :math:`C_{battery}` is the battery capacity [kWh], :math:`r_{discharge}` is the discharge depth [%],
and :math:`F_{ttw}` is the tank-to-wheel energy consumption [kWh/km].


Similarly, plug-in hybrid vehicles are dimensioned to obtain an energy
storage capacity of the battery that corresponds with the capacity of
models available today. The sizing of the battery is similar to what is
described above for battery electric vehicles. The energy storage
capacity of the battery is particularly important for plugin hybrid
vehicles, as it conditions the electric utility factor (the share of
kilometers driven in battery-depleting mode) which calculation is
described in the next section.

Table 10 Parameters for battery sizing for plug-in hybrid vehicles using
NMC battery chemistry

+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+---------------+------------+
|                                                                    | Unit      | Small                | Medium                               | Large         | Large SUV  |
+====================================================================+===========+======================+======================================+===============+============+
| Battery storage capacity (reference)                               | kWh       | 9                    | 13                                   | 18                         |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Commercial models with similar electric and fuel storage capacity  |           | Kia Niro, Kia Xceed  | Skoda Octavia, VW Golf, Cupra Leon   | Suzuki Across, VW Touareg  |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Battery mass (system)                                              | Kilogram  | 80                   | 105                                  | 160                        |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Battery cell mass                                                  | %         | 60%                                                                                      |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Battery cell mass                                                  | Kilogram  | 48                   | 63                                   | 96                         |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Balance of Plant mass                                              | Kilogram  | 32                   | 42                                   | 64                         |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Energy density                                                     | kWh/kg    | 0.2                  |                                                                   |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Battery storage capacity                                           | kWh       | 9                    | 13                                   | 19                         |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Depth of discharge                                                 | %         | 80%                  |                                                                   |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Battery storage capacity (available)                               | kWh       | 7.2                  | 10.4                                 | 15.6                       |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+
| Fuel tank storage capacity                                         | L         | 45                   | 52                                   | 64                         |
+--------------------------------------------------------------------+-----------+----------------------+--------------------------------------+----------------------------+


Note that carculator only considers NMC batteries for plugin hybrid
vehicles.

Electric utility factor
-----------------------

Diesel and gasoline plugin hybrid vehicles are modeled as a composition
of an ICE vehicle and a battery electric vehicle to the extent
determined by the share of km driven in battery-depleting mode (also
called "electric utility factor"). This electric utility factor is
calculated based on a report from the ICCT [24]_, which provides measured
electricity utility factors for 6'000 PHEV *private* owners in Germany
in relation to the vehicle range in battery-depleting mode.

A first step consists in determining the energy consumption of the PHEV
in electric mode as well as its battery size, in order to know its range
autonomy. When the range autonomy is known, the electric utility factor
is interpolated based on the data presented in Table 11.

Table 11 Data points used to interpolate the electric utility factor

+----------------------------------+----------------------------------+
| **Range in battery-depleting     | **Observed electric utility      |
| mode [km]**                      | factor [%]**                     |
+==================================+==================================+
| 20                               | 30                               |
+----------------------------------+----------------------------------+
| 30                               | 41                               |
+----------------------------------+----------------------------------+
| 40                               | 50                               |
+----------------------------------+----------------------------------+
| 50                               | 58                               |
+----------------------------------+----------------------------------+
| 60                               | 65                               |
+----------------------------------+----------------------------------+
| 70                               | 71                               |
+----------------------------------+----------------------------------+
| 80                               | 75                               |
+----------------------------------+----------------------------------+

Once the electric utility factor :math:`U` is known, it is used as a partitioning
ratio to compose the vehicle between the PHEV in combustion mode, and the PHEV in electric mode,
where:

.. math::

    F_{ttw_phev} = (F_{ttw_phev_e} \times U) + (F_{ttw_phev_c} \times (1 - U))

.. math::

    m_{curb_phev} = (m_{curb_phev_e} \times U) + (m_{curb_phev_c} \times (1 - U))

where :math:`F_{ttw_phev_e}` is the tank-to-wheel energy consumption [kWh/km] of the electric PHEV,
:math:`F_{ttw_phev_c}` is the tank-to-wheel energy consumption [kk/km] of the combustion PHEV,
:math:`m_{curb_phev_e}` is the curb weight [kg] of the electric PHEV, and :math:`m_{curb_phev_c}` is the
curb weight [kg] of the combustion PHEV.



Inventory modelling
*******************

Once the vehicles are modeled, the calculated parameters of each of them
is passed to the inventory.py calculation module to derive inventories.
When the inventories for the vehicle and the transport are calculated,
they can be normalized by the kilometric lifetime (i.e., vehicle-kilometer)
or by the kilometric multiplied by the passenger occupancy (i.e., passenger-kilometer).

Road demand
-----------

The demand for construction and maintenance of roads and road-related
infrastructure is calculated on the following basis:

-  Road construction: 5.37e-7 meter-year per kg of vehicle mass per km.

-  Road maintenance: 1.29e-3 meter-year per km, regardless of vehicle
   mass.

The driving mass of the vehicle consists of the mass of the vehicle in
running condition (including fuel) in addition to the mass of passengers
and cargo, if any. Unless changed, the passenger mass is 75 kilograms,
and the average occupancy is 1.6 persons per vehicle.

The demand rates used to calculate the amounts required for road
construction and maintenance (based on vehicle mass per km and per km,
respectively) are taken from [25]_.

Because roads are maintained by removing surface layers older than those
that are actually discarded, road infrastructure disposal is modeled in
ecoinvent as a renewal rate over the year in the road construction
dataset.

Fuel properties
---------------

For all vehicles with an internal combustion engine, carbon dioxide
(CO\ :sub:`2`) and sulfur dioxide (SO\ :sub:`2`) emissions are
calculated based on the fuel consumption of the vehicle and the carbon
and sulfur concentration of the fuel observed in Switzerland and Europe.
Sulfur concentration values are sourced from HBEFA 4.1 [26]_. Lower
heating values and CO\ :sub:`2` emission factors for fuels are sourced
from p.86 and p.103 of [27]_. The fuel properties shown in Table 12 are
used for fuels purchased in Switzerland.

Table 12 Fuels characteristics

+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
|                                       | Volumetric mass density [kg/l]  | Lower heating value [MJ/kg]  | CO2 emission factor [kg CO2/kg]  | SO2 emission factor [kg SO2/kg]  |
+=======================================+=================================+==============================+==================================+==================================+
| Gasoline                              | 0.75                            | 42.6                         | 3.14                             | 1.6e-5                           |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Bioethanol                            | 0.75                            | 26.5                         | 1.96                             | 1.6e-5                           |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Synthetic gasoline                    | 0.75                            | 43                           | 3.14                             | 0                                |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Diesel                                | 0.85                            | 43                           | 3.15                             | 8.85e-4                          |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Biodiesel                             | 0.85                            | 38                           | 2.79                             | 8.85e-4                          |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Synthetic diesel                      | 0.85                            | 43                           | 3.15                             | 0                                |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Natural gas                           |                                 | 47.5                         | 2.68                             |                                  |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Bio-methane                           |                                 | 47.5                         | 2.68                             |                                  |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+
| Synthetic methane                     |                                 | 47.5                         | 2.68                             |                                  |
+---------------------------------------+---------------------------------+------------------------------+----------------------------------+----------------------------------+


Because large variations are observed in terms of sulfur concentration
in biofuels, similar values than that of conventional fuels are used.

Exhaust emissions
-----------------

Emissions of regulated and non-regulated substances during driving are
approximated using emission factors from HBEFA 4.1 [26]_. Emission
factors are typically given in gram per km. Emission factors
representing free flowing driving conditions and urban and rural traffic
situations are used. Additionally, cold start emissions as well as
running, evaporation and diurnal losses are accounted for, also sourced
from HBEFA 4.1 [26]_.


For vehicles with an internal combustion engine, the sulfur
concentration values in the fuel can slightly differ across regions -
although this remains rather limited within Europe. The values provided
by HBEFA 4.1 are used for Switzerland, France, Germany, Austria and
Sweden. For other countries, values from [28]_ are used.

Table 13 Sulfur concentration values examples for on-road fuel in
Switzerland and average Europe

========================= =============== ==========
**Sulfur [ppm/fuel wt.]** **Switzerland** **Europe**
========================= =============== ==========
Gasoline                  8               8
Diesel                    10              8
========================= =============== ==========

The amount of sulfur dioxide released by the vehicle over one km [kg/km] is calculated as:

.. math::

        SO_2 = r_{S} \times F_{fuel} \times (64/32)

where :math:`r_{S}` is the sulfur content per kg of fuel [kg SO2/kg fuel],
:math:`F_{fuel}` is the fuel consumption of the vehicle [kg/km],
and :math:`64/32` is the ratio between the molar mass of SO2 and the molar mass of O2.

Country-specific fuel blends are sourced from the IEA's Extended World
Energy Balances database [29]_. By default, the biofuel used is assumed
to be produced from biomass residues (i.e., second-generation fuel):
fermentation of crop residues for bioethanol, esterification of used
vegetable oil for biodiesel and anaerobic digestion of sewage sludge for
bio-methane.

Table 14 Specification examples of fuel blends for Switzerland and
average Europe

========================= =============== ==========
**Biofuel share [% wt.]** **Switzerland** **Europe**
========================= =============== ==========
Gasoline blend            1.2             4
Diesel blend              4.8             6
Compressed gas blend      22              9
========================= =============== ==========

A number of fuel-related emissions other than CO\ :sub:`2` and
SO\ :sub:`2` are considered, using the HBEFA 4.1 database [30]_.

Six sources of emissions are considered:

-  Exhaust emissions: emissions from the combustion of fuel during
   operation. Their concentration relates to the fuel consumption and
   the emission standard of the vehicle.

-  Cold start emissions: emissions when starting the engine. The factor
   is given in grams per engine start. 2.3 engine starts per day are
   considered [27]_ and an annual mileage of 12'000 km.

-  Diurnal emissions: evaporation of the fuel due to a temperature
   increase of the vehicle. The factor is given in grams per day.
   Emissions are distributed evenly along the driving cycle, based on an
   annual mileage of 12'000 km per year.

-  Hot soak emissions: evaporative emissions occurring after the vehicle
   has been used. The factor is given in grams per trip. The emission is
   added at the end of the driving cycle.

-  In addition, running loss emissions: emissions related to the
   evaporation of fuel (i.e., not combusted) during operation. The
   factor is given in grams per km. Emissions are distributed evenly
   along the driving cycle.

-  Other non-exhaust emissions: brake, tire road wear and re-suspended
   road dust emissions, as well as emissions of refrigerant.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image20.png
   :width: 4.38333in
   :height: 2.325in

Figure 11 Representation of the different sources of emission other than
exhaust emissions

For exhaust emissions, factors based on the fuel consumption are derived
by comparing emission data points for different traffic situations
(i.e., grams emitted per vehicle-km) for in a free flowing driving
situation, with the fuel consumption corresponding to each data point
(i.e., MJ of fuel consumed per km), as illustrated in Figure 12 for a
diesel-powered engine. The aim is to obtain emission factors expressed
in grams of substance emitted per MJ of fuel consumed, to be able to
model emissions of passenger cars of different sizes and fuel efficiency
and for different driving cycles.

Hence, the emission of substance i at second s of the driving cycle is
calculated as follows:

.. math::

    E(i,s) = F_ttw(s) \times X(i, e)

where :math:`E(i,s)` is the emission of substance i at second s of the driving cycle,
:math:`F_ttw(s)` is the fuel consumption of the vehicle at second s,
and :math:`X(i, e)` is the emission factor of substance i in the given driving conditions.

To that, we add the following terms:

    - Cold start emissions on the first second of the driving cycle
    - Evaporation emissions: on the last second of the driving cycle
    - Diurnal and running losses: distributed evenly over the driving cycle


**Important remark**: the degradation of anti-pollution systems for
diesel and gasoline cars (i.e., catalytic converters) is accounted for
as indicated by HBEFA, by applying a degradation factor on the emission
factors for CO, HC and NO\ :sub:`x` for gasoline cars, as well as on CO
and NO\ :sub:`x` for diesel cars. These factors are shown in Table 15
for passenger cars with a mileage of 200'000 km, which is the default
lifetime value in *carculator*. The degradation factor corresponding to
half of the vehicle kilometric lifetime is used, to obtain a
lifetime-weighted average degradation factor.

Table 15 Degradation factors at 200'000 km for passenger cars

+-----------------------------------+--------------------------+-------+------+------------------------+------+
| Degradation factor at 200 000 km  | Gasoline passenger cars  |       |      | Diesel passenger cars  |      |
+===================================+==========================+=======+======+========================+======+
|                                   | CO                       | HC    | NOx  | CO                     | NOx  |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-1                            | 1.9                      | 1.59  | 2.5  |                        |      |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-2                            | 1.6                      | 1.59  | 2.3  |                        | 1.25 |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-3                            | 1.75                     | 1.02  | 2.9  |                        | 1.2  |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-4                            | 1.9                      | 1.02  | 2    | 1.3                    | 1.06 |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-5                            | 2                        |       | 2.5  | 1.3                    | 1.03 |
+-----------------------------------+--------------------------+-------+------+------------------------+------+
| EURO-6                            | 1.3                      |       | 1.3  | 1.4                    | 1.15 |
+-----------------------------------+--------------------------+-------+------+------------------------+------+


.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image21.png
   :width: 6.27014in
   :height: 7.84756in

Figure 12 Relation between emission factor and fuel consumption for a
diesel-powered passenger car. Dots represent HBEFA 4.1 emission factors
for different traffic situation for a diesel engine, for different
emission standards.

However, as Figure 12 shows, the relation between amounts emitted and
fuel consumption is not always obvious and using a linear relation
between amounts emitted and fuel consumption can potentially be
incorrect. In addition, emissions of ammonia (NH\ :sub:`3`) and Nitrous
oxides (N\ :sub:`2`\ O) seem to be related to the emission standard
(e.g., use of urea solution) and engine temperature rather than the fuel
consumption.

To confirm that such approach does not yield kilometric emissions too
different from the emission factors per vehicle-kilometer proposed by
HBEFA 4.1, Figure 13 compares the emissions obtained by *carculator*
using the WLTC driving cycle over 1 vehicle-km (red dots) with the
distribution of the emission factors for different traffic situations
(green box-and-whiskers) as well as the traffic situation-weighted
average emission factor (yellow dots) given by HBEFA 4.1 for different
emission standards for a medium diesel-powered passenger car.

There is some variation across traffic situations, but the emissions
obtained remain, for most substances, within the 50% of the distributed
HBEFA values across traffic situations. Also, the distance between the
modeled emission and the traffic-situation-weighted average is
reasonable.

**Important remark**: NO\ :sub:`x` emissions for emission standards
EURO-4 and 5 tend to be under-estimated compared to HBEFA's values. It
is also important to highlight that, in some traffic situations, HBEFA's
values show that emissions of CO, HC, NMHC and PMs for vehicles with
early emission standards can be much higher that what is assumed in
*carculator*. There is overall a good agreement between traffic
situation-weighted average emission factors and those used in
*carculator*.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image22.png
   :width: 6.27014in
   :height: 5.97153in

Figure 13 Validation of the exhaust emissions model with the emission
factors provided by HBEFA 4.1 for a medium size diesel-powered passenger
car. Box-and-whiskers: distribution of HBEFA's emission factors for
different traffic situations (box: 50% of the distribution, whiskers:
90% of the distribution). Yellow dots: traffic situation-weighted
average emission factors. Red dots: modeled emissions calculated by
*carculator* with the WLTC cycle, using the relation between fuel
consumption and amounts emitted.

NMHC speciation
_______________

After NMHC emissions are quantified, EEA/EMEP's 2019 Air Pollutant
Emission Inventory Guidebook provides factors to further specify some of
them into the substances listed in Table 16.

Table 16 NMVOC sub-species as fractions of the mass emitted

=================== ========================= =======================
\                   **All gasoline vehicles** **All diesel vehicles**
\                   *Wt. % of NMVOC*          *Wt. % of NMVOC*
Ethane              3.2                       0.33
Propane             0.7                       0.11
Butane              5.2                       0.11
Pentane             2.2                       0.04
Hexane              1.6                       0
Cyclohexane         1.1                       0.65
Heptane             0.7                       0.2
Ethene              7.3                       10.97
Propene             3.8                       3.6
1-Pentene           0.1                       0
Toluene             11                        0.69
m-Xylene            5.4                       0.61
o-Xylene            2.3                       0.27
Formaldehyde        1.7                       12
Acetaldehyde        0.8                       6.47
Benzaldehyde        0.2                       0.86
Acetone             0.6                       2.94
Methyl ethyl ketone 0.1                       1.2
Acrolein            0.2                       3.58
Styrene             1                         0.37
NMVOC, unspecified  50.8                      55
=================== ========================= =======================

Non-exhaust emissions
---------------------

A number of emission sources besides exhaust emissions are considered.
They are described in the following sub-sections.

Engine wear emissions
_____________________

Metals and other substances are emitted during the combustion of fuel
because of engine wear. These emissions are scaled based on the fuel
consumption, using the emission factors listed in Table 17, sourced from
[31]_.

Table 17 Emission factors for engine wear as fractions of the fuel mass
combusted

=========== ========================= =======================
\           **All gasoline vehicles** **All diesel vehicles**
\           *kg/MJ fuel*              *kg/MJ fuel*
PAH         8.19E-10                  1.32E-09
Arsenic     7.06E-12                  2.33E-12
Selenium    4.71E-12                  2.33E-12
Zinc        5.08E-08                  4.05E-08
Copper      9.88E-10                  4.93E-10
Nickel      3.06E-10                  2.05E-10
Chromium    3.76E-10                  6.98E-10
Chromium VI 7.53E-13                  1.40E-12
Mercury     2.05E-10                  1.23E-10
Cadmium     2.54E-10                  2.02E-10
=========== ========================= =======================

Abrasion emissions
__________________

We distinguish four types of abrasion emissions, besides engine wear
emissions:

-  brake wear emissions: from the wearing out of brake drums, discs and
   pads

-  tires wear emissions: from the wearing out of rubber tires on the
   asphalt

-  road wear emissions: from the wearing out of the road pavement

and re-suspended road dust: dust on the road surface that is
re-suspended as a result of passing traffic, "due either to shear forces
at the tire/road surface interface, or air turbulence in the wake of a
moving vehicle" [32]_.

[32]_ provides an approach for estimating the mass and extent of these
abrasion emissions. They propose to disaggregate the abrasion emission
factors presented in the EMEP's 2019 Emission inventory guidebook [31]_
for two-wheelers, passenger cars, buses and heavy good vehicles, to
re-quantify them as a function of vehicle mass, but also traffic
situations (urban, rural and motorway). Additionally, they present an
approach to calculate re-suspended road dust according to the method
presented in [33]_ - such factors are not present in the EMEP's 2019
Emission inventory guidebook - using representative values for dust load
on European roads.

The equation to calculate brake, tire, road and re-suspended road dust
emissions is the following:

.. math::

    EF=b.W^{\frac{1}{c}}

With:

-  :math:`EF` being the emission factor, in mg per vehicle-kilometer

-  :math:`W` being the vehicle mass, in tons

-  :math:`b` and :math:`c` being regression coefficients, whose values are presented
   in Table 18.

Table 18 Regression coefficients to estimate abrasion emissions

+--------+------------+------+--------+------+-----------+------+-------------+------+--------+------+-----------+------+------------+------+-------------------------+------+
|        | Tire wear  |      |        |      |           |      | Brake wear  |      |        |      |           |      | Road wear  |      | Re-suspended road dust  |      |
+========+============+======+========+======+===========+======+=============+======+========+======+===========+======+============+======+=========================+======+
|        | Urban      |      | Rural  |      | Motorway  |      | Urban       |      | Rural  |      | Motorway  |      |            |      |                         |      |
+--------+------------+------+--------+------+-----------+------+-------------+------+--------+------+-----------+------+------------+------+-------------------------+------+
|        | b          | c    | b      | c    | b         | c    | b           | c    | b      | c    | b         | c    | b          | c    | b                       | c    |
+--------+------------+------+--------+------+-----------+------+-------------+------+--------+------+-----------+------+------------+------+-------------------------+------+
| PM 10  | 5.8        | 2.3  | 4.5    | 2.3  | 3.8       | 2.3  | 4.2         | 1.9  | 1.8    | 1.5  | 0.4       | 1.3  | 2.8        | 1.5  | 2                       | 1.1  |
+--------+------------+------+--------+------+-----------+------+-------------+------+--------+------+-----------+------+------------+------+-------------------------+------+
| PM 2.5 | 8.2        | 2.3  | 6.4    | 2.3  | 5.5       | 2.3  | 11          | 1.9  | 4.5    | 1.5  | 1         | 1.3  | 5.1        | 1.5  | 8.2                     | 1.1  |
+--------+------------+------+--------+------+-----------+------+-------------+------+--------+------+-----------+------+------------+------+-------------------------+------+


The respective amounts of brake and tire wear emissions in urban, rural
and motorway driving conditions are weighted, to represent the driving
cycle used. The weight coefficients sum to 1 and the coefficients
considered are presented in Table *19*. They have been calculated by
analyzing the speed profile of each driving cycle, with the exception of
two-wheelers, for which no driving cycle is used (i.e., the energy
consumption is from reported values) and where simple assumptions are
made in that regard instead.

Table 19 Weighting coefficients to calculate representative abrasion
emissions given a type of use/driving cycle

============= ============= ===== ===== ========
\             Driving cycle Urban Rural Motorway
============= ============= ===== ===== ========
Passenger car WLTP          0.33  0.24  0.43
============= ============= ===== ===== ========

Finally, for electric and (plugin) hybrid vehicles, the amount of brake
wear emissions is reduced. This reduction is calculated as the ratio
between the sum of energy recuperated by the regenerative braking system
and the sum of negative resistance along the driving cycle. The logic is
that the amount of negative resistance that could not be met by the
regenerative braking system needs to be met with mechanical brakes. This
is illustrated in Figure 14, where the distance between the recuperated
energy and the total negative motive energy corresponds to the amount of
energy that needs to be provided by mechanical brakes. Table 20 lists
such reduction actors for the different powertrains.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image23.png
   :width: 4.1546in
   :height: 2.85in

Figure 14 Negative motive energy and recuperated energy between second
300 and 450 of the WLTC driving cycle

Table 20 Approximate reduction factors for brake wear emissions. Values
differ slightly across size classes.

+-------------+-------------+-------------+-------------+-------------+
|             | Driving     | Reduction   | Reduction   | Reduction   |
|             | cycle       | factor for  | factor for  | factor for  |
|             |             | hybrid      | plugin      | battery and |
|             |             | vehicles    | hybrid      | fuel cell   |
|             |             |             | vehicles    | electric    |
|             |             |             |             | vehicles    |
+=============+=============+=============+=============+=============+
| Passenger   | WLTP        | -72%        | -73%        | -76%        |
| car         |             |             |             |             |
+-------------+-------------+-------------+-------------+-------------+

The sum of PM 2.5 and PM 10 emissions is used as the input for the
ecoinvent v.3.x LCI datasets indicated in Table 21.

Table 21 LCI datasets used to approximate PM emissions composition and
emissions to air, soil and water

+-------------+-------------+-------------+-------------+-------------+
|             | Tire wear   | Brake wear  | Road wear   | R           |
|             |             |             |             | e-suspended |
|             |             |             |             | road dust   |
+=============+=============+=============+=============+=============+
| Passenger   | Tire wear   | Brake wear  | Road wear   |             |
| car         | emissions,  | emissions,  | emissions,  |             |
|             | passenger   | passenger   | passenger   |             |
|             | car         | car         | car         |             |
+-------------+-------------+-------------+-------------+-------------+

Finally, we assume that the composition of the re-suspended road dust is
evenly distributed between brake, road and tire wear particles.

Figure 15 shows the calculated abrasion emissions for passenger cars in
mg per vehicle-kilometer, following the approach presented above.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image24.png
   :width: 6.26667in
   :height: 3.00136in

Figure 15 Total particulate matter emissions (<2.5 µm and 2.5-10 µm) in
mg per vehicle-kilometer for passenger cars.

Re-suspended road dust emissions are assumed to be evenly composed of
brake wear (33.3%), tire wear (33.3%) and road wear (33.3%) particles.

Refrigerant emissions
_____________________

The use of refrigerant for onboard air conditioning systems is
considered for passenger cars until 2021. The supply of refrigerant gas R134a is
accounted for. Similarly, the leakage of the refrigerant is also
considered. For this, the calculations from [34]_ are used. Such emission
is included in the transportation dataset of the corresponding vehicle.
The overall supply of refrigerant amounts to the initial charge plus the
amount leaked throughout the lifetime of the vehicle, both listed in
Table 22. This is an important aspect, as the refrigerant gas R134a has
a Global Warming potential of 2'400 kg CO\ :sub:`2`-eq./kg released in
the atmosphere.

Table 22 Use and loss of refrigerant gas for onboard air conditioning
systems

======================================== =============================
\                                        Passenger cars (except Micro)
Initial charge [kg per vehicle lifetime] 0.55
Lifetime loss [kg per vehicle lifetime]  0.75
======================================== =============================

**Important assumption**: it is assumed that electric and plug-in
electric vehicles also use a compressor-like belt-driven air
conditioning system, relying on the refrigerant gas R134a. In practice,
an increasing, but still minor, share of electric vehicles now use a
(reversible) heat pump to provide cooling.

**Important remark:** Micro cars do not have an air conditioning system.
Hence, no supply or leakage of refrigerant is considered for those.

**Important remark:** After 2021, R134a is no longer used.

Noise emissions
---------------

Noise emissions along the driving cycle of the vehicle are quantified
using the method developed within the CNOSSOS project [35]_, which are
expressed in joules, for each of the 8 octaves. Rolling and propulsion
noise emissions are quantified separately.

The sound power level of rolling noise is calculated using:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image25.png
   :width: 3.45in
   :height: 0.65in

With:

-  *V\ m* being the instant speed given by the driving cycle, in km/h

-  *V\ ref* being the reference speed of 70 km/h

And *A\ R,i,m* and *B\ R,i,m*\ are unitless and given in Table 23.

The propulsion noise level is calculated using:

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image26.png
   :width: 3.6in
   :height: 0.625in

With:

And *A\ P,i,m* and *B\ P,i,m*\ are unitless and given in Table 23.

Table 23 Noise level coefficients for passenger cars

================================= ====== ====== ====== ======
Octave band center frequency (Hz) *A\ R* *B\ R* *A\ P* *B\ P*
================================= ====== ====== ====== ======
63                                84     30     101    -1.9
125                               88.7   35.8   96.5   4.7
250                               91.5   32.6   98.8   6.4
500                               96.7   23.8   96.8   6.5
1000                              97.4   30.1   98.6   6.5
2000                              90.9   36.2   95.2   6.5
4000                              83.8   38.3   88.8   6.5
8000                              80.5   40.1   82.7   6.5
================================= ====== ====== ====== ======

A correction factor for battery electric and fuel cell electric vehicles
is applied, and is sourced from [36]_. Also, electric vehicles are added
a warning signal of 56 dB at speed levels below 20 km/h. Finally, hybrid
vehicles are assumed to use an electric engine up to a speed level of 30
km/h, beyond which the combustion engine is used.

The total noise level (in A-weighted decibels) is calculated using the
following equation:

.. math:: L_{W,\ dBA} = 10*\log\left( 10^{\frac{L_{W,R}}{10}} \right) + 10*log(10^{\frac{L_{W,P}}{10}})

The total sound power level is converted into Watts (or joules per
second), using the following equation:

.. math:: L_{W} = \ 10^{- 12}*10^{\frac{L_{W,\ dBA}}{10}}

The total sound power, for each second of the driving cycle, is then
distributed between the urban, suburban and rural inventory emission
compartments.

Figure 16.a illustrates a comparison of noise levels between an ICEV and
BEV as calculated by the tool, over the driving cycle WLTC. In this
figure, the noise levels at different frequency ranges have been summed
together to obtain a total noise level (in dB), and converted to dB(A)
using the A-weighting correction factor, to better represent the
"loudness" or discomfort to the human ear. Typically, propulsion noise
emissions dominate in urban environments (which corresponds to the
section 3.1 of the driving cycle), thereby justifying the use of
electric vehicles in that regard. This is represented by the difference
between the ICEV and BEV lines in the section 3.1 of the driving cycle
in Figure 16.a. The difference in noise level between the two
powertrains diminishes at higher speed levels (see sections 3.2, 3.3 and
3.4) as rolling noise emissions dominate above a speed level of
approximately 50 km/h. This can be seen in Figure 16.b, which sums up
the sound energy produced, in joules, over the course of the driving
cycle.

The study from Cucurachi and Heijungs [37]_ provides compartment-specific
noise emission characterization factors against midpoint and endpoint
indicators - expressed in Person-Pascal-second and Disability-Adjusted
Life Year, respectively.

.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image27.png
   :width: 60%
.. image:: https://github.com/romainsacchi/carculator/raw/master/docs/image28.png
   :width: 30%

Figure 16 a) Noise emission level comparison between ICEV and BEV, based
on the driving cycle WLTC. b) Summed sound energy comparison between
ICEV, BEV and PHEV, over the duration of the WLTC driving cycle.

Electricity mix calculation
---------------------------

Electricity supply mix are calculated based on the weighting from the
distribution the lifetime kilometers of the vehicles over the years of
use. For example, should a BEV enter the fleet in Poland in 2020, most
LCA models of passenger vehicles would use the electricity mix for
Poland corresponding to that year, which corresponds to the row of the
year 2020 in Table 24, based on ENTSO-E's TYNDP 2020 projections
(National Trends scenario) [38]_. *carculator* calculates instead the
average electricity mix obtained from distributing the annual kilometers
driven along the vehicle lifetime, assuming an equal number of
kilometers is driven each year. Therefore, with a lifetime of 200,000 km
and an annual mileage of 12,000 kilometers, the projected electricity
mixes to consider between 2020 and 2035 for Poland are shown in Table
24. Using the kilometer-distributed average of the projected mixes
between 2020 and 2035 results in the electricity mix presented in the
last row of Table 24. The difference in terms of technology contribution
and unitary GHG-intensity between the electricity mix of 2020 and the
electricity mix based on the annual kilometer distribution is
significant (-23%). The merit of this approach ultimately depends on
whether the projections will be realized or not.

It is also important to remember that the unitary GHG emissions of each
electricity-producing technology changes over time, as the background
database ecoinvent has been transformed by premise [39]_: for example,
photovoltaic panels become more efficient, as well as some of the
combustion-based technologies (e.g., natural gas). For more information
about the transformation performed on the background life cycle
database, refer to [39]_.

Table 24 Example of calculation of the carbon intensity of a
km-distributed electricity supply mix for Poland, along with the per kWh
GHG-intensity, for a vehicle first driven in 2020 and driven for the
next 16 years.

+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| year  | Biomass  | Coal  | Gas  | Gas CCGT  | Gas CHP  | Hydro  | Hydro, reservoir  | Lignite  | Nuclear  | Oil  | Solar  | Waste  | Wind  | Wind, offshore  | g CO2-eq./kWh  |
+=======+==========+=======+======+===========+==========+========+===================+==========+==========+======+========+========+=======+=================+================+
| 2020  | 3%       | 46%   | 2%   | 3%        | 0%       | 3%     | 1%                | 29%      | 3%       | 0%   | 0%     | 0%     | 9%    | 0%              | 863            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2021  | 2%       | 43%   | 2%   | 4%        | 1%       | 3%     | 1%                | 29%      | 2%       | 0%   | 1%     | 3%     | 9%    | 0%              | 841            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2022  | 2%       | 41%   | 1%   | 5%        | 1%       | 3%     | 1%                | 28%      | 2%       | 0%   | 2%     | 5%     | 9%    | 0%              | 807            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2023  | 1%       | 38%   | 1%   | 5%        | 2%       | 2%     | 1%                | 28%      | 1%       | 0%   | 3%     | 8%     | 10%   | 0%              | 781            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2024  | 1%       | 36%   | 0%   | 6%        | 2%       | 2%     | 0%                | 27%      | 1%       | 0%   | 3%     | 11%    | 10%   | 0%              | 745            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2025  | 0%       | 33%   | 0%   | 7%        | 3%       | 2%     | 0%                | 27%      | 0%       | 0%   | 4%     | 13%    | 10%   | 0%              | 724            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2026  | 0%       | 31%   | 0%   | 8%        | 3%       | 2%     | 0%                | 25%      | 0%       | 0%   | 5%     | 13%    | 11%   | 2%              | 684            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2027  | 0%       | 28%   | 0%   | 9%        | 4%       | 2%     | 0%                | 24%      | 0%       | 0%   | 6%     | 12%    | 12%   | 3%              | 652            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2028  | 0%       | 25%   | 0%   | 9%        | 5%       | 2%     | 0%                | 23%      | 0%       | 0%   | 6%     | 12%    | 13%   | 5%              | 614            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2029  | 0%       | 23%   | 0%   | 10%       | 6%       | 2%     | 0%                | 21%      | 0%       | 0%   | 7%     | 11%    | 14%   | 6%              | 580            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2030  | 0%       | 20%   | 0%   | 11%       | 6%       | 2%     | 0%                | 20%      | 0%       | 0%   | 8%     | 10%    | 15%   | 8%              | 542            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2031  | 0%       | 19%   | 0%   | 11%       | 7%       | 2%     | 0%                | 18%      | 1%       | 0%   | 9%     | 10%    | 16%   | 8%              | 514            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2032  | 0%       | 17%   | 0%   | 10%       | 8%       | 2%     | 0%                | 16%      | 3%       | 0%   | 9%     | 9%     | 17%   | 9%              | 470            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2033  | 0%       | 16%   | 0%   | 10%       | 8%       | 2%     | 0%                | 14%      | 4%       | 0%   | 10%    | 8%     | 17%   | 10%             | 437            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2034  | 0%       | 15%   | 0%   | 10%       | 9%       | 2%     | 0%                | 12%      | 5%       | 0%   | 10%    | 8%     | 18%   | 11%             | 408            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| 2035  | 0%       | 13%   | 0%   | 9%        | 10%      | 2%     | 0%                | 11%      | 7%       | 0%   | 11%    | 7%     | 19%   | 12%             | 377            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+
| Mix   | 0%       | 26%   | 0%   | 7%        | 5%       | 2%     | 0%                | 21%      | 2%       | 0%   | 6%     | 8%     | 13%   | 5%              | 668            |
+-------+----------+-------+------+-----------+----------+--------+-------------------+----------+----------+------+--------+--------+-------+-----------------+----------------+


Inventories for fuel pathways
-----------------------------

A number of inventories for fuel production and supply are used by
*carculator*. They represent an update in comparison to the inventories
used in the passenger vehicles model initially published by Cox et
al.[5]_. The fuel pathways presented in Table 25 are from the literature
and not present as generic ecoinvent datasets.

+-----------+---------------------------+---------------------------+
| Author(s) | Fuel type                 | Description               |
+===========+===========================+===========================+
| [40]_     | Bioethanol from forest    | Biofuels made from        |
|           | residues                  | biomass residues (e.g.,   |
|           |                           | wheat straw, corn starch) |
|           |                           | or energy crops (e.g.,    |
|           |                           | sugarbeet). For energy    |
|           |                           | crops biofuels, indirect  |
|           |                           | land use change is        |
|           |                           | included.                 |
+-----------+---------------------------+---------------------------+
|           | Bioethanol from wheat     |                           |
|           | straw                     |                           |
+-----------+---------------------------+---------------------------+
|           | Bioethanol from corn      |                           |
|           | starch                    |                           |
+-----------+---------------------------+---------------------------+
|           | Bioethanol from sugarbeet |                           |
+-----------+---------------------------+---------------------------+
| [41]_     | e-Gasoline                | Gasoline produced from    |
|           | (Methanol-to-Gasoline)    | methanol, via a           |
|           |                           | Methanol-to-Gasoline      |
|           |                           | process. The carbon       |
|           |                           | monoxide is provided by a |
|           |                           | reverse water gas shift   |
|           |                           | process, feeding on       |
|           |                           | carbon dioxide from       |
|           |                           | direct air capture. In    |
|           |                           | carculator, one can       |
|           |                           | choose the nature of the  |
|           |                           | heat needed for the       |
|           |                           | methanol distillation as  |
|           |                           | well as for regenerating  |
|           |                           | the DAC sorbent: natural  |
|           |                           | gas, waste heat, biomass  |
|           |                           | heat, or market heat      |
|           |                           | (i.e., a mix of natural   |
|           |                           | gas and fuel oil).        |
+-----------+---------------------------+---------------------------+
| [40]_     | Biodiesel from            | 2\ :sup:`nd` and          |
|           | micro-algae               | 3\ :sup:`rd` generation   |
|           |                           | biofuels made from        |
|           |                           | biomass residues or       |
|           |                           | algae.                    |
+-----------+---------------------------+---------------------------+
|           | Biodiesel from used       |                           |
|           | cooking oil               |                           |
+-----------+---------------------------+---------------------------+
| [42]_     | e-Diesel                  | Diesel produced from      |
|           | (Fischer-Tropsch)         | "blue crude" via a        |
|           |                           | Fischer-Tropsch process.  |
|           |                           | The H\ :sub:`2` is        |
|           |                           | produced via              |
|           |                           | electrolysis, while the   |
|           |                           | CO\ :sub:`2` comes from   |
|           |                           | direct air capture. Note  |
|           |                           | that in *carculator*, two |
|           |                           | allocation approaches at  |
|           |                           | the crude-to-fuel step    |
|           |                           | are possible between the  |
|           |                           | different co-products     |
|           |                           | (i.e., diesel, naphtha,   |
|           |                           | wax oil, kerosene):       |
|           |                           | energy or economic.       |
+-----------+---------------------------+---------------------------+
| [43]_     | Biomethane from sewage    | Methane produced from the |
|           | sludge                    | anaerobic digestion of    |
|           |                           | sewage sludge. The biogas |
|           |                           | is upgraded to biomethane |
|           |                           | (the CO\ :sub:`2` is      |
|           |                           | separated and vented out) |
|           |                           | to a vehicle grade        |
|           |                           | quality.                  |
+-----------+---------------------------+---------------------------+
|           | Synthetic methane         | Methane produced via an   |
|           |                           | electrochemical           |
|           |                           | methanation process, with |
|           |                           | H\ :sub:`2` from          |
|           |                           | electrolysis and          |
|           |                           | CO\ :sub:`2` from direct  |
|           |                           | air capture.              |
+-----------+---------------------------+---------------------------+
| [44, 45]_ | Hydrogen from             | The electricity           |
|           | electrolysis              | requirement to operate    |
|           |                           | the electrolyzer changes  |
|           |                           | over time: from 58 kWh    |
|           |                           | per kg of H\ :sub:`2` in  |
|           |                           | 2010, down to 44 kWh in   |
|           |                           | 2050, according to [46]_. |
+-----------+---------------------------+---------------------------+
| [45, 47]_ | Hydrogen from Steam       | Available for natural gas |
|           | Methane Reforming         | and biomethane, with and  |
|           |                           | without Carbon Capture    |
|           |                           | and Storage (CCS).        |
+-----------+---------------------------+---------------------------+
| [44]_     | Hydrogen from woody       | Available with and        |
|           | biomass gasification      | without Carbon Capture    |
|           |                           | and Storage (CCS).        |
+-----------+---------------------------+---------------------------+

Table 25 List of inventories for different fuel types

Inventories for energy storage components
-----------------------------------------

The source for the inventories used to model energy storage components
are listed in Table 26.

+-----------+---------------------------+---------------------------+
| Author(s) | Energy storage type       | Description               |
+===========+===========================+===========================+
| [48,49]_  | NMC-111/622/811 battery   | Originally from [48]_,    |
|           |                           | then updated and          |
|           |                           | integrated in ecoinvent   |
|           |                           | v.3.8 (with some errors), |
|           |                           | corrected and integrated  |
|           |                           | in *carculator*.          |
|           |                           | Additionally, these       |
|           |                           | inventories relied        |
|           |                           | exclusively on synthetic  |
|           |                           | graphite. This is has too |
|           |                           | been modified: the anode  |
|           |                           | production relies on a    |
|           |                           | 50:50 mix of natural and  |
|           |                           | synthetic graphite, as it |
|           |                           | seems to be the current   |
|           |                           | norm in the industry      |
|           |                           | [50]_. Inventories for    |
|           |                           | natural graphite are from |
|           |                           | [51]_.                    |
+-----------+---------------------------+                           |
|           | NCA battery               |                           |
+-----------+---------------------------+                           |
|           | LFP battery               |                           |
+-----------+---------------------------+---------------------------+
| [52]_     | Type IV hydrogen tank,    | Carbon fiber being one of |
|           | default                   | the main components of    |
|           |                           | Type IV storage tanks,    |
|           |                           | new inventories for       |
|           |                           | carbon fiber              |
|           |                           | manufacturing have been   |
|           |                           | integrated to             |
|           |                           | *carculator*, from [53]_. |
+-----------+---------------------------+---------------------------+
| [54]_     | Type IV hydrogen tank,    |                           |
|           | LDPE liner                |                           |
+-----------+---------------------------+---------------------------+
|           | Type IV hydrogen tank,    |                           |
|           | aluminium liner           |                           |
+-----------+---------------------------+---------------------------+

Table 26 List of inventories for different energy storage solutions

Life cycle impact assessment
****************************

To build the inventory of every vehicle, *carculator* populates a
three-dimensional array *A* (i.e., a tensor) such as:

.. math:: \ A = \left\lbrack a_{\text{ijk}} \right\rbrack,\ i = 1,\ \ldots,\ L,\ j = 1,\ \ldots,\ M,\ k = 1,\ \ldots,\ N

The second and third dimensions (i.e., *M* and *N*) have the same
length. They correspond to product and natural flow exchanges between
supplying activities (i.e., *M*) and receiving activities (i.e., *N*).
The first dimension (i.e., *L*) stores model iterations. Its length
depends on whether the analysis is static or if an uncertainty analysis
is performed (e.g., Monte Carlo).

Given a final demand vector *f* (e.g., 1 kilometer driven with a
specific vehicle, represented by a vector filled with zeroes and the
value 1 at the position corresponding to the index *j* of the driving
activity in dimension M) of length equal to that of the second dimension
of *A* (i.e., *M*), *carculator* calculates the scaling factor *s* so
that:

.. math:: s = A^{- 1}f

Finally, the scaling factor *s* is multiplied with a characterization
matrix *B.* This matrix contains midpoint characterization factors for a
number of impact assessment methods (as rows) for every activity in *A*
(as columns).

As described earlier, the tool chooses between several
characterization matrices *B*, which contain pre-calculated values for
activities for a given year, depending on the year of production of the
vehicle as well as the REMIND climate scenario considered (i.e.,
"SSP2-Baseline", "SSP2-PkBudg1150" or "SSP2-PkBudg500"). Midpoint and
endpoint (i.e., human health, ecosystem impacts and resources use)
indicators include those of the ReCiPe 2008 v.1.13 impact assessment
method, as well as those of ILCD 2018. Additionally, it is possible to
export the inventories in a format compatible with the LCA framework
Brightway2 [51] or Simapro [52], thereby allowing the characterization
of the results against a larger number of impact assessment methods.


References
==========

.. [1] Thiel C, Schmidt J, Van Zyl A, Schmid E. Cost and well-to-wheel implications of the vehicle fleet CO2 emission regulation in the European Union. Transp Res Part A Policy Pract 2014;63:25-42. https://doi.org/10.1016/j.tra.2014.02.018.

.. [2] Ducker Frontier. Aluminum Content in European Passenger Cars. Eur Alum 2019:13.

.. [3] European Commission. Monitoring of CO2 emissions from passenger cars - Regulation (EU) 2019/631 — European Environment Agency 2021. https://www.eea.europa.eu/data-and-maps/data/co2-cars-emission-19

.. [4] European Environment Agency. Monitoring of CO2 emissions from passenger cars - Regulation (EC) No 443/2009 - European Environment Agency 2019. https://www.eea.europa.eu/data-and-maps/data/co2-cars-emission-16

.. [5] Cox B, Bauer C, Mendoza Beltran A, van Vuuren DP, Mutel C. Life cycle environmental and cost comparison of current and future passenger cars under different energy scenarios. Appl Energy2 2020.

.. [6] Bauer C, Cox B, Heck T, Hirschberg S, Hofer J, Schenler W, et al. Opportunities and challenges for electric mobility: an interdisciplinary assessment of passenger vehicles Final report of the THELMA project in co-operation with the Swiss Competence Center for Energy Research "Efficient technologies and systems for mobility." 2016.

.. [7] Sacchi R, Bauer C, Cox B, When, Where and How can battery-electric vehicles help reduce greenhouse gas emissions? Renew Sustain Energy Rev 2022.

.. [8] Veronika Henze. China Dominates the Lithium-ion Battery Supply Chain, but Europe is on the Rise. BloombergNEF 2020. https://about.bnef.com/blog/china-dominates-the-lithium-ion-battery-supply-chain-but-europe-is-on-the-rise/

.. [9] Xinhua. China's CATL unveils cell-to-pack battery platform. 2019. http://www.xinhuanet.com/english/2019-09/13/c_138389934.htm (accessed November 14, 2021).

.. [10] Mark K. BYD's New Blade Battery Set to Redefine EV safety Standards. INSIDEEVs 2020:1-2.

.. [11] BatteryUniversity. BU-216: Summary Table of Lithium-based Batteries - Battery University 2021. https://batteryuniversity.com/article/bu-216-summary-table-of-lithium-based-batteries

.. [12] Yang X-G, Liu T, Wang C-Y. Thermally modulated lithium iron phosphate batteries for mass-market electric vehicles. Nat Energy 2021 62 2021;6:176-85. https://doi.org/10.1038/s41560-020-00757-7.

.. [13] Göhlich D, Fay TA, Jefferies D, Lauth E, Kunith A, Zhang X. Design of urban electric bus systems. Des Sci 2018;4. https://doi.org/10.1017/dsj.2018.10.

.. [14] Preger Y, Barkholtz HM, Fresquez A, Campbell DL, Juba BW, Romàn-Kustas J, et al. Degradation of Commercial Lithium-Ion Cells as a Function of Chemistry and Cycling Conditions. J Electrochem Soc 2020;167:120532. https://doi.org/10.1149/1945-7111/abae37.

.. [15] Cox B, Althaus H-J, Christian Bauer I, Sacchi R, Mutel C, Mireille Faist Emmenegger P, et al. Umweltauswirkungen von Fahrzeugen im urbanen Kontext Schlussbericht. 2020.

.. [16] Schwertner M, Weidmann U. Comparison of well-to-wheel efficiencies for different drivetrain configurations of transit buses. Transp Res Rec 2016;2539:55-64. https://doi.org/10.3141/2539-07.

.. [17] Simons A, Bauer C. A life-cycle perspective on automotive fuel cells. Appl Energy 2015;157:884-96. https://doi.org/10.1016/j.apenergy.2015.02.049.

.. [18] Eudy L, Post M. Fuel Cell Buses in U.S. Transit Fleets: Current Status 2020. 2020.

.. [19] Kurtz J, Sprik S, Saur G, Onorato S. Fuel Cell Electric Vehicle Durability and Fuel Cell Performance. 2018.

.. [20] Mock P. Footprint versus mass: How to best account for weight reduction in the european vehicle CO2 regulation. vol. 2020. 2017.

.. [21] Sebastian BM, Thimons MA. Life Cycle Greenhouse Gas and Energy Study of Automotive Lightweighting. 2017.

.. [22] Hottle T, Caffrey C, McDonald J, Dodder R. Critical factors affecting life cycle assessments of material choice for vehicle mass reduction. Transp Res Part D, Transp Environ 2017;56:241. https://doi.org/10.1016/J.TRD.2017.08.010.

.. [23] World Steel Association. STEEL IN THE CIRCULAR ECONOMY A life cycle perspective. Worldsteel Asscociation 2015:16.

.. [24] Plötz P, Moll C, Bieker G, Mock P, Li Y. REAL-WORLD USAGE OF PLUG-IN HYBRID ELECTRIC VEHICLES FUEL CONSUMPTION, ELECTRIC DRIVING, AND CO2 EMISSIONS. 2020.

.. [25] Spielmann M, Dones R, Bauer C. Life Cycle Inventories of Transport Services. Final report. Dübendorf and Villigen, CH: 2007.

.. [26] Notter B, Keller M, Cox B. Handbook emission factors for road transport 4.1 Quick reference. 2019.

.. [27] Swiss Federal Office for the Environment. Switzerland's National Inventory Report 2021. 2021.

.. [28] Miller J, Jin L. Global Progress Toward Soot-Free Diesel Vehicles in 2019. Icct, Ccac 2019:35. https://theicct.org/sites/default/files/publications/Global_progress_sootfree_diesel_2019_20190920.pdf

.. [29] International Energy Agency (IEA). Extended world energy balances 2021. https://doi.org/https://doi.org/10.1787/data-00513-en.

.. [30] INFRAS. Handbook Emission Factors for Road Transport 2019. https://www.hbefa.net/e/index.html.

.. [31] European Environment Agency. Air pollutant emission inventory guidebook 2019 2019. https://www.eea.europa.eu/publications/emep-eea-guidebook-2019/part-b-sectoral-guidance-chapters/1-energy/1-a-combustion/1-a-3-b-i/view.

.. [32] Beddows DCS, Harrison RM. PM10 and PM2.5 emission factors for non-exhaust particles from road vehicles: Dependence upon vehicle mass and implications for battery electric vehicles. Atmos Environ 2021;244:117886. https://doi.org/10.1016/J.ATMOSENV.2020.117886.

.. [33] US EPA. Emission Factor Documentation for AP-42, Section 13.2.1: Paved Roads. Measurement Policy Group, Office of Air Quality Planning and Standards. 2011.

.. [34] Stolz P, Messmer A, Frischknecht R. Life Cycle Inventories of Road and Non-Road Transport Services. 2016.

.. [35] Stylianos Kephalopoulos, Marco Paviotti FA-L. Common Noise Assessment Methods in Europe (CNOSSOS-EU). vol. 122. 2012.

.. [36] Pallas MA, Bérengier M, Chatagnon R, Czuka M, Conter M, Muirhead M. Towards a model for electric vehicle noise emission in the European prediction method CNOSSOS-EU. Appl Acoust 2016;113:89-101. https://doi.org/10.1016/j.apacoust.2016.06.012.

.. [37] Cucurachi S, Schiess S, Froemelt A, Hellweg S. Noise footprint from personal land-based mobility. J Ind Ecol 2019;23:1028-38. https://doi.org/10.1111/jiec.12837.

.. [38] Entso-e. ENTSO-E Ten-Year Network Development Plan 2020 - Main Report - November 2020 - Version for public consultation 2020.

.. [39] Sacchi R, Terlouw T, Siala K, Dirnaichner A, Bauer C, Cox B, et al. PRospective EnvironMental Impact asSEment (premise): a streamlined approach to producing databases for prospective Life Cycle Assessment using Integrated Assessment Models. Renew Sustain Energy Rev 2022.

.. [40] Cozzolini F. Life Cycle Assessment of Biofuels in EU/CH. 2018.

.. [41] Hank C, Lazar L, Mantei F, Ouda M, White RJ, Smolinka T, et al. Comparative well-to-wheel life cycle assessment of OME3-5 synfuel production via the power-to-liquid pathway. Sustain Energy Fuels 2019;3:3219-33. https://doi.org/10.1039/c9se00658c.

.. [42] Van Der Giesen C, Kleijn R, Kramer GJ. Energy and climate impacts of producing synthetic hydrocarbon fuels from CO2. Environ Sci Technol 2014;48:7111-21. https://doi.org/10.1021/es500191g.

.. [43] Zhang X, Witte J, Schildhauer T, Bauer C. Life cycle assessment of power-to-gas with biogas as the carbon source. Sustain Energy Fuels 2020. https://doi.org/10.1039/c9se00986h.

.. [44] Antonini C, Treyer K, Moioli E, Bauer C, Mazzotti M. Hydrogen from wood gasification with CCS - a techno-environmental analysis of production and use as transport fuel. ChemRxiv 2020. https://doi.org/10.26434/chemrxiv.13213553.v1.

.. [45] Antonini C, Treyer K, Streb A, van der Spek M, Bauer C, Mazzotti M. Hydrogen production from natural gas and biomethane with carbon capture and storage - A techno-environmental analysis. Sustain Energy Fuels 2020;4:2967-86. https://doi.org/10.1039/d0se00222d.

.. [46] Bauer C, Desai H, Heck T, Schneider S, Terlouw T, Treyer K, et al. Electricity storage and hydrogen: Technologies, costs and environmental burdens. 2021.

.. [47] Zhang X, Bauer C, Mutel CL, Volkart K. Life Cycle Assessment of Power-to-Gas: Approaches, system variations and their environmental implications. Appl Energy 2017;190:326-38. https://doi.org/10.1016/j.apenergy.2016.12.098.

.. [48] Dai Q, Kelly JC, Gaines L, Wang M. Life Cycle Analysis of Lithium-Ion Batteries for Automotive Applications. Batter 2019, Vol 5, Page 48 2019;5:48. https://doi.org/10.3390/BATTERIES5020048.

.. [49] Wernet G, Bauer C, Steubing B, Reinhard J, Moreno-Ruiz E, Weidema B. The ecoinvent database version 3 (part I): overview and methodology. Int J Life Cycle Assess 2016;21:1218-30.

.. [50] Supply Chain Looms as Serious Threat to Batteries' Green Reputation \| Greentech Media n.d. https://www.greentechmedia.com/articles/read/graphite-the-biggest-threat-to-batteries-green-reputation

.. [51] Engels P, Cerdas F, Dettmer T, Frey C, Hentschel J, Herrmann C, et al. Life cycle assessment of natural graphite production for lithium-ion battery anodes based on industrial primary data. J Clean Prod 2022;336:130474. https://doi.org/10.1016/J.JCLEPRO.2022.130474.

.. [52] Simons A, Bauer C. A life-cycle perspective on automotive fuel cells. Appl Energy 2015;157:884-96. https://doi.org/10.1016/J.APENERGY.2015.02.049.

.. [53] Benitez A, Wulf C, de Palmenaer A, Lengersdorf M, Röding T, Grube T, et al. Ecological assessment of fuel cell electric vehicles with special focus on type IV carbon fiber hydrogen tank. J Clean Prod 2021;278:123277. https://doi.org/10.1016/j.jclepro.2020.123277.

.. [54] Evangelisti S, Tagliaferri C, Brett DJL, Lettieri P. Life cycle assessment of a polymer electrolyte membrane fuel cell system for passenger vehicles. J Clean Prod 2017;142:4339-55. https://doi.org/10.1016/j.jclepro.2016.11.159.



