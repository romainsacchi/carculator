.. _validity:

Passenger-car calibration and validation
========================================

Start with :doc:`validation_examples` for bar charts and an explanation of the
evidence. This page records detailed checks and limitations; dated test totals
and before/after results refer to their stated software snapshots.

The 2025 inputs combine historical parameter trends, sourced component assumptions
and explicit engineering assumptions. They are not a newly fitted fleet of
2025 cars. Historical curb-mass calibration is described in :doc:`modeling`;
it does not independently validate current fuel or electricity consumption.

Installed PHEV battery capacity
-------------------------------

Both driving-mode cost calculations price the full electric-mode battery pack.
The combustion-mode intermediate's smaller operating storage no longer reduces
installed battery purchase cost when modes are combined. Purchase-dependent
financing and, for trucks, insurance and residual credit follow that correction.
Unit-price inputs, markup and explicit cost overrides retain their conventions.
The shared model also preserves nominal capacity and computes electric range
from usable capacity. Battery replacement counts are unchanged pending a
separate lifetime-use review.

.. _charging-cost-accounting:

Charging cost accounting
------------------------

Electricity running costs use grid purchases: ``electricity consumption`` in
kWh/km times the electricity tariff. Grid consumption already includes both
battery-charge and charger losses; neither efficiency is applied again when
billing that electricity. Previously the cost formula omitted charger losses.
At 90% charger efficiency it understated the electricity component by 10%; at
80% efficiency it understated it by 20%. This correction changes costs, while
preserving vehicle energy demand, inventory electricity exchanges and LCIA.

BEVs and PHEV electric intermediates use this grid-based calculation. Combined
PHEVs retain the utility-factor-weighted sum of electric and combustion costs,
with the electric share applied once. Fuel-mode costs retain their existing
convention. Tariffs and charging-efficiency assumptions have not been refitted.

Car costs remain per vehicle-km. The default Swiss 2025 ``Medium`` BEV on
WLTC uses 16.95 kWh/100 km from the grid; at the bundled tariff its electricity
cost is approximately EUR 3.40/100 km, corrected from EUR 3.06/100 km.

Completed model/inventory checks and the shared billing calculation are described
in the `shared charging-cost validation <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/validity.rst#charging-cost-accounting>`_.


Target-range sizing consistency
-------------------------------

The BEV range-sizing loop now converges battery pack mass together with driving
mass, power and energy demand. For a 2025 Medium BEV on WLTC with a 400 km target,
LFP requires 78.45 kWh nominal capacity and 17.70 kWh/100 km from the grid;
NMC-811 requires 74.15 kWh and 16.73 kWh/100 km. A frozen-vehicle energy
recalculation at each final mass confirms the requested range within the sizing
tolerance. Default runs without a target retain the recorded baseline outputs.

Regression tests cover four chemistries, two years and two load samples,
independent fixed-capacity runs, user overrides, mixed powertrains and bounded
convergence. Completed inventories verify both charging electricity and battery
production exchanges, with finite life cycle impacts. These establish numerical
and accounting consistency; they are not measured-vehicle validation.
See the `reproduction and repaired results <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/bev_target_range_issue.rst>`_.

Capacity and pack-mass changes were also tested in 24 complete model/inventory
runs: four chemistries, three levels per input direction, two years (2020/2025)
and two load samples, or 96 vehicle/year/sample cells. Increasing nominal capacity
from 40 to 80 kWh or pack mass from 250 to 550 kg increases energy demand and range
throughout this tested scope. Increasing pack mass also increases nominal capacity.
For a 2025 Medium NMC-811 BEV on WLTC, the capacity sweep changes range from
226.1 to 428.2 km and grid demand from 15.97 to 16.86 kWh/100 km. The mass sweep
changes capacity from 52.56 to 115.64 kWh, range from 291.9 to 591.1 km, and grid
demand from 16.25 to 17.66 kWh/100 km. Fixed curb mass, energy consumption or target
range were not imposed in these sweeps.

Mass/energy consistency and battery/electricity inventory exchanges pass the
regressions. A related repair restricts capacity overrides to selected vehicles,
preserving unrelated FCEV battery component masses in mixed runs. The shared
report linked above includes the full sweep results and source provenance.

Evidence reviewed
-----------------

The shared measurement catalog contains 19 paired car observations covering
18 configurations. Those baseline runs use WLTC, whereas the ADAC observations
use Ecotest. Driving masses are measured or reconstructed from documented test
payloads, with the basis retained in the catalog. Matching mass alone does not
match road load, cold starts, auxiliaries, temperature or the driving protocol.

For four mini BEVs, an additional experiment reconstructs ADAC's electric-cycle
ending from its published graph. It is not an official numerical trace. The
following are AC charging electricity, in kWh/100 km:

.. list-table:: Cycle sensitivity, without consumption fitting
   :header-rows: 1

   * - Vehicle
     - WLTC
     - Graph reconstruction
     - ADAC
   * - Dacia Spring 65
     - 13.28
     - 14.70
     - 16.70
   * - Fiat 500e RED
     - 14.19
     - 15.59
     - 15.90
   * - Fiat 500e La Prima
     - 14.88
     - 16.28
     - 17.10
   * - Fiat 500e Cabrio
     - 14.98
     - 16.38
     - 17.40

The graph trace exceeds rated shaft power briefly; a slower-acceleration probe
checks that limitation. These results support investigating cycle mismatch,
not fitting generic motor or battery efficiencies to the remaining residuals.
See `ADAC reconstruction and source protocol <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/adac_cycle_comparison.rst>`_.

For lower-medium petrol cars, the Golf diagnostic gives about 7.03 L/100 km
without hybrid assistance and 6.53 L/100 km with the source drag coefficient
and opt-in start-stop/fuel-cut controls, against the reported ADAC 5.6 L/100 km.
The protocol and control state remain imperfectly matched. This is a diagnostic,
not a class-wide petrol calibration. The controller is disabled by default;
its thresholds, buffer and restart costs are engineering assumptions. See
`petrol diagnostics <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/petrol_car_energy_diagnostics.rst>`_ and
`control accounting and public API <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/combustion_controls.rst>`_.

Petrol full hybrids and depleted petrol PHEVs use the pinned FASTSim 2016
Prius Atkinson engine curve at shaft load. The effective downstream efficiency
of 0.8 remains a class assumption, not a separately measured transmission
parameter. This transfer is not a time-resolved power-split controller.

Hybrid regeneration and battery boundaries have analytical regression tests.
Hybrid motor peak ratings are independent of combined system power; generic
ratios and component assumptions are not independently identified by consumption
alone. Neither these checks nor the historical fleet plots establish universal
validation of hybrids, plug-in hybrids or fuel-cell cars.

Energy boundaries and time trends
---------------------------------

``TtW energy`` is kJ/km. For a BEV it represents net stored-energy depletion;
``model.battery_terminal_energy`` reports net terminal DC energy separately.
``electricity consumption`` is charging electricity in kWh/km. Multiplying it
by 100 gives kWh/100 km. A meter boundary must be identified before comparing
these outputs. Regeneration and battery/charger losses must not be counted twice.

The 2025 motor/inverter (0.90), electric transmission (0.97), charger (0.90)
and symmetric battery one-way (sqrt(0.97)) values are component assumptions in their
documented scopes, not universally measured efficiencies. For relevant hybrid
scopes, the independent motor peak/system-power ratio is 0.65. The temporal
update preserves all 2025 scalar values and uncertainty distributions. Storage
and charger trends preserve relative legacy losses; newly explicit component
assumptions are extended across tabulated years to avoid interpolating from missing
zero values. Historical estimates and future projections therefore change.

The family audit completes 546 annual cases (21 configurations, 2015–2040),
including availability-masked historical cells. The former inputs caused 20
sizing failures in this grid. All 40 existing 2025 measurement-comparison runs
retain their energy use and driving mass exactly. These are consistency and
regression checks, not 546 empirical validations. See
`temporal method, plots and limitations <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/temporal_energy.rst>`_.

Reproducibility
---------------

The installed package includes :download:`2025 record provenance
<../carculator/data/defaults_2025_provenance.json>` and
:download:`temporal provenance and original affected records
<../carculator/data/temporal_energy_provenance.json>`. Overrides should use measured
vehicle-specific inputs where available. Retain source, year, cycle, driving
mass, meter boundary and uncertainty assumptions with each comparison.

With matching Python 3.12 sibling checkouts, run from ``carculator_utils``::

   python scripts/validate_energy_measurements.py --output /tmp/measurements-new
   python scripts/audit_energy_time_trends.py --output /tmp/temporal-new

The shared `measurement catalog and outputs <https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/energy_measurements.rst>`_
record excluded observations as well as paired values. Multiple cycles of one
vehicle and AC/DC measurements from one run are not independent vehicles.

Plug-in hybrid physical packs
-----------------------------

The 2020 PHEV calibration uses matched manufacturer energy/mass pairs or explicit
size/density proxies, recorded in ``data/phev_2020_calibration.json``. Capacities
range from 8.8 to 24 kWh. Both electric and combustion intermediates carry the
same complete physical pack during sizing; electric driving share does not
reduce installed capacity or its purchase cost. The existing 80% usable-energy
window and NMC-622/cell-allocation assumptions remain modelling proxies.

The approved battery life is 200,000 vehicle-km with the existing 18-year limit.
The broader gross-throughput/cycling review was skipped at campaign closure;
these assumptions do not establish measured battery durability.

Fuel-cell car buffer batteries
------------------------------

From model year 2020, FCEV batteries use a matched complete-pack reference from
Toyota's second-generation Mirai announcement of 30 November 2020: 310.8 V,
4 Ah (1.2432 kWh), 44.6 kg, and 31.5 kW discharge for 10 seconds. See
`Toyota's specification <https://newsroom.toyota.eu/introducing-the-all-new-toyota-mirai/>`_.
Capacity and mass scale together with the required buffer power inside the
existing mass/power iteration. This is a class-level proxy, not a calibration
of each vehicle to a measured FCEV. The default Medium car under the Swiss WLTC
assumptions has about 1.60 kWh and 57.3 kg, using 0.894 kg H2/100 km.

The same pack ratios are retained in future model years; no improvement is
inferred from this single reference. Historical inputs before 2020 are retained.
The existing cell-mass fraction only allocates the complete pack internally,
and the existing chemistry and cost inputs remain proxies. Override the matched
reference through ``CarModel(..., fuel_cell_battery={"nominal_kwh": ...,
"pack_kg": ..., "discharge_10s_kw": ...})``. Explicit capacity overrides retain
that capacity and must supply the required 10-second buffer power.

The discharge rating does not establish continuous power, charging capability,
SOC/temperature limits or lifetime. Charging remains unconstrained, and the
cycle calculation assumes instantaneous fuel-cell response. This change cannot
validate fuel-cell transients or sustained battery assistance. For the Premise
adapter, batteries still use the approved MIX capacity market in nominal kWh;
its ordinary battery mass per kWh differs from this high-power physical pack.
Lower purchased kWh must therefore not be interpreted as a verified reduction
in high-power battery manufacturing burdens.

Additional methane leakage
--------------------------

Cars now split additional methane leakage between fossil and non-fossil flows according to the fuel blend, instead of classifying all leakage as fossil.
The shared calculation preserves the existing convention: loss in kg per km
is engine fuel times the ``CNG pump-to-tank leakage`` ratio, and purchased fuel
is engine fuel plus that loss. Fossil/non-fossil methane follows the blend.
Both generic-air methane flows now enter non-exhaust impacts and exports;
combustion CO2 and HBEFA exhaust emissions are unchanged.

The current default is **zero additional leakage** at all tabulated years.
The former 0.4% assumption combined several delivery and vehicle stages and did
not establish an extra loss after the selected supplier. Losses already present
in the supplier inventory and HBEFA exhaust factors remain included. Set
``CNG pump-to-tank leakage`` only when evidence supports an additional loss
outside that supplier's boundary; a value of zero does not mean the whole gas
supply chain is leak-free.
See the shared `methane leakage boundary and verification notes
<https://github.com/Laboratory-for-Energy-Systems-Analysis/carculator_utils/blob/master/docs/methane_leakage.rst>`_.

Completed regressions include Medium cars, diesel/BEV controls,
2020/2025/2030, two samples, fossil gas, sewage biomethane, biological synthetic
methane and mixed fuels. They check mass balance, unchanged upstream inventories,
ReCiPe/EF climate contributions, and repeated multi-year Brightway exports.
These establish accounting consistency, not measured leakage-rate validation.
