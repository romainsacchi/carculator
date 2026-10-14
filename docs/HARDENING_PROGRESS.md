# Hardening implementation record

> Historical development record from 6 October 2026. Its findings, proposed work,
> test counts and environment support describe that stage of development.
> For current behaviour and evidence, use the [user guide](usage.rst),
> [validation examples](validation_examples.rst) and [changelog](../CHANGELOG.md).


Date: 2026-10-06. This is an unreleased first implementation of the
[robustness roadmap](ROBUSTNESS_PLAN.md), following the
[seven-repository audit](REVIVAL_AUDIT.md). It is not a claim that every roadmap
stage or every scientific scenario has been qualified.

## Implemented changes

### Shared behavior and numerical contracts

- All vehicle input subclasses honor custom parameter dictionaries/files and
  extra labels. Parameters, scopes, model arrays, and model overrides are owned
  copies, preventing one calculation from changing the next one's inputs.
- The array builder's coordinate dictionaries describe actual array positions.
  Sensitivity values now have named `reference` and one-at-a-time 10% perturbation
  coordinates. Static and stochastic numeric value coordinates remain supported.
- `stochastic(n, seed=...)` uses a local stats_arrays RNG. Seeded samples are
  reproducible without changing NumPy's global RNG. Invalid sample counts fail.
- Nested model selection contexts restore their enclosing scope after exceptions.
- Battery physics is applied even when chemistry-specific cost metadata is
  absent. Density, mass share, and cycle life are required; generic costs remain
  the fallback, recorded on `battery_cost_fallbacks`. Bus battery defaults no
  longer replace explicit chemistry, origin, or capacity overrides.
- Sizing checks each vehicle/sample with absolute plus relative tolerances and a
  configurable 100-iteration limit. Errors identify the parameter and coordinates.
  Truck convergence ignores pre-2020 electrified powertrains and their PHEV
  intermediates, consistently with the existing availability policy. The policy
  and scientific dataset have not been broadened.
- Noise logarithms no longer leave masked output memory uninitialized at stops;
  sound-power exponentiation uses float64 and incorporates the scale in the exponent.
- Passenger/tonne-km loads align by size, powertrain, year, and sample. Active
  vehicles require positive finite loads. Export normalizes an independent copy,
  so repeated exports preserve A and activity mappings. Multi-year Brightway
  file/importer output returns every year. Canonical passenger/tonne-kilometre
  labels are preserved in export; SimaPro maps passenger kilometers to personkm.
  Invalid software/format pairs and LCIA method/indicator names fail early.

### Tests and consumers

- Regression tests cover the repaired defects, seeded sampling, labelled failure
  diagnostics, load normalization, repeated export, and missing export extras.
- Bus/truck/two-wheeler model construction happens in fixtures, not collection.
  Mutable model instances are copied per test. Unasserted workbook generation is
  replaced by physical assertions; export files use temporary directories.
  Cars reuse a built model template while retaining an independent copy per test.
- Car sensitivity is tested through a real model and pkm impact calculation.
- The API PHEV-p range override now updates range instead of payload; a constructor
  boundary test checks both intermediate powertrains and preserves payload.
- The online application has pure helper tests. Serialization preserves booleans
  and handles NumPy arrays. Calculation data paths resolve from the application
  location rather than the worker's current directory. These are not substitutes
  for full Flask/Redis/database/browser integration.
- Real car export regressions cover vkm/pkm/tkm, two years, both file formats,
  repeated calls, and unchanged LCA results after export.

### Second batch: focused unit tests and cost arithmetic

The unit-test expansion adds 67 pytest cases across all five core packages:
14 shared financial cases, 11 car cases, 14 bus cases, 14 truck cases, and 14
two-wheeler cases. Generated tests run 30 deterministic mass-balance examples per
vehicle library and 50 discounted-payment examples in utils; these examples are
additional executions within the reported pytest case counts.

- Each vehicle suite checks curb/driving/cargo mass balance, lightweighting,
  passenger loads, and batteries using controlled inputs. Bus luggage and truck
  available payload have their own expectations.
- Reordered input dimensions and two labelled samples must preserve values.
  Independent model instances and their caller's input arrays must remain isolated.
- Car, bus, and two-wheeler cost tests use an independent sum of discounted annual
  payments, with 5/10-year lifetimes, two annual mileages, and 0/5% interest.
  They check purchase, midpoint replacement, charging efficiency, total costs,
  and the bus passenger-km divisor. These tests do not use golden spreadsheets.
- Truck infrastructure tests cover scalar zero and near-zero interest,
  negative-but-valid interest, zero demand, charger energy limits, and labelled
  sample/year broadcasting. Mixed diesel/BEV arrays explicitly mask absent
  chargers before financial validation; active chargers require a positive
  lifetime. The zero-demand surcharge policy remains zero.
- The shared capital-recovery helper is checked against present-value sums,
  including generated rates/lifetimes, domain validation, and zero-rate limits.
  Hypothesis is a test extra in every package, never a core runtime dependency.

The initial regressions reproduced bus/two-wheeler annuities calculated using
lifetime **kilometres** and replacement discounting that collapsed to zero.
Those methods now convert distance to lifetime **years**, use stable annual
capital recovery, and discount replacements at the already documented midpoint
of vehicle life. Zero-lifetime/zero-mileage cells accrue no capital annuity;
operating costs for unsupported combinations still require an availability policy.
Truck infrastructure uses the same helper: a scalar zero interest rate previously
raised `ZeroDivisionError`, and near-zero rates suffered cancellation error.
Passenger-car cost calculations already passed the independent oracle and were
left unchanged.

These tests exercise formulas and ownership; they do not establish scientific
validity for every default dataset or powertrain. Battery replacement policy,
PHEV boundaries, repeated full-model execution, and reviewed numerical reference
fixtures remain follow-up work.

### Packaging and CI

All five core packages use setuptools pyproject metadata, explicit scientific
resource inclusion, and one literal `_version.py::VERSION`. Existing public
`__version__` tuples are preserved. setup.py is a shim where it existed;
requirements.txt installs the project metadata. Conda recipes are aligned but
have not been built or installed in this work.

| Package | Candidate version |
| --- | --- |
| carculator_utils | 1.3.6.dev0 |
| carculator | 1.9.6.dev0 |
| carculator_bus | 0.1.1.dev0 |
| carculator_truck | 0.5.1.dev0 |
| carculator_two_wheeler | 0.1.1.dev0 |

The candidate Python range is 3.11–3.12 with NumPy <2. See the executed
verification section for the locally tested interpreters. Non-local platforms
remain qualification targets until their CI jobs execute successfully.
Core modeling needs neither Brightway nor Excel exporters. Use `[excel]` for
Brightway XLSX output and `[brightway]` for bw2io importers. The latter constrains
bw2io to 0.8.12, bw2data to 3.x and bw2calc to 1.x: a clean install demonstrated
that bw2io 0.8.12 plus bw2data 4.7 passes pip check but crashes on import because
it compares a string version with a tuple. Future Brightway support needs its own
migration and export tests.

Each library's CI builds actual artifacts, installs outside the checkout, checks
bundled-resource hashes, runs model/LCA smoke cases without calculation network
access, and runs tests against the installed wheels with export extras. The
matrix targets Linux, macOS, and Windows with Python 3.11/3.12. Repository
permissions are read-only. The former formatting commits, ordinary-push TestPyPI
and Anaconda publication, and embedded publishing jobs have been removed.
At that stage, no replacements published automatically. Release publishing has
since been restored for published GitHub releases, gated on the artifact
verification matrix; see [RELEASING.md](RELEASING.md) for triggers and recovery.

The harness is maintained in `carculator_utils/scripts/verify_installation.py`.
The car repository's script delegates to its sibling utils checkout. Example:

```sh
python scripts/verify_installation.py \
  --repositories ../carculator_utils . ../carculator_bus ../carculator_truck ../carculator_two_wheeler \
  --output /tmp/carculator-artifacts --run-tests
```

The output directory must be new. Dependency installation needs network access.
The harness uses temporary venvs and a private Brightway data directory; it does
not install into the caller's environment. Resource SHA-256 checks cover all 132
bundled JSON/YAML/CSV/NPZ/XLSX files in both direct wheels and wheels rebuilt from
sdists. Reports include resolved versions, artifact hashes, and numerical outputs.

Vehicle CI checks out utils from `CARCULATOR_UTILS_REF` (repository variable), or
master when unset. During development, point it at the matching candidate commit.
Until the utils changes exist on GitHub, these candidate vehicle workflows cannot
pass against the older master. Release order is utils, vehicles, then consumer
pin upgrades. Record exact coordinated commits before release. The Flask TCS Git
pins and online application's existing carculator pin deliberately remain unchanged.

Action configuration was checked against the official
[checkout](https://github.com/actions/checkout),
[setup-python](https://github.com/actions/setup-python), and
[artifact upload](https://github.com/actions/upload-artifact) documentation.

## Numerical change review

The same Python 3.11 environment ran original and candidate source side by side.
Cases below use 2020, CH, vkm; the truck uses Long haul. Climate totals are
kg CO2-eq/vkm. Complete sampled values are in `verification/model-comparison.json`.

| Case | Original TtW kJ/km | Candidate TtW | Original climate total | Candidate climate total |
| --- | ---: | ---: | ---: | ---: |
| Medium BEV car | 609.3123 | 609.3123 | 0.1376282 | 0.1376282 |
| 13m-city BEV-depot bus | 3594.1450 | 3814.0332 | 0.4215693 | 0.4831641 |
| 40t BEV truck | 5327.0063 | 5327.0063 | 1.3271250 | 1.3271250 |
| Bicycle <25 BEV | 25.0089 | 25.0089 | 0.0286840 | 0.0232903 |

The bus battery cell density changes from the erroneous zero to the existing
NMC-622 input, 0.24 kWh/kg. Its cell mass changes from 234.10 to 1035.09 kg as
sizing starts using the chemistry physics. The bicycle's density changes from
0.20 to the selected chemistry's 0.24; cell mass stays 2.31 kg. These are expected
consequences of fixing the partial-metadata battery setter, not new scientific
assumptions. They still require domain review before releasing changed results.
These four cases are regression evidence, not comprehensive scientific validation.
The only changed data resource is the SimaPro passenger-kilometre unit alias;
vehicle parameter datasets and background matrices are unchanged.

### Cost changes from the second batch

The same Python 3.11 environment compared the first hardening candidate with the
unit-test/cost correction candidate, using CH and 2020. Runtime input datasets
are unchanged. Values are currency units per stated distance unit.

| Case | Cost component | Before | After |
| --- | --- | ---: | ---: |
| 13m-city BEV-depot bus, per passenger-km | Annualized purchase | 0.0184612 | 0.0373004 |
| Same bus | Annualized replacement | 0.0000000 | 0.0079828 |
| Same bus | Total | 0.0648389 | 0.0916609 |
| Bicycle <25 BEV, per vehicle-km | Annualized purchase | -0.0042019 | -0.0108832 |
| Same bicycle | Annualized replacement | 0.0000000 | 0.0062946 |
| Same bicycle | Total | -0.0121237 | -0.0125105 |

The bicycle glider and purchase costs remain negative; the correction does not
validate those input coefficients. The strict expected-failure test remains.
Full input/output values are in [unit-cost-comparison.json](verification/unit-cost-comparison.json).

## Remaining work and known issues

1. **Scientific data/cost review:** the electric bicycle's glider cost is negative
   under current inputs, producing approximately -0.01251 currency units/km after
   correcting the annuity. A strict expected-failure test records this. The
   bus/two-wheeler annuity defects are repaired in the second batch, but replacement
   policies and unsupported two-wheeler operating costs still need review. Do not
   hide these by coercing all costs to zero.
2. **Fuel mapping integrity:** `methane - synthetic - electrochemical` references
   an activity missing from the shipped A-matrix index. It now fails visibly.
   The old car test changed caller fuel types but accidentally reused names
   cached by an earlier model. The regression now verifies rejection. Biological
   synthetic methane also maps to a sewage-sludge fuel; resolve scientific
   equivalence, compression, units, and provenance before changing mappings.
3. **Input/data schemas:** overlapping records still use legacy first-entry
   precedence, and missing static cells still become zero. Build an explicit
   availability mask and schema validation before replacing that behavior.
   Custom cycles, interpolation policy, year ranges, and fuel shares need broader
   property-based tests and contextual validation.
4. **Lifecycle and scientific regression:** characterize repeated `set_all()` and
   all stateful hooks; expand multi-sample, PHEV, zero-rate/zero-distance, cost and
   functional-unit cases. Establish reviewed small golden fixtures across fuels,
   scenarios and vehicle modes. Do not relax convergence to disguise divergence.
5. **Dependency/packaging qualification:** execute the committed OS/Python matrix,
   test minimum dependency bounds and editable installs, build Conda separately,
   check documentation examples, and add tag/version/release provenance checks.
   Pandas/xarray deprecation warnings remain in legacy methods. The new workflows
   have been inspected locally, not executed by GitHub Actions.
6. **Consumer integration:** validate real API responses for all vehicle types and
   both TCS/SwissCargo branches; build local Redis/database/worker tests, result
   expiry/download/browser tests, and runtime version labels. Non-finite-to-zero
   web serialization remains an explicitly documented legacy display policy.
7. **Release preparation:** review result changes and the xfail, pin exact family
   commits for consumers, publish tested artifacts through a separately gated
   release workflow, and exercise rollback. Nothing has been pushed or published.

## Executed verification

Executed on macOS arm64 with Python 3.11.12 and Python 3.12. Artifact builds,
clean installs, pip check, optional export imports, resource hashes, offline
model/LCA smoke calculations, and wheel-versus-sdist numerical comparisons passed
on both interpreters. Full tests ran against the installed wheels outside the
checkout. The four smoke-case outputs also agree across Python versions. They
match the first hardening batch's energy and climate results exactly. Initial truck tests
exposed the inactive-charger boundary; after correction, its wheel and sdist were
rebuilt and its suite rerun. Already passing suites retained their identical
artifacts. Final verification additionally compared all artifact Python source
bytes against the reviewed sources, as well as checking resource hashes.

| Suite | Python 3.11 installed wheels | Python 3.12 installed wheels |
| --- | ---: | ---: |
| utils | 54 passed | 54 passed |
| car | 48 passed | 48 passed |
| bus | 30 passed | 30 passed |
| truck | 34 passed | 34 passed |
| two-wheeler | 18 passed, 1 expected failure | 18 passed, 1 expected failure |
| **Core total** | **184 passed, 1 expected failure** | **184 passed, 1 expected failure** |

The earlier batch passed 8 API source tests and 2 online pure-helper tests on
Python 3.11 with explicit sibling source paths. Those consumer checks were not
rerun for this cost-arithmetic expansion and did not qualify the applications'
pinned production dependencies. The current core suite has 184 passing cases and
one documented strict expected failure. No export integration was skipped for a
missing dependency.

See [test-summary.json](verification/test-summary.json),
[Python 3.11 artifacts](verification/python311-artifacts.json), and
[Python 3.12 artifacts](verification/python312-artifacts.json). Adjacent environment
snapshots record the resolved minimal and export/test dependencies. They are
observations of this platform, not universal dependency locks.

Black and isort checks pass for the 10 Python files changed in this batch. Workflow
definitions have read-only repository permissions and no push/publication steps; GitHub has
not executed the matrix. Linux, Windows, Conda, full web services, and production
consumer pin upgrades remain unqualified. Deprecation warnings in legacy
background/presentation methods are recorded rather than treated as test success
without qualification.

## Source commit provenance

`../compatibility.json` records the committed candidate source revisions alongside
the pre-hardening baseline hashes. The car revision identifies its source/build
commit before this documentation commit. AGENTS.md files remain local and ignored;
previously tracked copies were removed from the index without deleting local files.

## FCEV electrical sizing correction (2026-10-14)

Fuel-cell and buffer ratings now cover motor shaft power divided by the existing
motor/inverter efficiency, plus auxiliary DC demand. The original percentage
split applies at that electrical boundary; stack own consumption is included
once. The shared implementation is in `carculator_utils`; the car sizing loop
already calculates auxiliaries before the fuel-cell components. Native tests
check the supply against the independently calculated energy trace and retain
the complete-pack mass/capacity identities. This does not establish battery
charge/discharge capability or replace the separate empirical pack review.
