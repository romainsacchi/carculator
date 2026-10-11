# Changelog

Notable user-facing changes to `carculator`.

## [1.9.7] - 2026-10-11

### Migration to ecoinvent 3.12

- Use IEA PVPS 2026 photovoltaic production mixes through the shared background: country-specific coefficients and generation-weighted GLO/EU27 aggregates replace the fixed Swiss solar proxy. An explicit legacy option remains; exports require matching premise PV suppliers.
- Correct HDPE- and aluminium-lined 700-bar hydrogen tank inventories to one kg of dry assembly and convert manufacturing electricity from MJ to kWh. FCEV tank sizing is unchanged; recalculate affected life-cycle impacts.
- Use the coordinated `carculator_utils 1.3.7` background with premise 2.5.4 and ecoinvent 3.12 cutoff. The foreground matrix, activity index and all 57 LCIA coefficient matrices were rebuilt together. Recalculate saved LCA results; results from the previous background are not directly interchangeable.
- Rename explicit battery chemistry selections and chemistry-specific custom parameter names from `NMC-523` to `NMC-532`. The new inventory uses the actual 5:3:2 Ni:Mn:Co chemistry. Existing engineering and price assumptions carry forward; this is not an alias for the former inventory.
- Inventory exports now default to ecoinvent 3.12. Match exported suppliers against the documented compatible background. Ecoinvent 3.9/3.10 remain explicit legacy targets, but suppliers without verified older links raise an error; requesting an older version does not convert the new background automatically.
- Brightpath supplies Brightway Excel, SimaPro CSV and foreground openLCA JSON-LD exports. Select one sample before export; all selected years are retained. Brightway requires background matching, and openLCA requires provider and elementary-flow mapping before calculation. SimaPro uses Latin-1 and retains unsupported noise quantities in audit comments.
- Python support remains `>=3.12,<3.13` with NumPy `>=1.26.4,<2`. Install matching family versions and rebuild saved models/arrays when applying the corrected inputs, costs or inventory accounting described below.

### Other changes

- Review the documentation against current code, explain units and assumptions
  in plain English, and add validation bar charts with sources and explicit
  limits. Correct outdated descriptions of exports, repeated runs and scientific
  accounting; preserve dated historical results as such.
- Resolve overlapping bundled parameter records into disjoint scopes while preserving every effective static value and uncertainty distribution. Archive original records in `data/overlap_resolution_provenance.json`; regenerate old seeded arrays after this record cleanup.
- Remove import-time global warning suppression and refresh API/documentation contracts; shared sulfur-table reductions use explicit pandas axis arguments.
- Default additional CNG pump-to-tank leakage to zero beyond the delivered-fuel supplier boundary. Remove the unqualified 0.4% overlay while retaining upstream and exhaust methane; explicit measured residuals remain supported.

### Compatibility and installation

- Require the stable `carculator_utils>=1.3.7` release, including its Brightpath runtime dependency.
- Keep core model/LCIA calculations independent of Brightway projects and imports. Export writers now come through Brightpath; `excel` remains a compatibility alias and `brightway` selects the legacy stack (`bw2io<0.9`, `bw2data<4`, `bw2calc<2`).

### Model and inventory changes

- Replace the NMC-523 option with NMC-532 for the ecoinvent 3.12 background refresh. Update explicit chemistry selections and custom parameter names; engineering priors retain their previous numerical values.
- Use the matching shared premise 2.5.4 / ecoinvent 3.12 coefficient bundle; exports default to 3.12 and reject unverified older supplier links.
- Make repeated `set_all()` calls rebuild from retained inputs, with stable costs/energy and retained PHEV components. Preserve explicit input edits and selected sample prices; see the shared repeat-run guide.
- Bill BEVs and PHEV electric operation from grid electricity consumption, including charger losses. Preserve fuel-mode costs, PHEV driving-share weighting and model-specific cost units; verify costs against completed inventory purchases. See [charging cost accounting](docs/validity.rst#charging-cost-accounting).
- Cars now split additional methane leakage between fossil and non-fossil flows according to the fuel blend, instead of classifying all leakage as fossil. Use the shared mass balance, include both origins in impacts/exports, and document the historical loss-rate boundary; see [validation](docs/validity.rst#additional-methane-leakage).
- Make projected costs reproducible with `stochastic(n, seed=...)`, retaining factors across sample/year selections and serialization without using NumPy's global RNG. Keep deterministic static/sensitivity factors and explicit battery prices. `stochastic(1)` now also samples cost factors; regenerate old stochastic cost results.
- Correct year/sample alignment in automatic component-cost projections. Multi-year sensitivity references now match static prices and costs; sampled factors remain attached to their samples across years. Preserve existing price curves, explicit battery prices and physical/inventory outputs.
- Preserve explicit generic and selected-chemistry battery prices through cost adjustment, including scoped zero and per-sample constructor inputs. Verify completed purchase and replacement costs and unchanged default pricing; see [usage](docs/usage.rst#battery-unit-costs).

### Inventory export

- Inherit Brightpath (`>=1.0.0a6,<1.1`, v1 alpha API) writers for Brightway Excel, SimaPro CSV and foreground-only openLCA JSON-LD from `carculator_utils`; document examples and return values in the [export guide](docs/inventory_export.rst).
- Default to ecoinvent 3.12 cut-off; older 3.9/3.10 targets reject suppliers without verified links. Retain one selected sample per export, every selected year and unchanged source inventories/impacts. Brightway importers still require background matching and writing.
- Document the SimaPro Latin-1 layout and identifier changes, warnings for omitted custom noise flows, and openLCA provider/elementary-flow mapping required before calculation. Remove obsolete presamples and uncertainty-export claims.

### Known limitations

- The 2025 inputs combine engineering priors and selected calibration evidence; they are not a measured validation of every vehicle configuration.
- ADAC cycle reconstruction is approximate. Remaining petrol-car residuals must be interpreted with road load, temperature and operating controls.
- Target range supersedes capacity, and capacity supersedes input pack mass. Fixed total mass or explicit consumption may suppress chemistry-dependent consumption changes.

See [release preparation](RELEASING.md) for coordinated versions and verification.

## [1.9.6] - 2026-10-08

### Compatibility and installation

- Require Python 3.12 (`>=3.12,<3.13`); older Python environments must be recreated.
- Use NumPy `>=1.26.4,<2` through the shared runtime.
- Require the stable `carculator_utils>=1.3.6` release, including its export extras.
- Build wheels and source distributions from centralized `pyproject.toml` metadata.
- Keep core model/LCIA use independent of Brightway; install `excel` or `brightway` extras for export. The Brightway extra targets the legacy stack (`bw2io<0.9`, `bw2data<4`, `bw2calc<2`).
- Align documentation versions with the package version and provide complete documentation-build dependencies.

### Model and inventory changes

- Converge BEV target range with battery mass, driving mass, motor power and energy demand, allowing chemistry to affect the required capacity and consumption.
- Verify capacity and pack-mass changes through completed models and inventories across four chemistries, two years and two load samples.
- Add native 2025 inputs and explicit component-efficiency priors, with temporal extensions that avoid artificial 2020/2025/2030 discontinuities.
- Expose opt-in petrol start-stop and fuel-cut energy controls; retain explicit energy and component-efficiency overrides.
- Apply corrected shared battery and hybrid energy accounting and consistently mask unavailable configurations across reported energy and supply outputs.
- Use bounded per-cell sizing, preserve custom inputs, and verify labelled samples, mass balance and annualized costs.
- Verify fuel blends, PHEV carbon accounting, hot pollutant translation and multi-year export through the updated shared dependency.

### Documentation and verification

- Add current installation and executable 2025 quick-start examples, migration notes and a release checklist.
- Record calibration scope, measurement boundaries and numerical consistency separately from empirical validation.
- Verify built wheels and sdist-built wheels, packaged resource hashes, installed tests with export extras and offline core-only model/LCIA smoke runs.

### Known limitations

- The 2025 inputs combine engineering priors and selected calibration evidence; they are not a measured validation of every vehicle configuration.
- ADAC cycle reconstruction is approximate. Remaining petrol-car residuals must be interpreted with road load, temperature and operating controls.
- Target range supersedes capacity, and capacity supersedes input pack mass. Fixed total mass or explicit consumption may suppress chemistry-dependent consumption changes.

See [validation](docs/validity.rst) and [release preparation](RELEASING.md) for scope and verification instructions.

## [1.9.5] - 2026-04-29

### Added

- Added `AGENTS.md` with repository-specific guidance for automated coding agents.
- Added this changelog.

### Changed

- Replaced the legacy `setup.py` packaging configuration with `pyproject.toml`.
- Declared supported Python versions and package data in `pyproject.toml`.
- Updated documentation release metadata to `1.9.5`.
- Made inventory tests use fresh model instances and temporary export directories.
- Converted the model reference-output test from a fixture-writing smoke test into assertions.

### Fixed

- Fixed `CarInputParameters` so caller-provided parameter dictionaries and file paths are honored.
- Fixed FCEV cost adjustment so fuel-cell stack costs update `fuel cell cost per kW` instead of overwriting hydrogen tank costs.
- Fixed vehicle purchase-cost amortisation to use lifetime in years, with zero-lifetime safeguards.
- Fixed component replacement discounting to use year-based lifetime discounting.
- Fixed CNG pump-to-tank methane leakage so direct methane emissions are added to the inventory.
- Fixed invalid impact indicator handling so typos raise `ValueError` before xarray result construction.
- Made export tests skip cleanly when the optional `bw2io` export stack is not importable in the active Python environment.
- Fixed the conda publication workflow to upload the package format produced by `conda-build`.
- Added `pip` and runtime `python` requirements to the conda recipe to support conda build publication.
- Aligned `carculator_utils>=1.3.5` across packaging metadata.
