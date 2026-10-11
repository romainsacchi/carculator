# Releasing carculator 1.9.7

Status: prepared, not published. Use Python 3.12 and review the matching family set:

| Package | Prepared version |
| --- | --- |
| `carculator_utils` | 1.3.7 |
| `carculator` | 1.9.7 |
| `carculator_truck` | 0.5.2 |
| `carculator_bus` | 0.1.2 |
| `carculator_two_wheeler` | 0.1.1 |


1. Review `CHANGELOG.md`, known limitations in `docs/validity.rst`, and the worktree.
   Preserve local notebooks and inventories; commit only the intended release files.
2. Check `carculator/_version.py`, `conda/meta.yaml` and dependency metadata together.
   Sphinx reads the package version automatically. All vehicle packages require
   `carculator_utils>=1.3.7`, whose runtime includes Brightpath export writers.
   The `excel` extra is a compatibility alias; `brightway` selects the legacy stack.
   Review [export formats and linking requirements](docs/inventory_export.rst).
3. Build and verify actual artifacts from this checkout (use a new output directory):

```bash
python ../carculator_utils/scripts/verify_installation.py --repositories ../carculator_utils . --output /tmp/carculator-release-check --run-tests
```

For the complete family, run from `carculator_utils`:

```bash
python scripts/verify_installation.py \
  --repositories . ../carculator ../carculator_truck ../carculator_bus ../carculator_two_wheeler \
  --output /tmp/carculator-family-release-check --run-tests
```

The verifier builds wheels and sdists, rebuilds wheels from sdists, compares
packaged resource hashes, runs installed tests with export extras and checks
offline core-only model/LCIA runs. Keep `report.json`, test XML, dependency freezes
and artifact hashes with the release record. A local pass does not establish
that Linux, Windows or hosted CI passed. Confirm the release commit's CI separately.

4. Install `docs/docs_requirements.txt`, build with `python -m sphinx -b html docs docs/_build/html`,
   and run the README example. Check wheel/sdist metadata before publication:

```bash
python -m pip install twine
python -m twine check --strict /tmp/carculator-family-release-check/wheels/carculator-1.9.7-py3-none-any.whl /tmp/carculator-family-release-check/wheels/carculator-1.9.7.tar.gz
```

5. Publish `carculator_utils` **first** and confirm that 1.3.7 is available before
   publishing the four vehicle packages. The same ordering applies to conda.
   Replace `Unreleased` with the actual publication date in the changelog, commit
   the final release metadata, and verify the final artifacts. Create an annotated
   `v.1.9.7` tag on that reviewed commit, following existing tag conventions.
   Do not reuse an existing tag or upload a different build under the same version.
6. Publish the GitHub release for the reviewed tag. The workflow verifies
   installed artifacts on Linux, macOS and Windows before publishing the
   verified wheel and sdist to PyPI. Attach the verification record.
7. The workflow separately builds/tests the conda package and uploads it
   to the `romainsacchi` channel. Ensure utils is available there first.

An empty `release_tag` in a manual workflow run only verifies. A nonempty
tag also publishes. Creating a tag alone does not publish.

The prepared metadata and README examples target this release; older published
packages may not provide the documented APIs or 2025 defaults.


## Preparation on 2026-10-11

Prepare the versions above together; publish utils before its consumers.
The previous utils/car/truck/bus versions were published on 2026-10-08
and must not be reused. Two-wheeler 0.1.1 remains unpublished.

The [candidate record](docs/_static/release_candidate_20261011.json) records
the source revisions, green GitHub Actions runs and fresh local artifact
verification. Historical reports below describe their original builds.
The solar-PV and dry-mass hydrogen-tank updates in the original utils
working tree are excluded from this release candidate. Preserve them
for a separate reviewed and CI-tested change.

After committing this preparation, verify all five final commits in CI.
Downstream verification must use the prepared utils revision via
`CARCULATOR_UTILS_REF`; the utils family job must use the matching vehicle
revisions via its `CARCULATOR_*_REF` variables. Record the resolved SHAs.
Restore temporary CI references after the coordinated release.
Conda build/publish jobs were skipped by the ordinary verification runs;
their success must be checked during publication.

## Historical verification record

On 2026-10-08, the five-package installed-artifact suites passed **497 tests**,
with one existing expected two-wheeler cost failure. Wheel and sdist-built
wheel resource checks, offline core-only model/LCIA runs, strict Twine metadata
checks, README execution and the documented inventory exports passed. All five
Sphinx sites built using the release wheels and their ``docs`` extras.

The [release verification record](docs/_static/release_verification.json)
contains versions, artifact hashes, test counts and qualifications. Builds were
local on macOS with Python 3.12; hosted CI and conda builds need separate
qualification. Existing documentation warnings are recorded. These checks
exercise packaging and software consistency; they do not establish physical
plausibility or replace the measurement evidence and limitations in [model validation](docs/validity.rst).
