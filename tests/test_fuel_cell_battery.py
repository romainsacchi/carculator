"""Matched reference packs, explicit capacities and unaffected powertrains."""

from copy import deepcopy

import numpy as np
import pytest
import xarray as xr

from carculator import CarInputParameters, CarModel, fill_xarray_from_input_parameters


def inputs(powertrains=("FCEV",), years=(2020,), sizes=("Medium",)):
    parameters = CarInputParameters()
    parameters.static()
    return fill_xarray_from_input_parameters(
        parameters,
        scope={
            "size": list(sizes),
            "powertrain": list(powertrains),
            "year": list(years),
        },
    )[1]


def test_2020_pack_matches_all_three_manufacturer_quantities():
    model = CarModel(inputs(), country="CH")
    model.set_all()
    power = model["battery power"].item()
    scale = power / 31.5
    assert model["electric energy stored"].item() == pytest.approx(
        310.8 * 4 / 1000 * scale
    )
    assert model["energy battery mass"].item() == pytest.approx(44.6 * scale)
    assert model["energy battery mass"].item() == pytest.approx(
        model["battery cell mass"].item() + model["battery BoP mass"].item()
    )
    # Captured independent offline prototype with coupled mass/power sizing.
    assert model["electric energy stored"].item() == pytest.approx(1.597572, rel=1e-5)
    assert model["energy battery mass"].item() == pytest.approx(57.31316, rel=1e-5)
    assert model["TtW energy"].item() / 1200 == pytest.approx(0.893862, rel=1e-5)
    # No fabricated charging envelope is activated by the pack calibration.
    assert np.max(model.ecm.regenerative_curtailed_dc_kw) == 0


def test_reference_boundary_retains_historical_inputs_and_future_pack_ratios():
    model = CarModel(inputs(years=(2010, 2020, 2050)))
    model.set_all()
    for year in (2020, 2050):
        capacity = model["electric energy stored"].sel(year=year).item()
        mass = model["energy battery mass"].sel(year=year).item()
        power = model["battery power"].sel(year=year).item()
        assert capacity / mass == pytest.approx(1.2432 / 44.6)
        assert power / mass == pytest.approx(31.5 / 44.6)
    assert model["battery cell energy density"].sel(year=2010).item() == pytest.approx(
        0.1
    )


def test_capacity_override_retains_mass_and_checks_power():
    override = {"capacity": {("FCEV", "Medium", 2020): 2.0}}
    original = deepcopy(override)
    model = CarModel(inputs(), energy_storage=override)
    model.set_all()
    assert override == original
    assert model["electric energy stored"].item() == pytest.approx(2)
    assert model["energy battery mass"].item() == pytest.approx(2 * 44.6 / 1.2432)
    too_small = CarModel(
        inputs(), energy_storage={"capacity": {("FCEV", "Medium", 2020): 0.1}}
    )
    with pytest.raises(ValueError, match="10-second buffer"):
        too_small.set_all()


def test_pack_override_and_repeat_runs_preserve_caller_inputs():
    array = inputs(powertrains=("FCEV", "BEV", "ICEV-p"))
    original = array.copy(deep=True)
    profile = {"nominal_kwh": 2.0, "pack_kg": 60.0, "discharge_10s_kw": 40.0}
    copy = deepcopy(profile)
    baseline = CarModel(array)
    baseline.set_all()
    model = CarModel(array, fuel_cell_battery=profile)
    model.set_all()
    xr.testing.assert_identical(array, original)
    assert profile == copy
    for field in ("TtW energy", "curb mass", "electric energy stored", "battery power"):
        xr.testing.assert_equal(
            model[field].sel(powertrain=["BEV", "ICEV-p"]),
            baseline[field].sel(powertrain=["BEV", "ICEV-p"]),
        )
    cell = model.array.sel(powertrain="FCEV")
    assert (
        cell.sel(parameter="energy battery mass").item()
        / cell.sel(parameter="electric energy stored").item()
    ) == pytest.approx(30)
    previous = model.array.sel(parameter=["TtW energy", "energy battery mass"]).copy()
    model.set_all()
    xr.testing.assert_allclose(
        model.array.sel(parameter=["TtW energy", "energy battery mass"]), previous
    )


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_pack_reference_rejected(value):
    with pytest.raises(ValueError, match="finite and positive"):
        CarModel(inputs(), fuel_cell_battery={"pack_kg": value})


def test_reference_pack_allocation_preserves_labelled_samples():
    array = inputs().isel(value=[0, 0]).assign_coords(value=["light", "heavy"])
    array.loc[dict(parameter="glider base mass")] *= xr.DataArray(
        [0.8, 1.2], dims="value", coords={"value": array.value}
    )
    model = CarModel(
        array.transpose("value", "year", "size", "parameter", "powertrain")
    )
    model.set_all()
    assert list(model.array.value.values) == ["light", "heavy"]
    np.testing.assert_allclose(
        model["energy battery mass"] / model["electric energy stored"],
        44.6 / 1.2432,
        rtol=2e-6,
    )
    assert (
        model["energy battery mass"].sel(value="heavy").item()
        > model["energy battery mass"].sel(value="light").item()
    )


def test_completed_fcev_retains_complementary_power_and_complete_pack():
    inputs = CarInputParameters()
    inputs.static()
    _, array = fill_xarray_from_input_parameters(
        inputs, scope={"size": ["Medium"], "powertrain": ["FCEV"], "year": [2020, 2050]}
    )
    model = CarModel(array, country="CH")
    model.set_all()
    net_stack = model["fuel cell power"] / model["fuel cell own consumption"]
    auxiliary_fields = ["auxiliary energy"]
    # Recover delivered auxiliary DC demand from the independently calculated
    # hydrogen-input trace; it bypasses motor conversion losses.
    auxiliary = (
        model.energy.sel(parameter=auxiliary_fields).sum("parameter").max("second")
        * model["fuel cell system efficiency"]
    ).transpose(*model["power"].dims)
    required_dc = model["power"] / model["electric motor efficiency"] + auxiliary
    np.testing.assert_allclose(
        net_stack + model["battery power"], required_dc, rtol=2e-6
    )
    np.testing.assert_allclose(
        net_stack / required_dc, model["fuel cell power share"], rtol=2e-6
    )
    assert np.all(model["battery power"] > 0)
    assert np.all(model["power battery cost"] > 0)
    np.testing.assert_allclose(
        model["battery cell mass"] + model["battery BoP mass"],
        model["energy battery mass"],
        rtol=2e-6,
    )
    np.testing.assert_allclose(
        model["battery cell mass"] / model["energy battery mass"],
        model["battery cell mass share"],
        rtol=2e-6,
    )
    np.testing.assert_allclose(
        model["electric energy stored"],
        model["battery cell mass"] * model["battery cell energy density"],
        rtol=2e-6,
    )
    # The cell power rating supplies the complement in the power-sized car.
    np.testing.assert_allclose(
        net_stack + model["battery cell mass"] * model["battery cell power density"],
        required_dc,
        rtol=2e-6,
    )
