"""A calibrated PHEV carries one complete pack in both driving modes."""

from copy import deepcopy

import numpy as np
import pytest
import xarray as xr

from carculator import (
    CarInputParameters,
    CarModel,
    InventoryCar,
    fill_xarray_from_input_parameters,
)

TARGETS = {
    "Mini": (8.8, 120),
    "Small": (8.8, 120),
    "Lower medium": (13, 128.375),
    "Medium": (12, 118.5),
    "Large": (12, 118.5),
    "Medium SUV": (18.1, 178.7375),
    "Large SUV": (24, 226.8),
    "Van": (13.6, 134.3),
}


def inputs(sizes=("Medium",), years=(2020,)):
    parameters = CarInputParameters()
    parameters.static()
    return fill_xarray_from_input_parameters(
        parameters,
        scope={
            "size": list(sizes),
            "powertrain": ["PHEV-p", "PHEV-d"],
            "year": list(years),
        },
    )[1]


def test_2020_calibration_matches_approved_nominal_capacity_and_pack_mass():
    array = inputs(TARGETS)
    original = array.copy(deep=True)
    model = CarModel(array, drop_hybrids=False)
    model.set_all()
    for size, (capacity, mass) in TARGETS.items():
        cell = model.array.sel(size=size, year=2020)
        # These independently transcribed source/proxy values include the whole
        # pack. A plausible kWh value alone must not conceal an undersized mass.
        np.testing.assert_allclose(
            cell.sel(parameter="electric energy stored"), capacity, rtol=2e-6
        )
        np.testing.assert_allclose(cell.sel(parameter="energy battery mass"), mass)
        np.testing.assert_allclose(
            cell.sel(parameter="battery cell mass")
            + cell.sel(parameter="battery BoP mass"),
            mass,
            rtol=2e-6,
        )
    xr.testing.assert_identical(array, original)


@pytest.mark.parametrize("mode", ["PHEV-e", "PHEV-c-p", "PHEV-p"])
def test_capacity_override_on_any_mode_reaches_both_modes_and_energy(mode):
    array = inputs()
    options = {"capacity": {(mode, "Medium", 2020): 24}}
    original = deepcopy(options)
    baseline = CarModel(array, drop_hybrids=False)
    baseline.set_all()
    model = CarModel(array, drop_hybrids=False, energy_storage=options)
    model.set_all()
    assert options == original
    np.testing.assert_allclose(model["electric energy stored"], 24, rtol=2e-6)
    np.testing.assert_allclose(model["energy battery mass"], 237, rtol=2e-6)
    for component in ("PHEV-e", "PHEV-c-p", "PHEV-c-d"):
        assert (
            model["driving mass"].sel(powertrain=component)
            > baseline["driving mass"].sel(powertrain=component)
        ).all()
        assert (
            model["TtW energy"].sel(powertrain=component)
            > baseline["TtW energy"].sel(powertrain=component)
        ).all()


@pytest.mark.parametrize(
    "kind,values", [("capacity", (12, 24)), ("electric", ("LFP", "NMC-622"))]
)
def test_conflicting_physical_packs_are_rejected(kind, values):
    with pytest.raises(ValueError, match="Conflicting PHEV battery"):
        CarModel(
            inputs(),
            energy_storage={
                kind: {
                    ("PHEV-e", "Medium", 2020): values[0],
                    ("PHEV-c-p", "Medium", 2020): values[1],
                }
            },
        )


def test_chemistry_override_preserves_one_pack_and_caller_choices():
    choices = {"electric": {("PHEV-c-d", "Medium", 2020): "LFP"}}
    original = deepcopy(choices)
    model = CarModel(inputs(), energy_storage=choices, drop_hybrids=False)
    model.set_all()
    assert choices == original
    for powertrain in model.array.powertrain.values:
        assert model.energy_storage["electric"][(powertrain, "Medium", 2020)] == "LFP"
    np.testing.assert_allclose(
        model["electric energy stored"],
        model["energy battery mass"]
        * model["battery cell mass share"]
        * model["battery cell energy density"],
        rtol=2e-6,
    )


def test_sampled_pack_is_shared_without_losing_years_or_samples():
    array = inputs(years=(2020, 2030)).sel(year=[2030, 2020])
    array = array.isel(value=[0, 0]).assign_coords(value=["reference", "heavier"])
    array.loc[
        dict(powertrain="PHEV-e", parameter="energy battery mass", value="heavier")
    ] *= 1.2
    model = CarModel(array, drop_hybrids=False)
    model.set_all()
    for parameter in (
        "energy battery mass",
        "battery cell mass",
        "electric energy stored",
    ):
        for mode in ("PHEV-c-p", "PHEV-c-d", "PHEV-p", "PHEV-d"):
            np.testing.assert_allclose(
                model[parameter].sel(powertrain=mode),
                model[parameter].sel(powertrain="PHEV-e"),
                rtol=2e-6,
            )
        np.testing.assert_allclose(
            model[parameter].sel(value="heavier")
            / model[parameter].sel(value="reference"),
            1.2,
            rtol=2e-6,
        )
    first = model.array.copy(deep=True)
    model.set_all()
    xr.testing.assert_allclose(first, model.array)


def test_calibrated_phev_inventory_and_lcia_are_finite():
    model = CarModel(inputs())
    model.set_all()
    inventory = InventoryCar(model, scenario="static", functional_unit="vkm")
    impacts = inventory.calculate_impacts()
    assert np.isfinite(impacts).all()
    assert impacts.sel(impact_category="climate change").sum().item() > 0


def test_combustion_mode_sustains_charge_and_purchases_all_fuel_energy():
    model = CarModel(inputs(TARGETS), drop_hybrids=False)
    model.set_all()
    for final, combustion in (("PHEV-p", "PHEV-c-p"), ("PHEV-d", "PHEV-c-d")):
        mode = model.array.sel(powertrain=combustion)
        purchased = (
            mode.sel(parameter="fuel consumption")
            * mode.sel(parameter="fuel density per kg")
            * mode.sel(parameter="LHV fuel MJ per kg")
            * 1000
        )
        np.testing.assert_allclose(
            purchased, mode.sel(parameter="TtW energy"), rtol=2e-6
        )
        # Carrying the full pack does not supply net grid energy in combustion
        # mode; the final vehicle receives electricity only through PHEV-e.
        assert (mode.sel(parameter="electric energy stored") > 0).all()
        assert (mode.sel(parameter="battery DoD") == 0).all()
        aggregate = model.array.sel(powertrain=final)
        electric = model.array.sel(powertrain="PHEV-e")
        grid = aggregate.sel(parameter="electricity consumption")
        fuel = (
            aggregate.sel(parameter="fuel consumption")
            * aggregate.sel(parameter="fuel density per kg")
            * aggregate.sel(parameter="LHV fuel MJ per kg")
            * 1000
        )
        stored = (
            grid
            * electric.sel(parameter="battery charge efficiency")
            * electric.sel(parameter="charger efficiency")
            * 3600
        )
        np.testing.assert_allclose(
            fuel + stored, aggregate.sel(parameter="TtW energy"), rtol=2e-6
        )
