"""Run the native sizing loop to guard FCEV power and pack accounting."""

import numpy as np

from carculator import CarInputParameters, CarModel, fill_xarray_from_input_parameters


def test_completed_fcev_retains_complementary_power_and_complete_pack():
    inputs = CarInputParameters()
    inputs.static()
    _, array = fill_xarray_from_input_parameters(
        inputs, scope={"size": ["Medium"], "powertrain": ["FCEV"], "year": [2020, 2050]}
    )
    model = CarModel(array, country="CH")
    model.set_all()
    net_stack = model["fuel cell power"] / model["fuel cell own consumption"]
    np.testing.assert_allclose(
        net_stack + model["battery power"], model["power"], rtol=2e-6
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
        model["power"],
        rtol=2e-6,
    )
