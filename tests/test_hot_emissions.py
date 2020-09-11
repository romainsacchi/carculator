from carculator.driving_cycles import get_standard_driving_cycle
from carculator.hot_emissions import HotEmissionsModel
import numpy as np
import pytest

dc = get_standard_driving_cycle()
dc_name = "WLTC"


def test_wrong_powertrain():
    hem = HotEmissionsModel(dc, dc_name)
    with pytest.raises(TypeError) as wrapped_error:
        hem.get_emissions_per_powertrain("electric")
    assert wrapped_error.type == TypeError


def test_output_emissions():
    hem = HotEmissionsModel(dc, dc_name)
    em = hem.get_emissions_per_powertrain("diesel", euro_class=6.3)

    # Carbon monoxide emission, diesel
    assert em[:, 3:6, :, :].sum() > 0.000013
    assert em[:, 3:6, :, :].sum() < 0.000016

    # Particulate matter emission, diesel
    assert em[:, 9:12, :, :].sum() > 3.7e-7
    assert em[:, 9:12, :, :].sum() < 4.5e-7

    # Euro-6d NOx emission limit, diesel
    assert em[:, 6:9, :, :].sum() < 8e-5

    hem = HotEmissionsModel(dc, dc_name)
    em = hem.get_emissions_per_powertrain("petrol", euro_class=6.3)

    # Carbon monoxide emission, petrol
    assert em[:, 3:6, :, :].sum() > 2.5e-4
    assert em[:, 3:6, :, :].sum() < 3.5e-4

    # Particulate matter emission, petrol
    assert em[:, 9:12, :, :].sum() > 9.3e-7
    assert em[:, 9:12, :, :].sum() < 9.9e-7

    # Euro-6d NOx emission limit, petrol
    assert em[:, 6:9, :, :].sum() < 6e-5
