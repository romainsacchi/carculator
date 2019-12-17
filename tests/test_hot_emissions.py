import pytest
from carculator.hot_emissions import HotEmissionsModel
from carculator.driving_cycles import get_standard_driving_cycle

def test_wrong_powertrain():
    dc = get_standard_driving_cycle()
    hem = HotEmissionsModel(dc)
    with pytest.raises(TypeError) as wrapped_error:
        hem.get_emissions_per_powertrain("electric")
    assert wrapped_error.type == TypeError

def test_output_emissions():
    dc = get_standard_driving_cycle()
    hem = HotEmissionsModel(dc)
    urban = hem.get_emissions_per_powertrain("diesel")[0]

    assert len(urban) == 11
    assert urban.sum() > 2.4e-5
    assert urban.sum() < 2.5e-5

