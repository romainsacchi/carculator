import pytest

from carculator.driving_cycles import get_standard_driving_cycle
from carculator.noise_emissions import NoiseEmissionsModel

dc = get_standard_driving_cycle()
dc_name = "WLTC"


def test_wrong_powertrain():
    nem = NoiseEmissionsModel(dc, dc_name)
    with pytest.raises(TypeError) as wrapped_error:
        nem.get_sound_power_per_compartment("foo")
    assert wrapped_error.type == TypeError


def test_output_emissions():
    nem = NoiseEmissionsModel(dc, dc_name)
    urban = nem.get_sound_power_per_compartment("combustion")[0]

    assert len(urban) == 8
    assert urban.sum() > 0.09
    assert urban.sum() < 0.1
