import pytest
from carculator.noise_emissions import NoiseEmissionsModel
from carculator.driving_cycles import get_standard_driving_cycle

def test_wrong_powertrain():
    dc = get_standard_driving_cycle()
    nem = NoiseEmissionsModel(dc)
    with pytest.raises(TypeError) as wrapped_error:
        nem.get_sound_power_per_compartment("foo")
    assert wrapped_error.type == TypeError

def test_output_emissions():
    dc = get_standard_driving_cycle()
    nem = NoiseEmissionsModel(dc)
    urban = nem.get_sound_power_per_compartment("combustion")[0]

    assert len(urban) == 8
    assert urban.sum() > 0.1
    assert urban.sum() < 0.2

