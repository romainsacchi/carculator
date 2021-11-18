import pytest
import numpy as np

from carculator.driving_cycles import get_standard_driving_cycle
from carculator.noise_emissions import NoiseEmissionsModel



def test_wrong_powertrain():
    dc = get_standard_driving_cycle()
    dc = np.repeat(dc.reshape(1, -1), 8, axis=0).T
    dc_name = "WLTC"


    nem = NoiseEmissionsModel(dc, dc_name)
    with pytest.raises(TypeError) as wrapped_error:
        nem.get_sound_power_per_compartment("foo")
    assert wrapped_error.type == TypeError


def test_output_emissions():
    dc = get_standard_driving_cycle()
    dc = np.repeat(dc.reshape(1, -1), 8, axis=0).T
    dc_name = "WLTC"


    nem = NoiseEmissionsModel(dc, dc_name)
    noise = nem.get_sound_power_per_compartment("combustion")

    assert len(np.squeeze(noise[0])) == 24
    assert noise[0].sum() > 0.9
    assert noise[0].sum() < 1.1



def test_urban_output_emissions():
    dc2 = get_standard_driving_cycle()
    dc2 = np.repeat(dc2.reshape(1, -1), 8, axis=0).T
    dc2_name = "WLTC 3.1"  # <- urban driving cycle

    nem = NoiseEmissionsModel(dc2, dc2_name)
    urban_noise = nem.get_sound_power_per_compartment("combustion")

    assert urban_noise[0, 0, 9:].sum() == 0