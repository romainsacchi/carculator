from carculator_utils.background_systems import BackgroundSystemModel
import numpy as np

def test_size_dictionary():
    bs = BackgroundSystemModel()
    assert len(bs.electricity_mix) == 92
    assert int(bs.electricity_mix.sel(country="FR").values.sum()) == 10
    assert int(bs.electricity_mix.sel(country="FR", year=2015).values.sum()) == 1


def test_cumulative_losses():
    bs = BackgroundSystemModel()
    assert len(bs.losses) == 146
    assert float(bs.losses["AL"]["LV"]) > 1.1


def test_share_biodiesel():
    bs = BackgroundSystemModel()
    share = bs.get_share_biofuel("biodiesel", "DE", 2015)
    assert share == np.array(0.052300778)
    share = bs.get_share_biofuel("biodiesel", "DE", [2015, 2016, 2017])
    assert np.allclose(share, np.array([0.052300778, 0.050711505, 0.050703788]))
    share = bs.get_share_biofuel("biodiesel", "DE", 2045)
    assert share == np.array(0.054851109)

def test_fuel_blend():
    bs = BackgroundSystemModel()
    blend = bs.define_fuel_blends(["ICEV-d", "HEV-d", "ICEV-g", "ICEV-p"], "DE", [2015, 2016, 2017])

    assert isinstance(blend, dict)
    assert len(blend) == 3
    assert all(i in list(blend.keys()) for i in ["diesel", "cng", "petrol"])
    assert blend["diesel"]["primary"]["type"] == "diesel"
    assert np.allclose(blend["diesel"]["primary"]["share"], np.array([1 - 0.052300778, 1 - 0.050711505, 1 - 0.050703788]))
    assert blend["diesel"]["primary"]["lhv"] == 43
    assert blend["diesel"]["primary"]["CO2"] == 3.15

    assert blend["diesel"]["secondary"]["type"] == "biodiesel - cooking oil"

    assert np.allclose(
        blend["diesel"]["secondary"]["share"],
        np.array([0.052300778, 0.050711505, 0.050703788]),
        rtol=1e-3
    )
    assert blend["diesel"]["secondary"]["lhv"] == 38
    assert blend["diesel"]["secondary"]["CO2"] == 2.79


