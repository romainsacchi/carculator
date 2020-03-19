from carculator.background_systems import BackgroundSystemModel
import pytest


def test_size_dictionary():
    bs = BackgroundSystemModel()
    assert len(bs.electricity_mix) == 88
    assert int(bs.electricity_mix.sel(country="FR").values.sum()) == 7
    assert int(bs.electricity_mix.sel(country="FR", year=2015).values.sum()) == 1


def test_cumulative_losses():
    bs = BackgroundSystemModel()
    assert len(bs.losses) == 146
    assert float(bs.losses["AL"]["LV"]) > 1.1
