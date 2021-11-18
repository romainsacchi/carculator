import pytest

from carculator.driving_cycles import get_standard_driving_cycle
from carculator.hot_emissions import HotEmissionsModel

dc = get_standard_driving_cycle()
dc_name = "WLTC"


def test_wrong_powertrain():
    hem = HotEmissionsModel(dc, dc_name)
    with pytest.raises(TypeError) as wrapped_error:
        hem.get_hot_emissions(
            powertrain_type=["ICEV-d", "BEV"],
            euro_class=[5],
            lifetime_km=[200000, 200000],
            energy_consumption=[[10, 10, 10], [10, 10, 10]],
            yearly_km=[12000, 12000],
        )
    assert wrapped_error.type == TypeError
