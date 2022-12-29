import numpy as np
import pytest

from carculator import (
    VehicleInputParameters,
    fill_xarray_from_input_parameters,
)
from carculator.model import CarModel


def test_type_cip():
    with pytest.raises(TypeError) as wrapped_error:
        fill_xarray_from_input_parameters("bla")
    assert wrapped_error.type == TypeError


def test_format_array():
    cip = VehicleInputParameters()
    cip.static()
    dcts, array = fill_xarray_from_input_parameters(cip)

    assert np.shape(array)[0] == len(dcts[0])
    assert np.shape(array)[1] == len(dcts[1])
    assert np.shape(array)[2] == len(dcts[2])
    assert np.shape(array)[3] == len(dcts[3])


def test_modify_array():
    cip = VehicleInputParameters()
    cip.static()
    _, array = fill_xarray_from_input_parameters(cip)

    array.loc[dict(
        powertrain="ICEV-d",
        size="Large",
        year=2020,
        parameter="lifetime kilometers",
    )] = 150000

    cm = CarModel(array)
    cm.set_all()
    assert cm.array.sel(
        powertrain="ICEV-d",
        size="Large",
        year=2020,
        parameter="lifetime kilometers",
    ).values == 150000


def test_scope():
    """Test that the use of scope dictionary works as intended"""
    cip = VehicleInputParameters()
    cip.static()
    scope = {"powertrain": ["ICEV-d"], "size": ["Lower medium"]}
    _, array = fill_xarray_from_input_parameters(cip, scope=scope)

    assert "BEV" not in array.coords["powertrain"].values
    assert "Large" not in array.coords["size"].values
