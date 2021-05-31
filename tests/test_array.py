from carculator import *
import numpy as np
import pytest


def test_type_cip():
    with pytest.raises(TypeError) as wrapped_error:
        fill_xarray_from_input_parameters("bla")
    assert wrapped_error.type == TypeError


def test_format_array():
    cip = CarInputParameters()
    cip.static()
    dcts, array = fill_xarray_from_input_parameters(cip)

    assert np.shape(array)[0] == len(dcts[0])
    assert np.shape(array)[1] == len(dcts[1])
    assert np.shape(array)[2] == len(dcts[2])
    assert np.shape(array)[3] == len(dcts[3])


def test_modify_array():
    cip = CarInputParameters()
    cip.static()
    _, array = fill_xarray_from_input_parameters(cip)

    dict_param = {
        ("Driving", "all", "all", "lifetime kilometers", "none"): {
            (2020, "loc"): 150000,
            (2040, "loc"): 150000,
        }
    }

    modify_xarray_from_custom_parameters(dict_param, array)
    assert (
        array.sel(
            powertrain="ICEV-d",
            size="Large",
            year=2020,
            parameter="lifetime kilometers",
        ).sum()
        == 150000
    )


def test_wrong_param_modify_array():
    cip = CarInputParameters()
    cip.static()
    _, array = fill_xarray_from_input_parameters(cip)

    dict_param = {
        ("Driving", "all", "all", "foo", "none"): {
            (2020, "loc"): 150000,
            (2040, "loc"): 150000,
        }
    }

    modify_xarray_from_custom_parameters(dict_param, array)
    with pytest.raises(KeyError) as wrapped_error:
        array.sel(powertrain="ICEV-d", size="Large", year=2020, parameter="foo")
    assert wrapped_error.type == KeyError

def test_scope():
    """ Test that the use of scope dictionary works as intended"""
    cip = CarInputParameters()
    cip.static()
    scope = {"powertrain":["ICEV-d"], "size": ["Lower medium"]}
    _, array = fill_xarray_from_input_parameters(cip, scope=scope)

    assert "BEV" not in array.coords["powertrain"].values
    assert "Large" not in array.coords["size"].values
