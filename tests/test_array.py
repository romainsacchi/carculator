import numpy as np
import pytest
from carculator import *

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