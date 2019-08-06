import pytest
import carculator.car_input_parameters as cip

def test_retrieve_list_powertrains():
    assert isinstance(cip.CarInputParameters().powertrains, list)
    assert len(cip.CarInputParameters().powertrains)>0