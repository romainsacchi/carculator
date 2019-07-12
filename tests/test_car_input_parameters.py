import pytest
import carculator.car_input_parameters as cip

def test_type_parameters():
    with pytest.raises(SystemExit) as wrapped_error:
        cip.CarInputParameters("Foo", 28)
    assert wrapped_error.type == SystemExit
    assert wrapped_error.value.code ==1

def test_retrieve_list_powertrains():
    assert isinstance(cip.CarInputParameters().powertrains, list)
    assert len(cip.CarInputParameters().powertrains)>0