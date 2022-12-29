import numpy as np
import pytest

from carculator_utils.driving_cycles import get_standard_driving_cycle_and_gradient


def test_cycle_retrieval_default_car():
    dc, grad = get_standard_driving_cycle_and_gradient(
        vehicle_type="car",
        vehicle_sizes=["Small", "Medium", "Large"],
        name="WLTC",
    )
    assert isinstance(dc, np.ndarray)
    assert dc.shape == (3144, 3)
    assert dc.shape == grad.shape

    dc, grad = get_standard_driving_cycle_and_gradient(
        vehicle_type="car",
        vehicle_sizes=["Small", "Medium", "Large"],
        name="NEDC",
    )
    assert isinstance(dc, np.ndarray)
    assert dc.shape == (3144, 3)
    assert dc.shape == grad.shape

def test_cycle_retrieval_default_bus():

    with pytest.raises(KeyError) as wrapped_error:
        dc, grad = get_standard_driving_cycle_and_gradient(
            vehicle_type="bus",
            vehicle_sizes=["Small", "Medium", "Large"],
            name="bus"
        )

    dc, grad = get_standard_driving_cycle_and_gradient(
        vehicle_type="bus",
        vehicle_sizes=["13m-city", "13m-coach-double", "18m"],
        name="bus"
    )
    assert isinstance(dc, np.ndarray)
    assert dc.shape == (17918, 3)
    assert dc.shape == grad.shape

def test_missing_cycle():
    with pytest.raises(KeyError) as wrapped_error:

        get_standard_driving_cycle_and_gradient(
            vehicle_type="car",
            vehicle_sizes=["Small", "Medium", "Large"],
            name="not a cycle",
        )
    assert wrapped_error.type == KeyError
