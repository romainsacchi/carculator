import numpy as np
import pytest

from carculator.driving_cycles import get_standard_driving_cycle


def test_cycle_retrieval_default():
    dc = get_standard_driving_cycle()
    assert isinstance(dc, np.ndarray)
    assert dc.sum() == 83744.6


def test_cycle_retrieval_wltc():
    dc = get_standard_driving_cycle("WLTC")
    assert isinstance(dc, np.ndarray)
    assert dc.sum() == 83744.6


def test_cycle_retrieval_nedc():
    dc = get_standard_driving_cycle("NEDC")
    assert isinstance(dc, np.ndarray)
    assert dc.sum() == 39353.0


def test_cycle_retrieval_cadc():
    dc = get_standard_driving_cycle("CADC")
    assert isinstance(dc, np.ndarray)
    assert dc.sum() == 186074.2


def test_missing_cycle():
    with pytest.raises(SystemExit) as wrapped_error:
        get_standard_driving_cycle("Foo")
    assert wrapped_error.type == SystemExit
    assert wrapped_error.value.code == 1
