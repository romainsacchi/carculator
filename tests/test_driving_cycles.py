import os, sys;
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import pandas
import pytest
from carculator.driving_cycles import get_standard_driving_cycle



def test_cycle_retrival_nedc():
    dc = get_standard_driving_cycle("NEDC")
    assert isinstance(dc, pandas.core.series.Series)
    assert dc.sum() == 39353.0
    assert dc.index.min() == 0
    assert dc.index.max() == 1200

def test_cycle_retrival_cadc():
    dc = get_standard_driving_cycle("CADC")
    assert isinstance(dc, pandas.core.series.Series)
    assert dc.sum() == 186074.2
    assert dc.index.min() == 0
    assert dc.index.max() == 3143

def test_missing_cycle():
    with pytest.raises(KeyError):
        get_standard_driving_cycle("Foo")
