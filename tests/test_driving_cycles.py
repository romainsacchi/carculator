import sys, os

import pandas
import pytest
from carculator.driving_cycles import get_standard_driving_cycle


def test_cycle_retrival_default():
    dc = get_standard_driving_cycle()
    assert isinstance(dc, pandas.core.series.Series)
    assert dc.sum() == 88744.6
    assert dc.index.min() == 0
    assert dc.index.max() == 1801

def test_cycle_retrival_wltc():
    dc = get_standard_driving_cycle("WLTC")
    assert isinstance(dc, pandas.core.series.Series)
    assert dc.sum() == 88744.6
    assert dc.index.min() == 0
    assert dc.index.max() == 1801


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
