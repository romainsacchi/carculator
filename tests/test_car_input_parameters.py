import json
from pathlib import Path

import pytest

import carculator.car_input_parameters as cip

DEFAULT = Path(__file__, "..").resolve() / "fixtures" / "default_test.json"
EXTRA = Path(__file__, "..").resolve() / "fixtures" / "extra_test.json"


def test_retrieve_list_powertrains():
    assert isinstance(cip.CarInputParameters().powertrains, list)
    assert len(cip.CarInputParameters().powertrains) > 5


def test_can_pass_directly():
    d, e = json.load(open(DEFAULT)), set(json.load(open(EXTRA)))
    e.remove("foobazzle")
    assert len(cip.CarInputParameters(d, e).powertrains) == 5
    assert len(cip.CarInputParameters(d, e).parameters) == 12


def test_alternate_filepath():
    assert len(cip.CarInputParameters(DEFAULT, EXTRA).powertrains) == 5
    assert len(cip.CarInputParameters(DEFAULT, EXTRA).parameters) == 13
