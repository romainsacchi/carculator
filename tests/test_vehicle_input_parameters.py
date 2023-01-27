import json
from pathlib import Path

import carculator.car_input_parameters as vip

DEFAULT = Path(__file__, "..").resolve() / "fixtures" / "default_test.json"
EXTRA = Path(__file__, "..").resolve() / "fixtures" / "extra_test.json"


def test_can_pass_directly():
    d, e = json.load(open(DEFAULT)), set(json.load(open(EXTRA)))
    e.remove("foobazzle")
    assert len(vip.CarInputParameters(d, e).powertrains) == 12
    assert len(vip.CarInputParameters(d, e).parameters) == 333


def test_alternate_filepath():
    assert len(vip.CarInputParameters(DEFAULT, EXTRA).powertrains) == 12
    assert len(vip.CarInputParameters(DEFAULT, EXTRA).parameters) == 333
