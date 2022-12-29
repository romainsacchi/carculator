import json
from pathlib import Path

import carculator_utils.vehicle_input_parameters as vip

DEFAULT = Path(__file__, "..").resolve() / "fixtures" / "default_test.json"
EXTRA = Path(__file__, "..").resolve() / "fixtures" / "extra_test.json"


def test_can_pass_directly():
    d, e = json.load(open(DEFAULT)), set(json.load(open(EXTRA)))
    e.remove("foobazzle")
    assert len(vip.VehicleInputParameters(d, e).powertrains) == 5
    assert len(vip.VehicleInputParameters(d, e).parameters) == 12


def test_alternate_filepath():
    assert len(vip.VehicleInputParameters(DEFAULT, EXTRA).powertrains) == 5
    assert len(vip.VehicleInputParameters(DEFAULT, EXTRA).parameters) == 13
