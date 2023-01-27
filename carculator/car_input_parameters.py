"""
car_input_parameters.py contains the VehicleInputParameters class, which, after instantiated, contains
all definitions (metadata) of input and calculated parameters, along with their values.
Also, it inherits methods from `klausen`, which exposes the methods .static() and .stochastic(),
which generate single or random values for input parameters.
"""

import json
from pathlib import Path
from typing import Union

from carculator_utils.vehicle_input_parameters import VehicleInputParameters


def load_parameters(obj: Union[str, Path, list]) -> list:
    """
    Returns a json object containing parameters' definitions
    :param obj: A filepath to a json file, or a json object
    :return: Returns a json object containing parameters' definitions
    """
    if isinstance(obj, (str, Path)):
        assert Path(obj).exists(), f"Can't find this filepath {obj}."
        return json.load(open(obj, encoding="utf-8"))

    # Already in correct form, just return
    return obj


class CarInputParameters(VehicleInputParameters):
    """ """

    DEFAULT = Path(__file__, "..").resolve() / "data" / "default_parameters.json"
    EXTRA = Path(__file__, "..").resolve() / "data" / "extra_parameters.json"

    def __init__(
        self,
        parameters: Union[str, Path, list] = None,
        extra: Union[str, Path, list] = None,
    ) -> None:
        """Create a `klausen <https://github.com/cmutel/klausen>`__ model with the car input parameters."""
        super().__init__(None)
