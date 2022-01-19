"""
car_input_parameters.py contains the CarInputParameters class, which, after instantiated, contains
all definitions (metadata) of input and calculated parameters, along with their values.
Also, it inherits methods from `klausen`, which exposes the methods .static() and .stochastic(),
which generate single or random values for input parameters.
"""

import json
from pathlib import Path
from typing import Union

from klausen import NamedParameters

DEFAULT = Path(__file__, "..").resolve() / "data" / "default_parameters.json"
EXTRA = Path(__file__, "..").resolve() / "data" / "extra_parameters.json"


def load_parameters(obj: Union[str, Path, list]) -> list:
    """
    Returns a json object containing parameters' definitions
    :param obj: A filepath to a json file, or a json object
    :return: Returns a json object containing parameters' definitions
    """
    if isinstance(obj, (str, Path)):
        assert Path(obj).exists(), "Can't find this filepath"
        return json.load(open(obj, encoding="utf-8"))

    # Already in correct form, just return
    return obj


class CarInputParameters(NamedParameters):
    """
    A class used to represent vehicles with associated type, size, technology, year and parameters.

    This class inherits from NamedParameters, located in the *klausen* package.
    It sources default parameters for all vehicle types from a dictionary in
    default_parameters and format them into an array following the structured described
    in the *klausen* package.

    :ivar sizes: List of string items e.g., ['Large', 'Lower medium', 'Medium', 'Mini', 'Medium SUV', 'Small', 'Van']
    :ivar powertrains: List of string items e.g., ['BEV', 'FCEV', 'HEV-p', 'ICEV-d', 'ICEV-g', 'ICEV-p', 'PHEV-c', 'PHEV-e']
    :ivar parameters: json or filepath to json containing input parameter's definitions
    :ivar extra: json or filepath to json containing calculated parameter's definitions
    :ivar years: List of integers e.g., [2000, 2010, 2020, 2040]
    :ivar metadata: Dictionary for metadata.
    :ivar values: Dictionary for storing values, of format {'param':[value]}.
    :ivar iterations: Number of iterations executed by the method :func:`~car_input_parameters.CarInputParameters.stochastic`.
        None if :func:`~car_input_parameters.CarInputParameters.static` used instead.
    :vartype iterations: int


    """

    def __init__(
        self,
        parameters: Union[str, Path, list] = None,
        extra: Union[str, Path, list] = None,
    ) -> None:
        """Create a `klausen <https://github.com/cmutel/klausen>`__ model with the car input parameters."""
        super().__init__(None)

        parameters = load_parameters(DEFAULT if parameters is None else parameters)
        extra = set(load_parameters(EXTRA if extra is None else extra))

        if not isinstance(parameters, dict):
            raise ValueError(
                f"Parameters are not correct type (expected `dict`, got `{type(parameters)}`)"
            )
        if not isinstance(extra, set):
            raise ValueError(
                f"Extra parameters are not correct type (expected `set`, got `{type(extra)}`)"
            )
        self.sizes = sorted(
            {size for o in parameters.values() for size in o.get("sizes", [])}
        )
        self.powertrains = sorted(
            {pt for o in parameters.values() for pt in o.get("powertrain", [])}
        )
        self.parameters = sorted(
            {o["name"] for o in parameters.values()}.union(set(extra))
        )

        # keep a list of input parameters, for sensitivity purpose
        self.input_parameters = sorted({o["name"] for o in parameters.values()})

        self.years = sorted({o["year"] for o in parameters.values()})
        self.add_car_parameters(parameters)

    def add_car_parameters(self, parameters: dict) -> None:
        """
        Split data and metadata according to ``klausen`` convention.

        The parameters are split into the *metadata* and *values* attributes
        of the CarInputParameters class by the add_parameters() method of the parent class.

        :param parameters: A dictionary that contains parameters.
        :type parameters: dict


        """
        keys = {"kind", "uncertainty_type", "amount", "loc", "minimum", "maximum"}

        reformatted = {}
        for key, dct in parameters.items():
            reformatted[key] = {k: v for k, v in dct.items() if k in keys}
            reformatted[key]["metadata"] = {
                k: v for k, v in dct.items() if k not in keys
            }

        self.add_parameters(reformatted)
