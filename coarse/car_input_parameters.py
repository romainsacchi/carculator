from .default_parameters import DEFAULT, EXTRA
from klausen import NamedParameters
import itertools
import numpy as np


class CarInputParameters(NamedParameters):
    def __init__(self, parameters=None, extra=None):
        """Create a `klausen <https://github.com/cmutel/klausen>`__ model with the car input parameters."""
        super().__init__(None)

        if parameters is None:
            parameters = DEFAULT
        if extra is None:
            extra = EXTRA

        self.sizes = sorted(
            {size for o in parameters.values() for size in o.get("sizes", [])}
        )
        self.powertrains = sorted(
            {pt for o in parameters.values() for pt in o.get("powertrain", [])}
        )
        self.parameters = sorted(
            {o["name"] for o in parameters.values()}.union(set(extra))
        )
        self.years = sorted({o["year"] for o in parameters.values()})

        self.add_car_parameters(parameters)

    def add_car_parameters(self, parameters):
        """Split data and metadata according to ``klausen`` convention."""
        KEYS = {"kind", "uncertainty_type", "amount", "loc", "minimum", "maximum"}

        reformatted = {}
        for key, dct in parameters.items():
            reformatted[key] = {k: v for k, v in dct.items() if k in KEYS}
            reformatted[key]["metadata"] = {
                k: v for k, v in dct.items() if k not in KEYS
            }

        self.add_parameters(reformatted)
