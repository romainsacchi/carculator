"""
.. module: car_input_parameters.py

"""

from .default_parameters import DEFAULT, EXTRA
from klausen import NamedParameters

class CarInputParameters(NamedParameters):
    """
    A class used to represent vehicles with associated type, size, technology, year and parameters.

    This class inherits from NamedParameters, located in the *klausen* package.
    It sources default parameters for all vehicle types from a dictionary in
    default_parameters and format them into an array following the structured described
    in the *klausen* package.

    ...

    :param parameters: A dictionary that contains parameters.
        If left unspecified, default parameters found in default_parameters.py are used.
    :type parameters: dict
    :param extra: A dictionary that contains additional parameters.
        If left unspecified, default additional parameters found in default_parameters.py are used.
    :type extra: dict

    :ivar sizes: List of string items e.g., ['Large', 'Lower medium', 'Medium', 'Mini', 'SUV', 'Small', 'Van']
    :vartype sizes: list
    :ivar powertrains: List of string items e.g., ['BEV', 'FCEV', 'HEV-p', 'ICEV-d', 'ICEV-g', 'ICEV-p', 'PHEV-c', 'PHEV-e']
    :vartype powertrains: list
    :ivar parameters: List of string items e.g., ['Benzene', 'CH4', 'CNG tank mass intercept',...]
    :vartype parameters: list
    :ivar years: List of integers e.g., [2017, 2040]
    :vartype years: list
    :ivar metadata: Dictionary for metadata.
    :vartype metadata: dict
    :ivar values: Dictionary for storing values, of format {'param':[value]}.
    :vartype values: dict
    :ivar iterations: Number of iterations executed by the method stochastic().
        None if static() used instead.
    :vartype iterations: int


    """
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
        """
        Split data and metadata according to ``klausen`` convention.

        The parameters are split into the *metadata* and *values* attributes
        of the CarInputParameters class by the add_parameters() method of the parent class.

        :param parameters: A dictionary that contains parameters.
        If left unspecified, default parameters found in default_parameters.py are used.
        :type parameters: dict


        """
        KEYS = {"kind", "uncertainty_type", "amount", "loc", "minimum", "maximum"}

        reformatted = {}
        for key, dct in parameters.items():
            reformatted[key] = {k: v for k, v in dct.items() if k in KEYS}
            reformatted[key]["metadata"] = {
                k: v for k, v in dct.items() if k not in KEYS
            }

        self.add_parameters(reformatted)
