"""

Submodules
==========

.. autosummary::
    :toctree: _autosummary


"""

__all__ = (
    "CarInputParameters",
    "fill_xarray_from_input_parameters",
    "CarModel",
    "InventoryCar",
)
__version__ = (1, 8, 0)

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"

from carculator_utils.array import fill_xarray_from_input_parameters

from .car_input_parameters import CarInputParameters
from .inventory import InventoryCar
from .model import CarModel
