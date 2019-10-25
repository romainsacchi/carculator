"""

Submodules
==========

.. autosummary::
    :toctree: _autosummary


"""


__all__ = (
    "CarInputParameters",
    "fill_xarray_from_input_parameters",
    "modify_xarray_from_custom_parameters",
    "get_standard_driving_cycle",
    "CarModel",
    "NoiseEmissionsModel",
    "HotEmissionsModel",
    "InventoryCalculation",
    "BackgroundSystemModel"
)
__version__ = (0, 0, 1)

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


from .array import (
    fill_xarray_from_input_parameters,
    modify_xarray_from_custom_parameters,
)
from .car_input_parameters import CarInputParameters
from .noise_emissions import NoiseEmissionsModel
from .hot_emissions import HotEmissionsModel
from .driving_cycles import get_standard_driving_cycle
from .model import CarModel
from .inventory import InventoryCalculation
from .background_systems import BackgroundSystemModel
