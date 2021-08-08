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
    "BackgroundSystemModel",
    "ExportInventory",
    "InternalNoiseModel",
)
__version__ = (1, 5, 4)

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


from .array import (
    fill_xarray_from_input_parameters,
    modify_xarray_from_custom_parameters,
)
from .background_systems import BackgroundSystemModel
from .car_input_parameters import CarInputParameters
from .driving_cycles import get_standard_driving_cycle
from .export import ExportInventory
from .hot_emissions import HotEmissionsModel
from .internal_noise import InternalNoiseModel
from .inventory import InventoryCalculation
from .model import CarModel
from .noise_emissions import NoiseEmissionsModel
