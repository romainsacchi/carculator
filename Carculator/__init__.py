"""

Submodules
==========

.. autosummary::
    :toctree: _autosummary


"""



__all__ = (
    "CarInputParameters",
    "fill_xarray_from_input_parameters",
    "get_standard_driving_cycle",
    "CarModel",
)
__version__ = (0, 0, 1)

# For relative imports to work in Python 3.6
import os, sys;
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


from .array import fill_xarray_from_input_parameters
from .car_input_parameters import CarInputParameters
from .driving_cycles import get_standard_driving_cycle
from .model import CarModel
