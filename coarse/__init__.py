__all__ = (
    "get_standard_driving_cycle",
    "CarInputParameters",
)
__version__ = (0, 0, 1)


from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


from .driving_cycles import get_standard_driving_cycle
from .car_input_parameters import CarInputParameters
