"""
driving_cycles.py loads a driving cycle based on
the name specific by the user.
The driving cycle returned is a numpy
array with speed levels (in km/h) for each
second of driving.
"""

import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import yaml

from . import DATA_DIR

FILEPATH_DC_SPECS = DATA_DIR / "driving cycle" / "dc_specs.yaml"


def detect_vehicle_type(vehicle_sizes: List[str]) -> str:
    """
    Detect the type of vehicle based on the size of the vehicle.
    """

    dc = get_driving_cycle_specs()

    for vehicle_type in dc["columns"]:
        for dc_name in dc["columns"][vehicle_type]:
            for size in dc["columns"][vehicle_type][dc_name]:
                if size in vehicle_sizes:
                    return vehicle_type

    raise ValueError("The vehicle size is not in the list of available vehicle sizes.")


def get_driving_cycle_specs() -> dict:
    """Get driving cycle specifications.

    :returns: A dictionary with driving cycle specifications.
    :rtype: dict

    """

    with open(FILEPATH_DC_SPECS, "r") as f:
        return yaml.safe_load(f)


def get_dc_column_number(
    vehicle_type: str, vehicle_size: List[str], dc_name: str
) -> List[int]:
    """
    Loads YAML file that contains the column number.
    Return the column number given a vehicle type and driving cycle name.
    """

    dc_specs = get_driving_cycle_specs()

    if vehicle_type not in dc_specs["columns"]:
        raise KeyError(
            f"Vehicle type {vehicle_type} is not in the list of "
            f"available vehicle types: {list(dc_specs['columns'].keys())}"
        )

    if dc_name not in dc_specs["columns"][vehicle_type]:
        raise KeyError(
            f"Driving cycle {dc_name} is not in the list of "
            f"available driving cycles: {list(dc_specs['columns'][vehicle_type].keys())}"
        )

    if not all(
        vehicle in dc_specs["columns"][vehicle_type][dc_name]
        for vehicle in vehicle_size
    ):
        raise KeyError(
            f"Vehicle size(s) {vehicle_size} is not in the list of "
            f"available vehicle sizes: {list(dc_specs['columns'][vehicle_type][dc_name].keys())}"
        )

    return [dc_specs["columns"][vehicle_type][dc_name][s] for s in vehicle_size]


def get_data(
    filepath: Path, vehicle_type: str, vehicle_sizes: List[str], name: str
) -> np.ndarray:

    try:
        col = get_dc_column_number(vehicle_type, vehicle_sizes, name)
        arr = np.genfromtxt(filepath, delimiter=";")
        # we skip the headers
        dc = arr[1:, col]
        return dc

    except KeyError as err:
        print(err, "The specified driving cycle could not be found.")
        raise


def get_standard_driving_cycle_and_gradient(
    vehicle_type: str, vehicle_sizes: List[str], name: str
) -> Tuple[np.ndarray, np.ndarray]:

    """Get driving cycle and gradient data as a Pandas `Series`.

    Driving cycles are given as km/h per second up to 3200 seconds.

    :param name: The name of the driving cycle.
    e.g., WLTC (Worldwide harmonized Light vehicles Test Cycles)
    :type name: str

    :returns: A pandas DataFrame object with driving time
    (in seconds) as index, and velocity (in km/h) as values.
    :rtype: panda.Series

    """

    filepath_dc = DATA_DIR / "driving cycle" / f"{vehicle_type}.csv"
    filepath_gradient = DATA_DIR / "gradient" / f"{vehicle_type}.csv"

    # definition of columns to select in the CSV file
    # each column corresponds to a size class
    # since the driving cycle is simulated for each size class
    return (
        get_data(filepath_dc, vehicle_type, vehicle_sizes, name),
        get_data(filepath_gradient, vehicle_type, vehicle_sizes, name),
    )
