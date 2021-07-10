import sys

import numpy as np

from . import DATA_DIR


def get_standard_driving_cycle(name="WLTC"):

    """Get driving cycle data as a Pandas `Series`.

    Driving cycles are given as km/h per second up to 3200 seconds.

    :param name: The name of the driving cycle. WLTC (Worldwide harmonized Light vehicles Test Cycles) is chosen by default if :param name: left unspecified.
    :type name: str

    ``name`` should be one of:

    * WLTC
    * WLTC 3.1
    * WLTC 3.2
    * WLTC 3.3
    * WLTC 3.4
    * CADC Urban
    * CADC Road
    * CADC Motorway
    * CADC Motorway 130
    * CADC
    * NEDC

    :returns: A pandas DataFrame object with driving time (in seconds) as index,
        and velocity (in km/h) as values.
    :rtype: panda.Series


    """
    dict_dc_names = {
        "WLTC": 1,
        "WLTC 3.1": 2,
        "WLTC 3.2": 3,
        "WLTC 3.3": 4,
        "WLTC 3.4": 5,
        "CADC Urban": 6,
        "CADC Road": 7,
        "CADC Motorway": 8,
        "CADC Motorway 130": 9,
        "CADC": 10,
        "NEDC": 11,
        "SORT 1": 12,
        "SORT 2": 13,
        "SORT 3": 14,
    }
    try:
        arr = np.genfromtxt(DATA_DIR / "driving_cycles.csv", delimiter=";")
        dc = arr[1:, dict_dc_names[name]]
        dc = dc[~np.isnan(dc)]

        return dc
    except KeyError:
        print("The specified driving cycle could not be found.")
        sys.exit(1)
