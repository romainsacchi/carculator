from . import DATA_DIR
import pandas as pd
import sys


def get_standard_driving_cycle(name="WLTC"):

    """Get driving cycle data as a Pandas `Series`.

    Driving cycles are given as km/h per second up to 3200 seconds.


    :param name: The name of the driving cycle. WLTC (Worldwide harmonized Light vehicles Test Cycles) is chosen by default if ``name`` left unspecified.
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
    try:
        return (
            pd.read_excel(DATA_DIR / "driving_cycles.xlsx", sheet_name="Driving cycles")
            .set_index("Time (s)")[name]
            .dropna()
        )
    except KeyError:
        print("The specified driving cycle could not be found.")
        sys.exit(1)


