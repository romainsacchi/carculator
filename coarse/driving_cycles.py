from . import DATA_DIR
import pandas as pd


def get_standard_driving_cycle(name):
    """Get driving cycle data as a Pandas ``Series``.

    Driving cycles are given as km/h per second up to 3200 seconds.

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

    """
    return pd.read_excel(DATA_DIR / "driving_cycles.xlsx", sheet_name='Driving cycles'
        ).set_index('Time (s)')[name].dropna()
