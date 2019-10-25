"""
.. module: background_systems.py

"""

import numpy as np
import pandas as pd
import xarray as xr
from . import DATA_DIR

class BackgroundSystemModel:
    """
    Retrieve and calculate important information to model in the background system: electricity mixes, etc.

    """

    def __init__(self):
        self.electricity_mix = self.get_electricity_mix()

    def get_electricity_mix(self):
        """
        Retreive electricity mixes and shape them into an xarray.
        Electricity mixes from EU-STEM for Europe, World Markal for outside Europe,
        from 2015 to 2050.

        :returns: An axarray with 'country' and 'year' as dimensions
        :rtype: xarray.core.dataarray.DataArray

        """
        filename = "electricity_mixes.csv"
        filepath = (DATA_DIR / filename)
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary that contains electricity mixes could not be found."
            )


        df = pd.read_csv(filepath, sep=";", index_col=['Country code', 'Year'])

        country_code = df.index.get_level_values('Country code').unique()
        year = df.index.get_level_values('Year').unique()
        tech = df.columns.unique()

        array = xr.DataArray(
            np.zeros((len(country_code), len(year), len(tech), 1)),
            coords=[country_code, year, tech, np.arange(1)],
            dims=["country", "year", "technology", "value"],
        )
        for r in country_code:
            val = df.loc[(df.index.get_level_values("Country code") == r), :]
            array.loc[dict(country=r, value=0)] = val

        return array

