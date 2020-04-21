from . import DATA_DIR
import numpy as np
import pandas as pd
import xarray as xr


class BackgroundSystemModel:
    """
    Retrieve and build dictionaries that contain important information to model in the background system:

        * gross electricity production mixes from nearly all countries in the world, from 2015 to 2050.
        * cumulative electricity transformation/transmission/distribution losses from high voltage to medium and low voltage.
        * share of biomass-derived fuel in the total consumption of liquid fuel in the transport sector. Source: REMIND.
    """

    def __init__(self):
        self.electricity_mix = self.get_electricity_mix()
        self.losses = self.get_electricity_losses()
        self.region_map = self.get_region_mapping()
        self.biofuel = self.get_biofuel_share()

    def get_electricity_losses(self):
        """
        Retrieve cumulative electricity losses from high to medium and low voltage.
        Source: `ecoinvent v.3.6 <https://www.ecoinvent.org/>`_.

        :returns: dictionary
        :rtype: dict

        """
        filename = "cumulative_electricity_losses.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The CSV file that contains electricity mixes could not be found."
            )
        with open(filepath) as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]

        (_, *header), *data = csv_list
        csv_dict = {}
        for row in data:
            key, *values = row
            csv_dict[key] = {key: value for key, value in zip(header, values)}
        return csv_dict

    def get_region_mapping(self):
        """
        Retrieve mapping between ISO country codes and REMIND regions.

        :returns: dictionary
        :rtype: dict

        """
        filename = "region_mapping.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The CSV file that contains correspondences between REMIND region names and ISO country codes "
                "could not be found."
            )
        with open(filepath) as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]

        (_, *header), *data = csv_list
        csv_dict = {}
        for row in data:
            key, *values = row
            csv_dict[key] = {key: value for key, value in zip(header, values)}
        return csv_dict

    def get_electricity_mix(self):
        """
        Retrieve electricity mixes and shape them into an xarray.
        Source:
            * for European countries (`EU Reference Scenario 2016 <https://ec.europa.eu/energy/en/data-analysis/energy-modelling/eu-reference-scenario-2016>`_),
            * for African countries (`TEMBA <http://www.osemosys.org/temba.html>`_ model)
            * and for other countries (`IEA World Energy outlook 2017 <https://www.iea.org/reports/world-energy-outlook-2017>`_)

        :returns: An axarray with 'country' and 'year' as dimensions
        :rtype: xarray.core.dataarray.DataArray

        """
        filename = "electricity_mixes.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The CSV file that contains electricity mixes could not be found."
            )

        df = pd.read_csv(filepath, sep=";", index_col=["country", "year"])
        df = df.reset_index()

        array = df.melt(id_vars=["country", "year"],
                    value_name="value")\
                    .groupby(["country","year", "variable"])["value"]\
                    .mean()\
                    .to_xarray()
        array = array.interpolate_na(dim="year", method="linear", fill_value='extrapolate').clip(0,1)
        array /= array.sum(axis=2)

        return array

    def get_biofuel_share(self):
        """
        Retrieve shares of biofuel consumption from REMIND and shape them into an xarray.

        :param country: Country to return the biofuel share for.
        :type country: str. 2-digit ISO country code.
        :return: An axarray with 'country' and 'year' as dimensions
        :rtype: xarray.core.dataarray.DataArray
        """
        filename = "biofuel_share.csv"
        filepath = DATA_DIR / filename

        if not filepath.is_file():
            raise FileNotFoundError(
                "The CSV file that contains biofuel shares could not be found."
            )
        df = pd.read_csv(filepath, sep=";")

        country_code = df["Region"].unique()
        scenario = df["Scenario"].unique()
        year = df["Year"].unique()
        tech = df.columns[3:]
        array = xr.DataArray(
            np.zeros((len(country_code), len(year), len(scenario), len(tech), 1)),
            coords=[country_code, year, scenario, tech, np.arange(1)],
            dims=["region", "year", "scenario", "fuel_type", "value"],
        )
        for r in country_code:
            for s in scenario:
                val = df.loc[
                    (df["Region"] == r) & (df["Scenario"] == s), "Biomass fuel":
                ]
                array.loc[dict(region=r, scenario=s, value=0)] = val
        return array
