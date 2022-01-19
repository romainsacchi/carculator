"""
background_systems.py contains BackgroundSystem, which contains different arrays relating to
the evolution of fuel blends, sulfur content in fuels, and electricity mixes in different
countries over time.
"""


from typing import Dict

import numpy as np
import pandas as pd
import xarray as xr

from . import DATA_DIR


def data_to_dict(csv_list: list) -> dict:
    """
    Returns a dictionary from a sequence of items.
    :param csv_list: list
    :return: dict
    """

    (_, *header), *data = csv_list
    csv_dict = {}
    for row in data:
        key, *values = row
        csv_dict[key] = dict(zip(header, values))

    return csv_dict


def get_electricity_losses() -> Dict[str, float]:
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
    with open(filepath, encoding="utf-8") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]

    return data_to_dict(csv_list)


def get_region_mapping() -> Dict[str, str]:
    """
    Retrieve mapping between ISO country codes and REMIND regions.

    :returns: dictionary
    :rtype: dict

    """
    filename = "region_mapping.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains correspondences between REMIND region "
            "names and ISO country codes could not be found."
        )
    with open(filepath, encoding="utf-8") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]

    return data_to_dict(csv_list)


def get_electricity_mix() -> xr.DataArray:
    """
    Retrieve electricity mixes and shape them into a xarray.
    Source:
        * for European countries (`ENTSOE TYNDP 2020 scenarios <https://2020.entsos-tyndp-scenarios.eu/>`_),
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

    dataframe = pd.read_csv(filepath, sep=";", index_col=["country", "year"])
    dataframe = dataframe.reset_index()

    array = (
        dataframe.melt(id_vars=["country", "year"], value_name="value")
        .groupby(["country", "year", "variable"])["value"]
        .mean()
        .to_xarray()
    )
    array = array.interpolate_na(
        dim="year", method="linear", fill_value="extrapolate"
    ).clip(0, 1)
    array /= array.sum(axis=2)

    return array


def get_biofuel_share() -> xr.DataArray:
    """
    Retrieve shares of biofuel consumption from REMIND and shape them into an xarray.

    :return: An axarray with 'country' and 'year' as dimensions
    :rtype: xarray.core.dataarray.DataArray
    """
    filename = "biofuel_share.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biofuel shares could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")

    country_codes = dataframe["Region"].unique()
    scenarios = dataframe["Scenario"].unique()
    year = dataframe["Year"].unique()
    tech = dataframe.columns[3:]
    array = xr.DataArray(
        np.zeros((len(country_codes), len(year), len(scenarios), len(tech), 1)),
        coords=[country_codes, year, scenarios, tech, np.arange(1)],
        dims=["region", "year", "scenario", "fuel_type", "value"],
    )
    for country_code in country_codes:
        for scenario in scenarios:
            val = dataframe.loc[
                (dataframe["Region"] == country_code)
                & (dataframe["Scenario"] == scenario),
                "Biomass fuel":,
            ]
            array.loc[dict(region=country_code, scenario=scenario, value=0)] = val
    return array


def get_biomethane_share() -> xr.DataArray:
    """

    :return: Returns an xarray with share of biomethane in the fuel blend, per country, per year.
    """
    filename = "share_bio_cng.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biomethane shares could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")

    return dataframe.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_biogasoline_share() -> xr.DataArray:
    """

    :return: Returns an xarray with share of bioethanol in the fuel blend, per country, per year.
    """
    filename = "share_bio_gasoline.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biogasoline shares could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")

    return dataframe.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_biodiesel_share() -> xr.DataArray:
    """

    :return: Returns an xarray with share of biodiesel in the fuel blend, per country, per year.
    """
    filename = "share_bio_diesel.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biodiesel shares could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")

    return dataframe.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_sulfur_content_in_fuel() -> xr.DataArray:
    """
    Retrieve sulfur content per kg of petrol and diesel.
    For CH, DE, FR, AU and SE, the concentration values come from HBEFA 4.1, from 1909 to 2020 (extrapolated to 2050).

    For the other countries, values come from
    Miller, J. D., & Jin, L. (2019). Global progress toward soot-free diesel vehicles in 2019.
    International Council on Clean Transportation.
    https://www.theicct.org/publications/global-progress-toward-soot-free-diesel-vehicles-2019

    There is an assumption made: countries that have high-sulfur content fuels (above 50 ppm in 2019) are assumed to
    improve over time to reach 50 ppm by 2050.

    :return: An axarray with 'country' and 'year' as dimensions
    :rtype: xarray.core.dataarray.DataArray
    """
    filename = "S_concentration_fuel.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains sulfur concentration values could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")
    dataframe = dataframe.groupby(["country", "year"]).sum().unstack()
    dataframe.loc[:, ("diesel", 1990)] = dataframe["diesel"].max(1)
    dataframe.loc[:, ("petrol", 1990)] = dataframe["petrol"].max(1)

    dataframe.loc[dataframe[("diesel", 2019)] > 50 / 1e6, ("diesel", 2050)] = 50 / 1e6
    dataframe.loc[dataframe[("petrol", 2019)] > 50 / 1e6, ("petrol", 2050)] = 50 / 1e6
    dataframe.loc[:, ("diesel", 2050)] = dataframe["diesel"].min(1)
    dataframe.loc[:, ("petrol", 2050)] = dataframe["petrol"].min(1)

    dataframe = dataframe.interpolate(axis=1)
    dataframe = dataframe.unstack().reset_index()
    dataframe = dataframe.rename(columns={"level_0": "fuel"})
    arr = dataframe.groupby(["country", "year", "fuel"]).sum()[0].to_xarray()
    return arr


class BackgroundSystemModel:
    """
    Retrieve and build dictionaries that contain important information to model in the background system:

        * gross electricity production mixes from nearly all countries in the world, from 2015 to 2050.
        * cumulative electricity transformation/transmission/distribution losses from high voltage to medium and low voltage.
        * share of biomass-derived fuel in the total consumption of liquid fuel in the transport sector. Source: REMIND.
        * share of bioethanol, biodiesel and biomethane, for each country, for different years.
        * share of sulfur in gasoline and diesel, for different countries and years.
    """

    def __init__(self) -> None:
        self.electricity_mix = get_electricity_mix()
        self.losses = get_electricity_losses()
        self.region_map = get_region_mapping()
        self.biofuel = get_biofuel_share()
        self.sulfur = get_sulfur_content_in_fuel()
        self.biomethane = get_biomethane_share()
        self.biogasoline = get_biogasoline_share()
        self.biodiesel = get_biodiesel_share()

    def __str__(self):
        return self.__class__.__name__
