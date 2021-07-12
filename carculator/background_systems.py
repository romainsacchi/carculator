import numpy as np
import pandas as pd
import xarray as xr

from . import DATA_DIR


def data_to_dict(csv_list):
    """
    Returns a dictionary from a sequence of items.
    :param data: list
    :return: dict
    """

    (_, *header), *data = csv_list
    csv_dict = {}
    for row in data:
        key, *values = row
        csv_dict[key] = {key: value for key, value in zip(header, values)}

    return csv_dict


def get_electricity_losses():
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

    return data_to_dict(csv_list)


def get_region_mapping():
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

    return data_to_dict(csv_list)


def get_electricity_mix():
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

    array = (
        df.melt(id_vars=["country", "year"], value_name="value")
        .groupby(["country", "year", "variable"])["value"]
        .mean()
        .to_xarray()
    )
    array = array.interpolate_na(
        dim="year", method="linear", fill_value="extrapolate"
    ).clip(0, 1)
    array /= array.sum(axis=2)

    return array


def get_biofuel_share():
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
            val = df.loc[(df["Region"] == r) & (df["Scenario"] == s), "Biomass fuel":]
            array.loc[dict(region=r, scenario=s, value=0)] = val
    return array


def get_biomethane_share():
    filename = "share_bio_cng.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biomethane shares could not be found."
        )
    df = pd.read_csv(filepath, sep=";")

    return df.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_biogasoline_share():
    filename = "share_bio_gasoline.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biogasoline shares could not be found."
        )
    df = pd.read_csv(filepath, sep=";")

    return df.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_biodiesel_share():
    filename = "share_bio_diesel.csv"
    filepath = DATA_DIR / filename

    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biodiesel shares could not be found."
        )
    df = pd.read_csv(filepath, sep=";")

    return df.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_sulfur_content_in_fuel():
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
    df = pd.read_csv(filepath, sep=";")
    df = df.groupby(["country", "year"]).sum().unstack()
    df.loc[:, ("diesel", 1990)] = df["diesel"].max(1)
    df.loc[:, ("petrol", 1990)] = df["petrol"].max(1)

    df.loc[df[("diesel", 2019)] > 50 / 1e6, ("diesel", 2050)] = 50 / 1e6
    df.loc[df[("petrol", 2019)] > 50 / 1e6, ("petrol", 2050)] = 50 / 1e6
    df.loc[:, ("diesel", 2050)] = df["diesel"].min(1)
    df.loc[:, ("petrol", 2050)] = df["petrol"].min(1)

    df = df.interpolate(axis=1)
    df = df.unstack().reset_index()
    df = df.rename(columns={"level_0": "fuel"})
    arr = df.groupby(["country", "year", "fuel"]).sum()[0].to_xarray()
    return arr


class BackgroundSystemModel:
    """
    Retrieve and build dictionaries that contain important information to model in the background system:

        * gross electricity production mixes from nearly all countries in the world, from 2015 to 2050.
        * cumulative electricity transformation/transmission/distribution losses from high voltage to medium and low voltage.
        * share of biomass-derived fuel in the total consumption of liquid fuel in the transport sector. Source: REMIND.
    """

    def __init__(self):
        self.electricity_mix = get_electricity_mix()
        self.losses = get_electricity_losses()
        self.region_map = get_region_mapping()
        self.biofuel = get_biofuel_share()
        self.sulfur = get_sulfur_content_in_fuel()
        self.biomethane = get_biomethane_share()
        self.biogasoline = get_biogasoline_share()
        self.biodiesel = get_biodiesel_share()
