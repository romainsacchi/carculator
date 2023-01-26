"""
background_systems.py contains BackgroundSystem, which contains different arrays relating to
the evolution of fuel blends, sulfur content in fuels, and electricity mixes in different
countries over time.
"""


from typing import Dict, List

import numpy as np
import pandas as pd
import xarray as xr
import yaml

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
    filepath = DATA_DIR / "electricity" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains electricity " "mixes could not be found."
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
    filepath = DATA_DIR / "electricity" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains " "electricity mixes could not " "be found."
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


def get_biofuel_share(filepath) -> xr.DataArray:
    """
    :return: Returns a xarray with share of biodiesel
    in the fuel blend, per country, per year.
    """
    if not filepath.is_file():
        raise FileNotFoundError(
            "The CSV file that contains biofuel shares " "shares could not be found."
        )
    dataframe = pd.read_csv(filepath, sep=";")

    return dataframe.groupby(["country", "year"]).sum().to_xarray().to_array()


def get_sulfur_content_in_fuel() -> xr.DataArray:
    """
    Retrieve sulfur content per kg of petrol and diesel.
    For CH, DE, FR, AU and SE, the concentration values come
    from HBEFA 4.1, from 1909 to 2020 (extrapolated to 2050).

    For the other countries, values come from
    Miller, J. D., & Jin, L. (2019). Global progress toward
    soot-free diesel vehicles in 2019.
    International Council on Clean Transportation.
    https://www.theicct.org/publications/global-progress-toward-soot-free-diesel-vehicles-2019

    There is an assumption made: countries that have
    high-sulfur content fuels (above 50 ppm in 2019) are assumed to
    improve over time to reach 50 ppm by 2050.

    :return: An axarray with 'country' and 'year' as dimensions
    :rtype: xarray.core.dataarray.DataArray
    """
    filename = "S_concentration_fuel.csv"
    filepath = DATA_DIR / "fuel" / filename

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


def get_default_fuels() -> dict:
    """
    Import default fuels from `default_fuels.yaml`
    """
    filename = "default_fuels.yaml"
    filepath = DATA_DIR / "fuel" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The YAML file that contains default fuels could not be found."
        )
    with open(filepath, "r") as file:
        default_fuels = yaml.safe_load(file)

    return default_fuels


def get_fuels_specs() -> dict:
    """
    Import fuel specifications from `fuel_specs.yaml`
    Contains names, LHV, CO2 emission factors.
    """
    with open(DATA_DIR / "fuel" / "fuel_specs.yaml", "r", encoding="utf-8") as stream:
        fuel_specs = yaml.safe_load(stream)

    return fuel_specs


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
        self.sulfur = get_sulfur_content_in_fuel()
        self.biomethane = get_biofuel_share(DATA_DIR / "fuel" / "share_bio_cng.csv")
        self.bioethanol = get_biofuel_share(
            DATA_DIR / "fuel" / "share_bio_gasoline.csv"
        )
        self.biodiesel = get_biofuel_share(DATA_DIR / "fuel" / "share_bio_diesel.csv")
        self.default_fuels = get_default_fuels()
        self.fuel_specs = get_fuels_specs()

    def __str__(self):
        return self.__class__.__name__

    def get_share_biofuel(self, fuel: str, country: str, years: List[int]) -> np.array:
        """
        Returns average share of biodiesel according to historical IEA stats
        with an upper limit of 30% when extrapolating.
        """

        map_array = {
            "biodiesel": self.biodiesel,
            "bioethanol": self.bioethanol,
            "biomethane": self.biomethane,
        }

        share_biofuel = np.squeeze(
            np.clip(
                map_array[fuel]
                .sel(country=country)
                .interp(
                    year=years,
                    kwargs={"fill_value": "extrapolate"},
                )
                .values,
                0,
                0.3,
            )
        )
        return share_biofuel

    def find_fuel_shares(
        self, fuel_blend: dict, fuel_type: str, country: str, years: List[int]
    ) -> [str, str, np.array, np.array]:
        """
        Find the fuel shares of the fuel blend, given a country, a fuel type and a list of years.
        """
        if fuel_type in fuel_blend:
            primary = fuel_blend[fuel_type]["primary"]["type"]

            try:
                # See if a secondary fuel type has been specified
                secondary = fuel_blend[fuel_type]["secondary fuel"]["type"]

            except KeyError:
                # A secondary fuel has not been specified, set one by default
                # Check first if the default fuel is not similar to the primary fuel

                if self.default_fuels[fuel_type]["secondary"] != primary:
                    secondary = self.default_fuels[fuel_type]["secondary"]
                else:
                    secondary = [
                        f for f in self.default_fuels[fuel_type]["all"] if f != primary
                    ][0]

            secondary_share = np.array(1) - np.zeros_like(np.array(years))

        else:

            primary = self.default_fuels[fuel_type]["primary"]
            secondary = self.default_fuels[fuel_type]["secondary"]

            if primary == "electrolysis":
                secondary_share = np.zeros_like(np.array(years))

            else:
                if fuel_type == "diesel":
                    if country in self.biodiesel.country.values:
                        secondary_share = self.get_share_biofuel(
                            "biodiesel", country, years
                        )
                    else:
                        secondary_share = self.get_share_biofuel(
                            "biodiesel", "RER", years
                        )

                elif fuel_type == "cng":
                    if country in self.biomethane.country.values:
                        secondary_share = self.get_share_biofuel(
                            "biomethane", country, years
                        )
                    else:
                        secondary_share = self.get_share_biofuel(
                            "biomethane", "RER", years
                        )
                else:
                    if country in self.bioethanol.country.values:
                        secondary_share = self.get_share_biofuel(
                            "bioethanol", country, years
                        )
                    else:
                        secondary_share = self.get_share_biofuel(
                            "bioethanol", "RER", years
                        )

        primary_share = 1 - np.squeeze(np.array(secondary_share))
        secondary_share = np.squeeze(secondary_share)

        return primary, secondary, primary_share, secondary_share

    def define_fuel_blends(
        self, powertrains: List[str], country: str, years: List[int]
    ) -> dict:
        """
        This function defines fuel blends from what is passed in `fuel_blend`.
        It populates a dictionary `self.fuel_blends` that contains the
        respective shares, lower heating values
        and CO2 emission factors of the fuels used.

        Source for LHV: https://www.bafu.admin.ch/bafu/en/home/topics/climate/state/data/climate-reporting/references.html

        :return:
        """

        fuel_blend = dict()

        fuel_to_powertrains_map = {
            "diesel": ["ICEV-d", "HEV-d", "PHEV-d", "PHEV-c-d"],
            "petrol": ["ICEV-p", "HEV-p", "PHEV-p", "PHEV-c-p"],
            "cng": ["ICEV-g"],
            "hydrogen": ["FCEV"],
        }

        _arr = lambda x: np.asarray(x) if not isinstance(x, np.ndarray) else x

        for fuel_type, pwt in fuel_to_powertrains_map.items():
            if any(i in powertrains for i in pwt):
                (
                    primary,
                    secondary,
                    primary_share,
                    secondary_share,
                ) = self.find_fuel_shares(fuel_blend, fuel_type, country, years)

                primary_share = _arr(primary_share)
                secondary_share = _arr(secondary_share)

                # add a trailing dimension to the array if it is 0D
                if primary_share.ndim == 0:
                    primary_share = np.expand_dims(primary_share, axis=0)
                if secondary_share.ndim == 0:
                    secondary_share = np.expand_dims(secondary_share, axis=0)

                fuel_blend[fuel_type] = {
                    "primary": {
                        "type": primary,
                        "share": primary_share,
                        "lhv": self.fuel_specs[primary]["lhv"],
                        "CO2": self.fuel_specs[primary]["co2"],
                        "density": self.fuel_specs[primary]["density"],
                        "name": tuple(self.fuel_specs[primary]["name"]),
                        "biogenic share": self.fuel_specs[primary]["biogenic_share"],
                    },
                    "secondary": {
                        "type": secondary,
                        "share": secondary_share,
                        "lhv": self.fuel_specs[secondary]["lhv"],
                        "CO2": self.fuel_specs[secondary]["co2"],
                        "density": self.fuel_specs[secondary]["density"],
                        "name": tuple(self.fuel_specs[secondary]["name"]),
                        "biogenic share": self.fuel_specs[secondary]["biogenic_share"],
                    },
                }

        return fuel_blend
