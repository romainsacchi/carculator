"""
particulates_emissions.py contains ParticulatesEmissionsModel which calculated the amount of
abrasion particles emitted, given a driving cycle.
"""
from typing import Union

import numpy as np
import xarray as xr

from .hot_emissions import get_driving_cycle_compartments


def _(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, -1)
    return obj


class ParticulatesEmissionsModel:
    """
    Calculate particulates emissions based on the method described in:
    https://www.eea.europa.eu/ds_resolveuid/6USNA27I4D

    and further disaggregated in:
    https://doi.org/10.1016/j.atmosenv.2020.117886

    Include emission from:

    - brake wear
    - tire wear
    - road wear
    - re-suspended road dust

    by considering:

    - vehicle mass
    - driving situation (urban, rural, motorway)

    into the following fractions:

    - PM 2.5
    - PM 10

    Emissions are subdivided in compartments: urban, suburban and rural.

    :param velocity: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "Urban delivery", "Regional delivery", "Long haul".
    :param cycle_name: name of the driving cycle. Str.


    """

    def __init__(self, velocity: xr.DataArray, mass: xr.DataArray) -> None:

        self.mass = mass.values / 1000  # in tons
        self.velocity = velocity / 1000 * 3600  # in km/h
        self.distance = velocity.sum(dim="second") / 1000

    def get_abrasion_emissions(self) -> np.ndarray:

        (
            tire_pm10_urban,
            tire_pm10_rural,
            tire_pm10_motorway,
            tire_pm25_urban,
            tire_pm25_rural,
            tire_pm25_motorway,
        ) = self.get_tire_wear_emissions()

        (
            brake_pm10_urban,
            brake_pm10_rural,
            brake_pm10_motorway,
            brake_pm25_urban,
            brake_pm25_rural,
            brake_pm25_motorway,
        ) = self.get_brake_wear_emissions()

        road_pm10, road_pm25 = self.get_road_wear_emissions()
        dust_pm10, dust_pm25 = self.get_resuspended_road_dust()

        urban_share = self.velocity.where(self.velocity < 50, 0).sum(
            dim="second"
        ) / self.velocity.sum(dim="second")
        suburban_share = self.velocity.where(
            (self.velocity > 50) & (self.velocity <= 80), 0
        ).sum(dim="second") / self.velocity.sum(dim="second")
        rural_share = self.velocity.where(self.velocity > 80, 0).sum(
            dim="second"
        ) / self.velocity.sum(dim="second")

        tire_wear = (tire_pm10_urban + tire_pm25_urban) * urban_share.T
        tire_wear += (tire_pm10_rural + tire_pm25_rural) * suburban_share.T
        tire_wear += (tire_pm10_motorway + tire_pm25_motorway) * rural_share.T

        brake_wear = (brake_pm10_urban + brake_pm25_urban) * urban_share.T
        brake_wear += (brake_pm10_rural + brake_pm25_rural) * suburban_share.T
        brake_wear += (brake_pm10_motorway + brake_pm25_motorway) * rural_share.T

        road_wear = road_pm10 + road_pm25
        road_dust = dust_pm10 + dust_pm25

        res = np.concatenate(
            [
                _(tire_wear),
                _(brake_wear),
                _(road_wear),
                _(road_dust),
            ],
            axis=-1,
        )

        return res.transpose(0, 1, 4, 2, 3)

    def get_tire_wear_emissions(self):
        """
        Returns tire wear emissions.

        :return:
        """

        # converted to kg per vkm
        pm10_urban = (5.8 * self.mass ** (1 / 2.3)) / 1e6
        pm25_urban = (8.2 * self.mass ** (1 / 2.3)) / 1e6

        pm10_rural = (4.5 * self.mass ** (1 / 2.3)) / 1e6
        pm25_rural = (6.4 * self.mass ** (1 / 2.3)) / 1e6

        pm10_motorway = (3.8 * self.mass ** (1 / 2.3)) / 1e6
        pm25_motorway = (5.5 * self.mass ** (1 / 2.3)) / 1e6

        return (
            pm10_urban,
            pm10_rural,
            pm10_motorway,
            pm25_urban,
            pm25_rural,
            pm25_motorway,
        )

    def get_brake_wear_emissions(self):
        """
        Returns brake wear emissions.

        :return:
        """
        # converted to kg per vkm
        pm10_urban = (4.2 * np.power(self.mass, 1 / 1.9)) / 1e6
        pm25_urban = (11 * np.power(self.mass, 1 / 1.9)) / 1e6

        pm10_rural = (1.8 * np.power(self.mass, 1 / 1.5)) / 1e6
        pm25_rural = (4.5 * np.power(self.mass, 1 / 1.5)) / 1e6

        pm10_motorway = (0.4 * np.power(self.mass, 1 / 1.3)) / 1e6
        pm25_motorway = (1.0 * np.power(self.mass, 1 / 1.3)) / 1e6

        return (
            pm10_urban,
            pm10_rural,
            pm10_motorway,
            pm25_urban,
            pm25_rural,
            pm25_motorway,
        )

    def get_road_wear_emissions(self):
        """
        Returns road wear emissions.

        :return:
        """
        # converted to kg per vkm
        pm10 = (2.8 * np.power(self.mass, 1 / 1.5)) / 1e6
        pm25 = (5.1 * np.power(self.mass, 1 / 1.5)) / 1e6

        return (pm10, pm25)

    def get_resuspended_road_dust(self):
        """
        Returns re-suspended road dust emissions.

        :return:
        """
        # converted to kg per vkm
        pm10 = (2 * np.power(self.mass, 1 / 1.1)) / 1e6
        pm25 = (8.2 * np.power(self.mass, 1 / 1.1)) / 1e6

        return (pm10, pm25)
