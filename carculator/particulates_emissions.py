"""
particulates_emissions.py contains ParticulatesEmissionsModel which calculated the amount of
abrasion particles emitted, given a driving cycle.
"""
import numpy as np
import xarray as xr


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

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "Urban delivery", "Regional delivery", "Long haul".
    :param cycle_name: name of the driving cycle. Str.


    """

    def __init__(self, cycle: np.ndarray, cycle_name: str, mass: xr.DataArray) -> None:

        self.mass = mass.values / 1000  # in tons
        self.cycle = cycle
        self.cycle_name = cycle_name
        self.cycle_environment = {
            "WLTC": {
                "urban start": 0,
                "urban stop": 590,
                "suburban start": 591,
                "suburban stop": 1023,
                "rural start": 1024,
                "rural stop": 1801,
            },
            "WLTC 3.1": {"urban start": 0, "urban stop": 590},
            "WLTC 3.2": {"suburban start": 0, "suburban stop": 433},
            "WLTC 3.3": {"rural start": 0, "rural stop": 455},
            "WLTC 3.4": {"rural start": 0, "rural stop": 323},
            "CADC Urban": {"urban start": 0, "urban stop": 994},
            "CADC Road": {"suburban start": 0, "suburban stop": 1082},
            "CADC Motorway": {"rural start": 0, "rural stop": 1068},
            "CADC Motorway 130": {"rural start": 0, "rural stop": 1068},
            "CADC": {
                "urban start": 0,
                "urban stop": 994,
                "suburban start": 995,
                "suburban stop": 2077,
                "rural start": 2078,
                "rural stop": 3144,
            },
            "NEDC": {
                "urban start": 0,
                "urban stop": 780,
                "rural start": 781,
                "rural stop": 1180,
            },
            "custom": {"urban start": 0, "urban stop": cycle[0][-1]},
        }

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

        urban_start = self.cycle_environment[self.cycle_name].get("urban start", 0)
        urban_stop = self.cycle_environment[self.cycle_name].get("urban stop", 0)
        suburban_start = self.cycle_environment[self.cycle_name].get(
            "suburban start", 0
        )
        suburban_stop = self.cycle_environment[self.cycle_name].get("suburban stop", 0)
        rural_start = self.cycle_environment[self.cycle_name].get("rural start", 0)
        rural_stop = self.cycle_environment[self.cycle_name].get("rural stop", 0)

        cycle_length = self.cycle.shape[0]
        urban_share = (urban_stop - urban_start) / cycle_length
        suburban_share = (suburban_stop - suburban_start) / cycle_length
        rural_share = (rural_stop - rural_start) / cycle_length

        tire_wear = (tire_pm10_urban + tire_pm25_urban) * urban_share
        tire_wear += (tire_pm10_rural + tire_pm25_rural) * suburban_share
        tire_wear += (tire_pm10_motorway + tire_pm25_motorway) * rural_share

        brake_wear = (brake_pm10_urban + brake_pm25_urban) * urban_share
        brake_wear += (brake_pm10_rural + brake_pm25_rural) * suburban_share
        brake_wear += (brake_pm10_motorway + brake_pm25_motorway) * rural_share

        road_wear = road_pm10 + road_pm25
        road_dust = dust_pm10 + dust_pm25

        res = np.vstack(
            (
                tire_wear[None, ...],
                brake_wear[None, ...],
                road_wear[None, ...],
                road_dust[None, ...],
            )
        )

        return res.transpose(1, 2, 0, 3, 4)

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
