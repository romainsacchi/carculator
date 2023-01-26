"""
noise_emissions.py contains NoiseEmissionsModel which calculates noise emissions, in joules,
given a driving cycle and a powertrain type.
"""

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr

MAP_PWT = {
    "ICEV-p": "ICEV",
    "ICEV-d": "ICEV",
    "ICEV-g": "ICEV",
    "HEV-p": "ICEV",
    "HEV-d": "ICEV",
    "PHEV-c-p": "ICEV",
    "PHEV-c-d": "ICEV",
    "PHEV-d": "ICEV",
    "PHEV-p": "ICEV",
    "BEV": "BEV",
    "BEV-depot": "BEV",
    "BEV-opp": "BEV",
    "BEV-motion": "BEV",
    "FCEV": "BEV",
    "PHEV-e": "BEV",
    "Human": "BEV",
}

MAP_SIZES = {
    "3.5t": "medium",
    "7.5t": "medium",
    "18t": "medium",
    "26t": "heavy",
    "32t": "heavy",
    "40t": "heavy",
    "60t": "heavy",
    "9m": "medium",
    "13m-city": "heavy",
    "13m-city-double": "heavy",
    "18m": "heavy",
    "13m-coach": "heavy",
    "13m-coach-double": "heavy",
    "Kick-scooter": "bicycle",
    "Bicycle <25": "bicycle",
    "Bicycle <45": "bicycle",
    "Bicycle cargo": "bicycle",
    "Moped <4kW": "scooter",
    "Scooter <4kW": "scooter",
    "Scooter 4-11kW": "scooter",
    "Motorcycle 11-35kW": "motorcycle",
    "Motorcycle 4-11kW": "motorcycle",
    "Motorcycle >35kW": "motorcycle",
}


def _(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, -1)
    return obj


def get_noise_coefficients(filepath: Path) -> Union[None, xr.DataArray]:
    """Noise coefficients extracted for vehicles from CNOSSOS-EU 2018
    detailed by size, powertrain and EURO class for each octave.

    :param filepath: Path to the noise coefficients file.
    :type filepath: Path
    :return: Noise coefficients.
    :rtype: xr.DataArray
    """

    try:
        df = pd.read_csv(filepath, sep=",")
        cols = ["octave", "coefficient"]

        if "powertrain" in df.columns:
            cols.append("powertrain")

        if "size" in df.columns:
            cols.append("size")

        ef = df.groupby(cols)["value"].mean().to_xarray()

    except FileNotFoundError:
        ef = None

    return ef


class NoiseEmissionsModel:
    """
    Calculate propulsion and rolling noise emissions for combustion, hybrid and electric vehicles,
    based on CNOSSOS model.

    :param velocity: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type velocity: pandas.Series

    """

    def __init__(self, velocity: xr.DataArray, vehicle_type: str) -> None:

        self.velocity = velocity / 1000 * 3600  # km/h to m/s
        self.rolling_coefficients = get_noise_coefficients(
            Path(__file__).parent
            / "data"
            / "emission_factors"
            / vehicle_type
            / "rolling_noise_coefficients.csv"
        )
        self.propulsion_coefficients = get_noise_coefficients(
            Path(__file__).parent
            / "data"
            / "emission_factors"
            / vehicle_type
            / "propulsion_noise_coefficients.csv"
        )

    def rolling_noise(self) -> np.ndarray:
        """Calculate noise from rolling friction.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)

        :returns: A numpy array with rolling noise (dB)
        for each 8 octaves, per second of driving cycle
        :rtype: numpy.array

        """

        if self.rolling_coefficients is None:
            return np.zeros_like(_(self.velocity))

        _nz = lambda x: np.where(x < 1, 1, x)

        array = np.repeat(
            np.log10(_nz(_(self.velocity) / 70), where=(_(self.velocity) > 0)),
            8,
            axis=-1,
        )

        if "size" in self.rolling_coefficients.dims:
            coefficients = self.rolling_coefficients.sel(
                coefficient="a",
                size=[MAP_SIZES[s] for s in self.velocity.coords["size"].values],
            )
            constants = self.rolling_coefficients.sel(
                coefficient="b",
                size=[MAP_SIZES[s] for s in self.velocity.coords["size"].values],
            )

        else:
            coefficients = self.rolling_coefficients.sel(coefficient="a")
            constants = self.rolling_coefficients.sel(coefficient="b")

        array = array * coefficients.T.values + constants.T.values

        return array

    def propulsion_noise(self) -> np.ndarray:
        """Calculate noise from propulsion engine and gearbox.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)

        For electric cars, special coefficients are applied from
        (`Pallas et al. 2016 <https://www.sciencedirect.com/science/article/pii/S0003682X16301608>`_ )

        Also, for electric cars, a warning signal of 56 dB is added when the car drives at 20 km/h or lower.

        :returns: A numpy array with propulsion noise (dB) for all 8 octaves, per second of driving cycle
        :rtype: numpy.array

        """

        if self.propulsion_coefficients is None:
            return np.zeros_like(_(self.velocity))

        _nz = lambda x: np.where(x < 1, 1, x)

        array = np.repeat(
            np.log10(_nz((_(self.velocity) - 70) / 70), where=(_(self.velocity) > 0)),
            8,
            axis=-1,
        )

        coefficients = self.propulsion_coefficients.sel(coefficient="a")
        constants = self.propulsion_coefficients.sel(coefficient="b")

        if "size" in coefficients.dims:

            a = coefficients.sel(
                powertrain=[MAP_PWT[p] for p in self.velocity.powertrain.values],
                size=[MAP_SIZES[s] for s in self.velocity.coords["size"].values],
            )
            b = constants.sel(
                powertrain=[MAP_PWT[p] for p in self.velocity.powertrain.values],
                size=[MAP_SIZES[s] for s in self.velocity.coords["size"].values],
            )

        else:

            a = coefficients.sel(
                powertrain=[MAP_PWT[p] for p in self.velocity.powertrain.values],
            )
            b = constants.sel(
                powertrain=[MAP_PWT[p] for p in self.velocity.powertrain.values],
            )

        correction = np.where(a.powertrain == "BEV", 1, 0)

        correction = _(correction) * np.array((0, 1.7, 4.2, 15, 15, 15, 13.8, 0))

        if "size" in coefficients.dims:
            a = a.values.transpose(1, 2, 0)
            b = b.values.transpose(1, 2, 0)
        else:
            a = a.T.values[:, None, :]
            b = b.T.values[:, None, :]

        array = array * a + b

        array = array - correction[:, None, :]

        return array

    def get_sound_power_per_compartment(self) -> np.ndarray:
        """
        Calculate sound energy (in J/s) over the driving cycle duration from
        sound power (in dB). The sound energy sums are further divided into
        `geographical compartments`: urban, suburban and rural,
        based on the profile of the driving cycle.


        :return: Sound energy (in Joules) per km driven, per geographical compartment.
        :rtype: numpy.ndarray
        """

        # rolling noise, in dB, for each second of the driving cycle
        rolling = self.rolling_noise()
        # propulsion noise, in dB, for each second of the driving cycle
        propulsion = self.propulsion_noise()

        # sum of rolling and propulsion noise sources

        total_noise = np.where(
            _(self.velocity) > 0,
            np.log10((10 ** (rolling / 10))) + np.log10((10 ** (propulsion / 10))),
            0,
        )

        # convert dBs to Watts (or J/s)
        sound_power = (10**-12) * (10 ** (total_noise / 10))

        # If the driving cycle selected is one of the driving cycle for
        # which carculator_utils has specifications,
        # we use the driving cycle "official" road section types to
        # compartmentalize emissions.
        # If the driving cycle selected is instead specified by the
        # user (passed directly as an array), we used
        # speed levels to compartmentalize emissions.

        distance = (self.velocity / 3600).sum(axis=0)

        urban_noise = np.where(_(self.velocity) <= 50, sound_power, 0).sum(axis=0) / _(
            distance
        )

        suburban_noise = np.where(
            (_(self.velocity) > 50) & (_(self.velocity) <= 80), sound_power, 0
        ).sum(axis=0) / _(distance)

        rural_noise = np.where(_(self.velocity) > 80, sound_power, 0).sum(axis=0) / _(
            distance
        )

        res = np.concatenate(
            (
                urban_noise,
                suburban_noise,
                rural_noise,
            ),
            axis=-1,
        )

        return res.transpose(3, 2, -1, 1, 0)
