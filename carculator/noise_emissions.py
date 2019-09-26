"""
.. module: noise_emissions.py

"""

import numpy as np
import xarray


def pn(cycle, powertrain_type):
    cycle = np.array(cycle)

    # Noise sources are calculated for speeds above 20 km/h.
    if powertrain_type in ("combustion", "electric"):
        array = np.tile((cycle - 70) / 70, 8).reshape((8, -1))
        constants = np.array((94.5, 89.2, 88, 85.9, 84.2, 86.9, 83.3, 76.1)).reshape(
            (-1, 1)
        )
        coefficients = np.array((-1.3, 7.2, 7.7, 8, 8, 8, 8, 8)).reshape((-1, 1))
        array = array * coefficients + constants

        if powertrain_type == "electric":
            # For electric cars, we add correction factors
            # We also add a 56 dB loud sound signal when the speed is below 20 km/h.
            correction = np.array((0, 1.7, 4.2, 15, 15, 15, 13.8, 0)).reshape((-1, 1))
            array -= correction
            array[:, cycle < 20] = 56
        else:
            array[:, cycle < 20] = 0
    else:
        # For non plugin-hybrids, apply electric engine noise coefficient up to 30 km/h
        # and combustion engine noise coefficients above 30 km/h
        electric = pn(cycle, "electric")
        electric_mask = cycle < 30

        array = pn(cycle, "combustion")
        array[:, electric_mask] = electric[:, electric_mask]
    return array


class NoiseEmissionsModel:
    """
    Calculate propulsion and rolling noise emissions for combustion, hybrid and electric vehicles, based on CNOSSOS model.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(self, cycle):

        self.cycle = cycle

    def rolling_noise(self):
        """Calculate noise from rolling friction.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)


        :returns: A numpy array with rolling noise (dB) for each 8 octaves, per second of driving cycle
        :rtype: numpy.array

        """
        cycle = np.array(self.cycle)
        array = np.tile(
            np.log10(cycle / 70, out=np.zeros_like(cycle), where=(cycle != 0)), 8
        ).reshape((8, -1))

        constants = np.array((79.7, 85.7, 84.5, 90.2, 97.3, 93.9, 84.1, 74.3)).reshape(
            (-1, 1)
        )
        coefficients = np.array((30, 41.5, 38.9, 25.7, 32.5, 37.2, 39, 40)).reshape(
            (-1, 1)
        )
        array = array * coefficients + constants
        array[:, cycle < 20] = 0
        return array

    def propulsion_noise(self, powertrain_type):
        """Calculate noise from propulsion engine and gearbox.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)

        :returns: A numpy array with propulsion noise (dB) for each 8 octaves, per second of driving cycle
        :rtype: numpy.array

        """
        cycle = np.array(self.cycle)

        # Noise sources are calculated for speeds above 20 km/h.
        if powertrain_type in ("combustion", "electric"):
            array = np.tile((cycle - 70) / 70, 8).reshape((8, -1))
            constants = np.array(
                (94.5, 89.2, 88, 85.9, 84.2, 86.9, 83.3, 76.1)
            ).reshape((-1, 1))
            coefficients = np.array((-1.3, 7.2, 7.7, 8, 8, 8, 8, 8)).reshape((-1, 1))
            array = array * coefficients + constants

            if powertrain_type == "electric":
                # For electric cars, we add correction factors
                # We also add a 56 dB loud sound signal when the speed is below 20 km/h.
                correction = np.array((0, 1.7, 4.2, 15, 15, 15, 13.8, 0)).reshape(
                    (-1, 1)
                )
                array -= correction
                array[:, cycle < 20] = 56
            else:
                array[:, cycle < 20] = 0
        else:
            # For non plugin-hybrids, apply electric engine noise coefficient up to 30 km/h
            # and combustion engine noise coefficients above 30 km/h
            electric = pn(cycle, "electric")
            electric_mask = cycle < 30

            array = pn(cycle, "combustion")
            array[:, electric_mask] = electric[:, electric_mask]
        return array

    def get_sound_power_per_compartment(self, powertrain_type):
        """
        Calculate sound energy (in J/s) over the driving cycle duration from sound power (in dB).
        The sound energy sums are further divided into `geographical compartments`: urban, suburban and rural.
        * *urban*: from 0 to 50 km/k
        * *suburban*: from 51 km/h to 80 km/h
        * *rural*: above 80 km/h


        :return: Sound energy (in Joules) per km driven, per geographical compartment.
        :rtype: numpy.array
        """

        # rolling noise, in dB, for each second of the driving cycle
        rolling = self.rolling_noise()
        # propulsion noise, in dB, for each second of the driving cycle
        propulsion = self.propulsion_noise(powertrain_type)

        # sum of rolling and propulsion noise sources
        total_noise = np.where(
            self.cycle != 0,
            10 * np.log10((10 ** (rolling / 10)) + (10 ** (propulsion / 10))),
            0,
        )

        # convert dBs to Watts (or J/s)
        sound_power = (10 ** -12) * (10 ** (total_noise / 10))

        distance = self.cycle.sum() / 3600

        # sum sound power over duration (J/s * s --> J) and divide by distance (--> J / km) and further
        # divide into compartments
        urban_sound = np.where(self.cycle <= 50, sound_power, 0).sum(axis=1) / distance
        suburban_sound = (
            np.where((self.cycle > 50) & (self.cycle <= 80), sound_power, 0).sum(axis=1)
            / distance
        )
        rural_sound = np.where(self.cycle > 80, sound_power, 0).sum(axis=1) / distance

        return np.array([urban_sound, suburban_sound, rural_sound])
