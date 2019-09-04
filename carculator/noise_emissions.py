"""
.. module: noise_emissions.py

"""

import numpy as np
import xarray

def _(o):
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(o, (np.ndarray, xarray.DataArray)):
        return np.expand_dims(o, -1)
    else:
        return o


class NoiseEmissionsModel():
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

        rolling = np.zeros((8, self.cycle.shape[0]))
        rolling[0] = np.where(self.cycle >= 20, 79.7 + 30 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                    where=(self.cycle != 0)), 0)
        rolling[1] = np.where(self.cycle >= 20, 85.7 + 41.5 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                      where=(self.cycle != 0)), 0)
        rolling[2] = np.where(self.cycle >= 20, 84.5 + 38.9 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                      where=(self.cycle != 0)), 0)
        rolling[3] = np.where(self.cycle >= 20, 90.2 + 25.7 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                      where=(self.cycle != 0)), 0)
        rolling[4] = np.where(self.cycle >= 20, 97.3 + 32.5 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                      where=(self.cycle != 0)), 0)
        rolling[5] = np.where(self.cycle >= 20, 93.9 + 37.2 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                      where=(self.cycle != 0)), 0)
        rolling[6] = np.where(self.cycle >= 20, 84.1 + 39 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                    where=(self.cycle != 0)), 0)
        rolling[7] = np.where(self.cycle >= 20, 74.3 + 40 * np.log10(self.cycle / 70, out=np.zeros_like(self.cycle),
                                                                    where=(self.cycle != 0)), 0)

        return rolling

    def propulsion_noise(self, powertrain_type):
        """Calculate noise from propulsion engine and gearbox.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)

        :returns: A numpy array with propulsion noise (dB) for each 8 octaves, per second of driving cycle
        :rtype: numpy.array

        """

        propulsion = np.zeros((8, self.cycle.shape[0]))

        # Noise sources are calculated for speeds above 20 km/h.
        if powertrain_type == 'combustion':

            propulsion[0] = np.where(self.cycle >= 20, 94.5 - 1.3 * ((self.cycle - 70) / 70), 0)
            propulsion[1] = np.where(self.cycle >= 20, 89.2 + 7.2 * ((self.cycle - 70) / 70), 0)
            propulsion[2] = np.where(self.cycle >= 20, 88 + 7.7 * ((self.cycle - 70) / 70), 0)
            propulsion[3] = np.where(self.cycle >= 20, 85.9 + 8 * ((self.cycle - 70) / 70), 0)
            propulsion[4] = np.where(self.cycle >= 20, 84.2 + 8 * ((self.cycle - 70) / 70), 0)
            propulsion[5] = np.where(self.cycle >= 20, 86.9 + 8 * ((self.cycle - 70) / 70), 0)
            propulsion[6] = np.where(self.cycle >= 20, 83.3 + 8 * ((self.cycle - 70) / 70), 0)
            propulsion[7] = np.where(self.cycle >= 20, 76.1 + 8 * ((self.cycle - 70) / 70), 0)

        elif powertrain_type == 'electric':
            # For electric cars, we add correction factors
            # We also add a 56 dB loud sound signal when the speed is below 20 km/h.

            propulsion[0] = np.where(self.cycle >= 20, 94.5 - 1.3 * ((self.cycle - 70) / 70), 56)
            propulsion[1] = np.where(self.cycle >= 20, (89.2 + 7.2 * ((self.cycle - 70) / 70)) - 1.7, 56)
            propulsion[2] = np.where(self.cycle >= 20, (88 + 7.7 * ((self.cycle - 70) / 70))- 4.2, 56)
            propulsion[3] = np.where(self.cycle >= 20, (85.9 + 8 * ((self.cycle - 70) / 70))- 15, 56)
            propulsion[4] = np.where(self.cycle >= 20, (84.2 + 8 * ((self.cycle - 70) / 70))- 15, 56)
            propulsion[5] = np.where(self.cycle >= 20, (86.9 + 8 * ((self.cycle - 70) / 70))- 15, 56)
            propulsion[6] = np.where(self.cycle >= 20, (83.3 + 8 * ((self.cycle - 70) / 70))- 13.8, 56)
            propulsion[7] = np.where(self.cycle >= 20, 76.1 + 8 * ((self.cycle - 70) / 70), 56)

        else:

            # For non plugin-hybrids, apply electric engine noise coefficient up to 30 km/h
            # and combustion engine noise coefficients above 30 km/h
            propulsion[0] = np.where(self.cycle >= 30, 94.5 - 1.3 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), 94.5 - 1.3 * ((self.cycle - 70) / 70), 56))
            propulsion[1] = np.where(self.cycle >= 30, 89.2 + 7.2 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (89.2 + 7.2 * ((self.cycle - 70) / 70)) - 1.7, 56))
            propulsion[2] = np.where(self.cycle >= 30, 88 + 7.7 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (88 + 7.7 * ((self.cycle - 70) / 70)) - 4.2, 56))
            propulsion[3] = np.where(self.cycle >= 30, 85.9 + 8 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (85.9 + 8 * ((self.cycle - 70) / 70)) - 15, 56))
            propulsion[4] = np.where(self.cycle >= 30, 84.2 + 8 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (84.2 + 8 * ((self.cycle - 70) / 70)) - 15, 56))
            propulsion[5] = np.where(self.cycle >= 30, 86.9 + 8 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (86.9 + 8 * ((self.cycle - 70) / 70)) - 15, 56))
            propulsion[6] = np.where(self.cycle >= 30, 83.3 + 8 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), (83.3 + 8 * ((self.cycle - 70) / 70)) - 13.8, 56))
            propulsion[7] = np.where(self.cycle >= 30, 76.1 + 8 * ((self.cycle - 70) / 70),
                                     np.where((self.cycle >= 20)&(self.cycle < 30), 76.1 + 8 * ((self.cycle - 70) / 70), 56))

        return propulsion

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
        total_noise = np.where(self.cycle != 0, 10 * np.log10((10 ** (rolling / 10)) + (10 ** (propulsion / 10))), 0)

        # convert dBs to Watts (or J/s)
        sound_power = (10 ** -12) * (10 ** (total_noise / 10))

        distance = self.cycle.sum() / 3600

        # sum sound power over duration (J/s * s --> J) and divide by distance (--> J / km) and further
        # divide into compartments
        urban_sound = np.where(self.cycle <= 50, sound_power, 0).sum(axis=1) / distance
        suburban_sound = np.where((self.cycle > 50) & (self.cycle <= 80), sound_power, 0).sum(axis=1) / distance
        rural_sound = np.where(self.cycle > 80, sound_power, 0).sum(axis=1) / distance

        return np.array([urban_sound, suburban_sound, rural_sound])

