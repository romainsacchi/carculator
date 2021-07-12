import numexpr as ne
import numpy as np


class NoiseEmissionsModel:
    """
    Calculate propulsion and rolling noise emissions for combustion, hybrid and electric vehicles, based on CNOSSOS model.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(self, cycle, cycle_name):

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
        }

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

        return array

    def propulsion_noise(self, powertrain_type):
        """Calculate noise from propulsion engine and gearbox.
        Model from CNOSSOS-EU project
        (http://publications.jrc.ec.europa.eu/repository/bitstream/JRC72550/cnossos-eu%20jrc%20reference%20report_final_on%20line%20version_10%20august%202012.pdf)

        For electric cars, special coefficients are applied from
        (`Pallas et al. 2016 <https://www.sciencedirect.com/science/article/pii/S0003682X16301608>`_ )

        Also, for electric cars, a warning signal of 56 dB is added when the car drives at 20 km/h or lower.

        :returns: A numpy array with propulsion noise (dB) for all 8 octaves, per second of driving cycle
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

                # Warming signal for electric cars of 56 dB at 20 km/h or lower
                array[:, cycle < 20] = 56

        else:
            # For non plugin-hybrids, apply electric engine noise coefficient up to 30 km/h
            # and combustion engine noise coefficients above 30 km/h
            electric = self.propulsion_noise("electric")
            electric_mask = cycle < 30

            array = self.propulsion_noise("combustion")
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

        if powertrain_type not in ("combustion", "electric", "hybrid"):
            raise TypeError("The powertrain type is not valid.")

        # rolling noise, in dB, for each second of the driving cycle
        rolling = self.rolling_noise()
        # propulsion noise, in dB, for each second of the driving cycle
        propulsion = self.propulsion_noise(powertrain_type)

        # sum of rolling and propulsion noise sources
        c = self.cycle

        total_noise = ne.evaluate(
            "where(c != 0, 10 * log10((10 ** (rolling / 10)) + (10 ** (propulsion / 10))), 0)"
        )

        # convert dBs to Watts (or J/s)
        sound_power = ne.evaluate("(10 ** -12) * (10 ** (total_noise / 10))")

        # If the driving cycle selected is one of the driving cycles for which carculator has specifications,
        # we use the driving cycle "official" road section types to compartmentalize emissions.
        # If the driving cycle selected is instead specified by the user (passed directly as an array), we used
        # speed levels to compartmentalize emissions.

        if self.cycle_name in self.cycle_environment:
            distance = self.cycle.sum() / 3600

            if "urban start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["urban start"]
                stop = self.cycle_environment[self.cycle_name]["urban stop"]
                urban = np.sum(sound_power[:, start:stop], axis=1) / distance

            else:
                urban = np.zeros(8)

            if "suburban start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["suburban start"]
                stop = self.cycle_environment[self.cycle_name]["suburban stop"]
                suburban = np.sum(sound_power[:, start:stop], axis=1) / distance

            else:
                suburban = np.zeros(8)

            if "rural start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["rural start"]
                stop = self.cycle_environment[self.cycle_name]["rural stop"]
                rural = np.sum(sound_power[:, start:stop], axis=1) / distance

            else:
                rural = np.zeros(8)

        else:
            distance = self.cycle.sum() / 3600

            # sum sound power over duration (J/s * s --> J) and divide by distance (--> J / km) and further
            # divide into compartments
            urban = ne.evaluate("sum(where(c <= 50, sound_power, 0), 1)") / distance
            suburban = (
                ne.evaluate("sum(where((c > 50) & (c <= 80), sound_power, 0), 1)")
                / distance
            )
            rural = ne.evaluate("sum(where(c > 80, sound_power, 0), 1)") / distance

        return np.array([urban, suburban, rural])
