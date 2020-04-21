import numexpr as ne
import numpy as np
import xarray


def _(o):
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(o, (np.ndarray, xarray.DataArray)):
        return np.expand_dims(o, -1)
    else:
        return o


class HotEmissionsModel:
    """
    Calculate hot pollutants emissions based on HBEFA 4.1 data, function of speed (given by the driving cycle)
    for vehicles with a combustion engine.

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
                "rural stop": 3146,
            },
            "NEDC": {
                "urban start": 0,
                "urban stop": 780,
                "rural start": 781,
                "rural stop": 1180,
            },
        }

    def get_emissions_per_powertrain(self, powertrain_type, euro_class):
        """
        Calculate hot pollutants emissions given a powertrain type (i.e., diesel, petrol, CNG) and a EURO pollution class, per air sub-compartment
        (i.e., urban, suburban and rural).
        Note that Nh3 and N2O emissions do not depend on the speed level. FOr those, average values observed across
        different traffic situations are used instead.

        The emission sums are further divided into `air compartments`: urban, suburban and rural.

            * *urban*: from 0 to 50 km/k
            * *suburban*: from 51 km/h to 80 km/h
            * *rural*: above 80 km/h

        :param powertrain_type: "diesel", "petrol" or "CNG"
        :type powertrain_type: str
        :param euro_class: integer, corresponding to the EURO pollution class
        :type euro_class: int
        :return: Pollutants emission per km driven, per air compartment.
        :rtype: numpy.array
        """
        c = self.cycle

        em_arr = np.zeros((len(self.cycle), 11))

        if powertrain_type == "diesel":
            if euro_class == 3:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-6.771e-09 * c ** 3 + 4.125e-06 * c ** 2 - 0.0006655 * c + 0.04737"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "5.235e-08 * c ** 3 - 2.378e-06 * c ** 2 - 0.001704 * c + 0.1599"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(2.692e-06 * c ** 3 - 0.0004576 * c ** 2 + 0.02411 * c + 0.3662) + (9.421e-07 * c ** 3 - 0.0001602 * c ** 2 + 0.00844 * c + 0.1282)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate("3.662e-08 * c ** 3 - 3.196e-06 * c ** 2 - 0.0001999 * c + 0.04362")
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-6.771e-10 * c ** 3 + 4.125e-07 * c ** 2 - 6.655e-05 * c + 0.004737"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "-6.094e-09 * c ** 3 + 3.713e-06 * c ** 2 - 0.000599 * c + 0.04263"
                )
                # Pb
                # No Pb!
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "6.276e-10 * c ** 3 - 5.211e-08 * c ** 2 - 1.748e-06 * c + 0.001136"
                )
                # N2O
                em_arr[:, 8] = 0.004
                # NH3
                em_arr[:, 9] = 0.001
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "-1.131e-10 * c ** 3 + 6.889e-08 * c ** 2 - 1.111e-05 * c + 0.000791"
                )
            if euro_class == 4:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.138e-08 * c ** 3 + 3.998e-06 * c ** 2 - 0.0004721 * c + 0.02618"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.238e-08 * c ** 3 + 2.399e-06 * c ** 2 - 0.00114 * c + 0.09208"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.938e-06 * c ** 3 - 0.0002904 * c ** 2 + 0.01012 * c + 0.6372) + (8.348e-07 * c ** 3 - 0.0001251 * c ** 2 + 0.004358 * c + 0.2744)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate("-1.386e-08 * c ** 3 + 5.912e-06 * c ** 2 - 0.0006624 * c + 0.03824")
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-1.707e-09 * c ** 3 + 5.996e-07 * c ** 2 - 7.082e-05 * c + 0.003928"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "-9.675e-09 * c ** 3 + 3.398e-06 * c ** 2 - 0.0004013 * c + 0.02226"
                )
                # Pb
                # No Pb!
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "4.21e-10 * c ** 3 - 5.317e-09 * c ** 2 - 5.653e-06 * c + 0.001271"
                )
                # N2O
                em_arr[:, 8] = 0.0035
                # NH3
                em_arr[:, 9] = 0.001
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "-9.106e-11 * c ** 3 + 3.198e-08 * c ** 2 - 3.777e-06 * c + 0.0002095"
                )
            if euro_class == 5:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.098e-08 * c ** 3 + 3.824e-06 * c ** 2 - 4.702e-04 * c + 0.02442"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "9.384e-08 * c ** 3 - 2.133e-05 * c ** 2 + 0.001079 * c + 0.01645"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.634e-06 * c ** 3 - 0.0002633 * c ** 2 + 0.01063 * c + 0.5499) + (5.312e-07 * c ** 3 - 8.517e-05 * c ** 2 + 0.003398 * c + 0.1846)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate("-3.416e-09 * c ** 3 + 1.005e-06 * c ** 2 - 0.0001049 * c + 0.005359")
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-6.586e-09 * c ** 3 + 2.295e-06 * c ** 2 - 0.0002821 * c + 0.01465"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "-4.391e-09 * c ** 3 + 1.53e-06 * c ** 2 - 0.0001881 * c + 0.009768"
                )
                # Pb
                # No Pb!
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "3.508e-10 * c ** 3 + 1.216e-08 * c ** 2 - 7.032e-06 * c + 0.001306"
                )
                # N2O
                em_arr[:, 8] = 0.008
                # NH3
                em_arr[:, 9] = 0.002
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "-8.782e-11 * c ** 3 + 3.059e-08 * c ** 2 - 3.761e-06 * c + 0.0001954"
                )
            if euro_class == 6:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "1.42e-09 * c ** 3 - 3.905e-07 * c ** 2 + 6.247e-05 * c + 0.004804"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "2.6e-09 * c ** 3 + 9.256e-07 * c ** 2 - 0.0002761 * c + 0.02438"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.159e-07 * c ** 3 - 1.891e-05 * c ** 2 + 0.0008328 * c + 0.01852) + (4.056e-08 * c ** 3 - 6.618e-06 * c ** 2 + 0.0002915 * c + 0.006481)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = np.clip(
                    ne.evaluate(
                        "2.12e-09 * c ** 3 - 4.61e-07 * c ** 2 + 3.018e-05 * c - 0.0001433"
                    ),
                    0,
                    None,
                )
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "1.278e-09 * c ** 3 - 3.514e-07 * c ** 2  + 5.622e-05 * c + 0.004324"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "1.42e-10 * c ** 3 - 3.905e-08 * c ** 2 + 6.247e-06 * c + 0.0004804"
                )
                # Pb
                # No Pb!
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "2.631e-10 * c ** 3 + 1.122e-08 * c ** 2 - 6.082e-06 * c + 0.001071"
                )
                # N2O
                em_arr[:, 8] = 0.014
                # NH3
                em_arr[:, 9] = 0.006
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "1.136e-11 * c ** 3 - 3.124e-09 * c ** 2 + 4.998e-07 * c + 3.844e-05"
                )

        if powertrain_type == "petrol":
            if euro_class == 3:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "6.91e-08 * c ** 3 - 1.357e-05 * c ** 2 + 0.0009447 * c - 0.01034"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.55e-05 * c ** 3 - 0.003246 * c ** 2 + 0.2153 * c - 4.007"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.383e-07 * c ** 3 - 1.557e-05 * c ** 2 - 0.0004461 * c + 0.1016) + (6.916e-09 * c ** 3 - 7.783e-07 * c ** 2 - 2.23e-05 * c + 0.005081)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.73e-08 * c ** 3 - 2.854e-06 * c ** 2 + 0.0001556 * c - 0.001652"
                )
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.073e-08 * c ** 3 - 4.072e-06 * c ** 2 + 0.0002834 * c - 0.003101"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "4.837e-08 * c ** 3 - 9.502e-06 * c ** 2 + 0.0006613 * c - 0.007235"
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.606e-12 * c ** 3 + 1.045e-08 * c ** 2 - 1.348e-06 * c + 0.0001014"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-1.521e-10 * c ** 3 + 1.655e-07 * c ** 2 - 2.133e-05 * c + 0.001606"
                )
                # N2O
                em_arr[:, 8] = 2.5e-3
                # NH3
                em_arr[:, 9] = 0.05
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "8.935e-09 * c ** 3 - 1.755e-06 * c ** 2 + 0.0001222 * c - 0.001336"
                )
            if euro_class == 4:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "8.594e-08 * c ** 3 - 1.688e-05 * c ** 2 + 0.001022 * c - 0.01428"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.519e-05 * c ** 3 - 0.003153 * c ** 2 + 0.2025 * c - 3.702"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.427e-07 * c ** 3 - 2.47e-05 * c ** 2 + 0.001137 * c + 0.03733) + (7.135e-09 * c ** 3 - 1.235e-06 * c ** 2 + 5.683e-05 * c + 0.001866)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "5.603e-09 * c ** 3 - 6.869e-07 * c ** 2 + 3.006e-05 * c - 2.313e-05"
                )
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "3.437e-08 * c ** 3 - 6.753e-06 * c ** 2 + 0.0004088 * c - 0.005713"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "5.156e-08 * c ** 3 - 1.013e-05 * c ** 2 + 0.0006132 * c - 0.008569"
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.227e-11 * c ** 3 + 1.322e-08 * c ** 2 - 1.55e-06 * c + 0.0001023"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-3.525e-10 * c ** 3 + 2.092e-07 * c ** 2 - 2.453e-05 * c + 0.001619"
                )
                # N2O
                em_arr[:, 8] = 1.0e-3
                # NH3
                em_arr[:, 9] = 0.04
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "8.594e-10 * c ** 3 - 1.688e-07 * c ** 2 + 1.022e-05 * c - 0.0001428"
                )
            if euro_class == 5:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-6.38e-08 * c ** 3 - 1.259e-05 * c ** 2 + 0.0007827 * c - 0.009119"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "9.941e-06 * c ** 3 - 0.00207 * c ** 2 + 0.137 * c - 2.51"
                )
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(2.053e-08 * c ** 3 - 2.671e-06 * c ** 2 - 0.0002311 * c + 0.04269) + (1.026e-09 * c ** 3 - 1.336e-07 * c ** 2 - 1.156e-05 * c + 0.002134)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.273e-08 * c ** 3 - 2.491e-06 * c ** 2 + 0.0001522 * c - 0.001577"
                )
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.552e-08 * c ** 3 - 5.037e-06 * c ** 2 + 0.0003131 * c - 0.003647"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "3.828e-08 * c ** 3 - 7.756e-06 * c ** 2 + 0.0004696 * c - 0.005471"
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.028e-11 * c ** 3 + 1.193e-08 * c ** 2 - 1.411e-06 * c + 9.467e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-3.21e-10 * c ** 3 + 1.889e-07 * c ** 2 - 2.234e-05 * c + 0.001499"
                )
                # N2O
                em_arr[:, 8] = 8.5e-4
                # NH3
                em_arr[:, 9] = 0.019
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "-6.38e-10 * c ** 3 + 1.259e-07 * c ** 2 - 7.827e-06 * c + 9.119e-05"
                )
            if euro_class == 6:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-6.698e-09 * c ** 3 + 2.661e-06 * c ** 2 - 0.0002239 * c + 0.008068"
                )
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "5.385e-07 * c ** 3 - 0.0001104 * c ** 2 + 0.006821 * c + 0.1707"
                )
                # NO + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(1.379e-07 * c ** 3 - 2.708e-05 * c ** 2 + 0.001678 * c - 0.01762) + (6.894e-09 * c ** 3 - 1.354e-06 * c ** 2 + 8.388e-05 * c - 0.0008809)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-3.18e-10 * c ** 3 + 2.95e-07 * c ** 2 - 2.911e-05 * c + 0.001359"
                )
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-2.679e-09 * c ** 3 + 1.064e-06 * c ** 2 - 8.957e-05 * c + 0.003227"
                )
                # NMVOC
                em_arr[:, 5] = ne.evaluate(
                    "-4.019e-09 * c ** 3 + 1.597e-06 * c ** 2 - 0.0001344 * c + 0.004841"
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.062e-11 * c ** 3 + 9.669e-09 * c ** 2 - 1.099e-06 * c + 7.273e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-3.264e-10 * c ** 3 + 1.53e-07 * c ** 2 - 1.74e-05 * c + 0.001151"
                )
                # N2O
                em_arr[:, 8] = 8.5e-4
                # NH3
                em_arr[:, 9] = 0.019
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "-6.698e-11 * c ** 3 + 2.661e-08 * c ** 2 - 2.239e-06 * c + 8.068e-05"
                )

        if powertrain_type == "CNG":
            if euro_class == 3:
                # HC
                em_arr[:, 0] = np.clip(
                    ne.evaluate(
                        "6.91e-08 * c ** 3 - 1.357e-05 * c ** 2 + 0.0009447 * c - 0.01034"
                    ),
                    0,
                    None,
                )
                # CO
                em_arr[:, 1] = np.clip(
                    ne.evaluate(
                        "1.472e-05 * c ** 3 - 0.003084 * c ** 2 + 0.2045 * c - 3.806"
                    ),
                    0,
                    None,
                )
                # NO + NO2
                em_arr[:, 2] = np.clip(
                    ne.evaluate(
                    "(2.075e-07 * c ** 3 - 2.335e-05 * c ** 2 - 0.0006691 * c + 0.1524) + (1.452e-08 * c ** 3 - 1.634e-06 * c ** 2 - 4.684e-05 * c + 0.01067)"
                ), 0, None)

                # PM <= 2.5 um
                em_arr[:, 3] = np.clip(
                    ne.evaluate(
                    "1.73e-08 * c ** 3 - 2.854e-06 * c ** 2 + 0.0001556 * c - 0.001652"
                ), 0, None)
                # CH4
                em_arr[:, 4] = np.clip(
                    ne.evaluate(
                        "4.215e-08 * c ** 3 - 8.28e-06 * c ** 2 + 0.0005763 * c - 0.006305"
                    ),
                    0,
                    None,
                )
                # NMVOC
                em_arr[:, 5] = np.clip(
                    ne.evaluate(
                        "2.695e-08 * c ** 3 - 5.294e-06 * c ** 2 + 0.0003684 * c - 0.004031"
                    ),
                    0,
                    None,
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-4.803e-12 * c ** 3 + 5.227e-09 * c ** 2 - 6.738e-07 * c + 5.072e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-7.603e-11 * c ** 3 + 8.274e-08 * c ** 2 - 1.067e-05 * c + 0.0008028"
                )
                # N2O
                em_arr[:, 8] = 2.5e-4
                # NH3
                em_arr[:, 9] = 0.03
                # Benzene
                em_arr[:, 10] = np.clip(
                    ne.evaluate(
                    "4.468e-09 * c ** 3 - 8.776e-07 * c ** 2 + 6.108e-05 * c - 0.0006682"
                ), 0, None)
            if euro_class == 4:
                # HC
                em_arr[:, 0] = np.clip(
                    ne.evaluate(
                        "6.557e-09 * c ** 3 + 3.168e-06 * c ** 2 - 0.0004758 * c + 0.03329"
                    ),
                    0,
                    None,
                )
                # CO
                em_arr[:, 1] = np.clip(
                    ne.evaluate(
                        "6.843e-06 * c ** 3 - 0.001373 * c ** 2 + 0.08654 * c - 1.442"
                    ),
                    0,
                    None,
                )
                # NO + NO2
                em_arr[:, 2] = np.clip(
                    ne.evaluate(
                    "(1.775e-06 * c ** 3 - 0.0003727 * c ** 2 + 0.02381 * c - 0.3833) + (1.398e-07 * c ** 3 - 2.943e-05 * c ** 2 + 0.001887 * c - 0.03122)"
                ), 0, None)

                # PM <= 2.5 um
                em_arr[:, 3] = np.clip(
                    ne.evaluate(
                    "3.087e-08 * c ** 3 - 5.349e-06 * c ** 2 + 0.0002932 * c - 0.002924"
                ), 0, None)
                # CH4
                em_arr[:, 4] = np.clip(
                    ne.evaluate(
                        "-1.646e-08 * c ** 3 + 7.339e-06 * c ** 2 - 0.000706 * c + 0.03439"
                    ),
                    0,
                    None,
                )
                # NMVOC
                em_arr[:, 5] = np.clip(
                    ne.evaluate(
                        "2.302e-08 * c ** 3 - 4.17e-06 * c ** 2 + 0.0002302 * c - 0.001097"
                    ),
                    0,
                    None,
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.752e-12 * c ** 3 + 5.948e-09 * c ** 2 - 7.025e-07 * c + 4.688e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-1.544e-10 * c ** 3 + 9.415e-08 * c ** 2 - 1.112e-05 * c + 0.0007421"
                )
                # N2O
                em_arr[:, 8] = 2.5e-4
                # NH3
                em_arr[:, 9] = 0.03
                # Benzene
                em_arr[:, 10] = np.clip(
                    ne.evaluate(
                    "4.326e-10 * c ** 3 - 8.507e-08 * c ** 2 + 5.158e-06 * c - 7.231e-05"
                ), 0, None)
            if euro_class == 5:
                # HC
                em_arr[:, 0] = np.clip(
                    ne.evaluate(
                        "-1.029e-08 * c ** 3 + 8.099e-06 * c ** 2 - 0.000892 * c + 0.04839"
                    ),
                    0,
                    None,
                )
                # CO
                em_arr[:, 1] = np.clip(
                    ne.evaluate(
                        "4.389e-06 * c ** 3 - 0.0008576 * c ** 2 + 0.0548 * c - 0.8267"
                    ),
                    0,
                    None,
                )
                # NO + NO2
                em_arr[:, 2] = np.clip(
                    ne.evaluate(
                    "(2.835e-07 * c ** 3 - 5.778e-05 * c ** 2 + 0.003404 * c - 0.02752) + (2.238e-08 * c ** 3 - 4.585e-06 * c ** 2 + 0.000276 * c - 0.002851)"
                ), 0, None)

                # PM <= 2.5 um
                em_arr[:, 3] = np.clip(
                    ne.evaluate(
                    "3.568e-08 * c ** 3 - 6.356e-06 * c ** 2 + 0.000347 * c - 0.002842"
                ), 0, None)
                # CH4
                em_arr[:, 4] = np.clip(
                    ne.evaluate(
                        "-2.619e-08 * c ** 3 + 1.075e-05 * c ** 2 - 0.001026 * c + 0.0469"
                    ),
                    0,
                    None,
                )
                # NMVOC
                em_arr[:, 5] = np.clip(
                    ne.evaluate(
                        "1.59e-08 * c ** 3 - 2.653e-06 * c ** 2 + 0.000134 * c + 0.001487"
                    ),
                    0,
                    None,
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.722e-12 * c ** 3 + 5.618e-09 * c ** 2 - 6.622e-07 * c + 4.414e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-1.539e-10 * c ** 3 + 8.892e-08 * c ** 2 - 1.048e-05 * c + 0.0006987"
                )
                # N2O
                em_arr[:, 8] = 8.5e-4
                # NH3
                em_arr[:, 9] = 0.019
                # Benzene
                em_arr[:, 10] = np.clip(
                    ne.evaluate(
                    "3.216e-10 * c ** 3 - 6.349e-08 * c ** 2 + 3.948e-06 * c - 4.585e-05"
                ), 0, None)
            if euro_class == 6:
                # HC
                em_arr[:, 0] = np.clip(
                    ne.evaluate(
                        "5.502e-08 * c ** 3 - 1.239e-05 * c ** 2 + 0.000846 * c - 0.004044"
                    ),
                    0,
                    None,
                )
                # CO
                em_arr[:, 1] = np.clip(
                    ne.evaluate(
                        "2.678e-06 * c ** 3 - 0.0005585 * c ** 2 + 0.03632 * c - 0.5159"
                    ),
                    0,
                    None,
                )
                # NO + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(2.163e-08 * c ** 3 - 3.758e-06 * c ** 2 + 5.485e-05 * c + 0.02848) + (4.35e-09 * c ** 3 - 7.625e-07 * c ** 2 + 2.251e-05 * c + 0.002413)"
                )
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.997e-09 * c ** 3 - 2.116e-07 * c ** 2 - 1.277e-05 * c + 0.002135"
                )
                # CH4
                em_arr[:, 4] = np.clip(
                    ne.evaluate(
                        "4.568e-08 * c ** 3 - 1.052e-05 * c ** 2 + 0.0007284 * c - 0.003946"
                    ),
                    0,
                    None,
                )
                # NMVOC
                em_arr[:, 5] = np.clip(
                    ne.evaluate(
                        "9.338e-09 * c ** 3 - 1.872e-06 * c ** 2 + 0.0001175 * c - 9.78e-05"
                    ),
                    0,
                    None,
                )
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-1.38e-11 * c ** 3 + 6.464e-09 * c ** 2 - 7.366e-07 * c + 4.665e-05"
                )
                # SO2
                em_arr[:, 7] = ne.evaluate(
                    "-2.185e-10 * c ** 3 + 1.023e-07 * c ** 2 - 1.166e-05 * c + 0.0007384"
                )
                # N2O
                em_arr[:, 8] = 8.5e-4
                # NH3
                em_arr[:, 9] = 0.019
                # Benzene
                em_arr[:, 10] = ne.evaluate(
                    "9.493e-11 * c ** 3 - 1.694e-08 * c ** 2 + 9.587e-07 * c + 4.34e-06"
                )

        if powertrain_type not in ("petrol", "diesel", "CNG"):
            raise TypeError("The powertrain type is not valid.")

        # If the driving cycle selected is one of the driving cycles for which carculator has specifications,
        # we use the driving cycle "official" road section types to compartmentalize emissions.
        # If the driving cycle selected is instead specified by the user (passed directly as an array), we used
        # speed levels to compartmentalize emissions.

        if self.cycle_name in self.cycle_environment:
            distance = self.cycle.sum() / 3600

            if "urban start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["urban start"]
                stop = self.cycle_environment[self.cycle_name]["urban stop"]
                dist_urban = self.cycle[start:stop].sum() / 3600
                urban = np.mean(em_arr[start:stop, :], axis=0) * (dist_urban / distance)
                urban /= 1000  # going from grams to kg

            else:
                urban = np.zeros((11))

            if "suburban start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["suburban start"]
                stop = self.cycle_environment[self.cycle_name]["suburban stop"]
                dist_suburban = self.cycle[start:stop].sum() / 3600
                suburban = np.mean(em_arr[start:stop, :], axis=0) * (
                    dist_suburban / distance
                )
                suburban /= 1000  # going from grams to kg

            else:
                suburban = np.zeros((11))

            if "rural start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["rural start"]
                stop = self.cycle_environment[self.cycle_name]["rural stop"]
                dist_rural = self.cycle[start:stop].sum() / 3600
                rural = np.mean(em_arr[start:stop, :], axis=0) * (dist_rural / distance)
                rural /= 1000  # going from grams to kg

            else:
                rural = np.zeros((11))

        else:
            distance = self.cycle.sum() / 3600
            dist_urban = self.cycle[self.cycle <= 50].sum() / 3600
            dist_suburban = (
                self.cycle[(self.cycle > 50) & (self.cycle <= 80)].sum() / 3600
            )
            dist_rural = self.cycle[self.cycle > 80].sum() / 3600

            urban = np.mean(em_arr[self.cycle <= 50], axis=0) * (dist_urban / distance)
            urban /= 1000  # going from grams to kg
            suburban = np.mean(
                em_arr[(self.cycle > 50) & (self.cycle <= 80)], axis=0
            ) * (dist_suburban / distance)
            suburban /= 1000  # going from grams to kg
            rural = np.mean(em_arr[self.cycle > 80], axis=0) * (dist_rural / distance)
            rural /= 1000  # going from grams to kg

            suburban = np.nan_to_num(suburban)
            rural = np.nan_to_num(rural)

        return (
            np.hstack(
                (urban.reshape(11, -1), suburban.reshape(11, -1), rural.reshape(11, -1))
            )
        ).reshape(1, 33, 1, 1)
