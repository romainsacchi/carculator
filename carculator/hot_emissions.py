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
        Also includes cold start emissions. Cold start emissions are also given by HBEFA 4.1 and expressed in given in g/start.
        Cold start emissions are divided by the average trip length in Europe (20 km), to normalize them per km.

        The emission sums are further divided into `air compartments`: urban, suburban and rural.

            * *urban*: from 0 to 50 km/k
            * *suburban*: from 51 km/h to 80 km/h
            * *rural*: above 80 km/h

        :param powertrain_type: "diesel", "petrol" or "CNG"
        :type powertrain_type: str
        :param euro_class: integer, corresponding to the EURO pollution class
        :type euro_class: float
        :return: Pollutants emission per km driven, per air compartment.
        :rtype: numpy.array
        """
        c = self.cycle

        # average trip length in Europe https://core.ac.uk/download/pdf/82726264.pdf
        # to normalize cold start emissions per km
        #(start_per_day * 365 / 12000) = 20

        # We also add NMHC: hot emissions + running losses + soak + diurnal emissions emissions from evaporation during parking time, driving and resting for gasoline cars.

        # average number of start per day per vehicle, according to section 3.2.6.2.2 of
        # https://www.bafu.admin.ch/dam/bafu/de/dokumente/luft/fachinfo-daten/Switzerlands-Informative-Inventory-Report-2019.pdf.download.pdf/switzerlands-informative-rep-2020.pdf
        
        start_per_day = 2.35
        stops_per_day = 2.35

        

        em_arr = np.zeros((len(self.cycle), 10))

        if powertrain_type == "diesel":
            if euro_class == 0:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-2.776e-07 * c ** 3 + 0.0001241 * c ** 2 - 0.01622 * c + 0.7348"
                ) + (0.107 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "8.11e-07 * c ** 3 - 3.256e-05 * c ** 2 - 0.01559 * c + 1.565"
                ) + (0.54 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.194e-06 * c ** 3 - 0.0002084 * c ** 2 + 0.01159 * c + 0.459)"
                ) + (0.698 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-5.921e-08 * c ** 3 + 6.236e-05 * c ** 2 - 0.007949 * c + 0.3523"
                ) + (0.112 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-6.663e-09 * c ** 3 + 2.978e-06 * c ** 2 - 0.0003892 * c + 0.01764"
                ) + (0.0026 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-2.71e-07 * c ** 3 + 0.0001211 * c ** 2 - 0.01583 * c + 0.7172"
                ) + (0.105 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0
                # NH3
                em_arr[:, 8] = 0.001 + (0.001 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-4.636e-09 * c ** 3 + 2.072e-06 * c ** 2 - 0.0002708 * c + 0.01227"
                ) + (0.0018 * (start_per_day * 365 / 12000))

            if euro_class == 1:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-8.425e-08 * c ** 3 + 2.798e-05 * c ** 2 - 0.003277 * c + 0.1743"
                ) + (0.07 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "2.836e-07 * c ** 3 - 2.236e-05 * c ** 2 - 0.004932 * c + 0.6634"
                ) + (0.38 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.212e-06 * c ** 3 - 0.0001786 * c ** 2 + 0.007101 * c + 0.4974) "
                ) + (0.738 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-6.042e-08 * c ** 3 + 3.566e-05 * c ** 2 - 0.004057 * c + 0.2134"
                ) + (0.122 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-2.022e-09 * c ** 3 + 6.715e-07 * c ** 2 - 7.866e-05 * c + 0.004182"
                ) + (0.0017 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-8.223e-08 * c ** 3 + 2.731e-05 * c ** 2 - 0.003199 * c + 0.1701"
                ) + (0.0678 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.004 + (0.0035 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.001 + (0.001 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-1.407e-09 * c ** 3 + 4.672e-07 * c ** 2 - 5.473e-05 * c + 0.00291"
                ) + (0.00116 * (start_per_day * 365 / 12000))

            if euro_class == 2:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-5.856e-08 * c ** 3 + 1.935e-05 * c ** 2 - 0.002251 * c + 0.1117"
                ) + (0.04 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "-1.096e-07 * c ** 3 + 4.647e-05 * c ** 2 - 0.00699 * c + 0.4191"
                ) + (0.159 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.293e-06 * c ** 3 - 0.0001745 * c ** 2 + 0.00481 * c + 0.6594)"
                ) + (1.014 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-1.135e-07 * c ** 3 + 3.772e-05 * c ** 2 - 0.003626 * c + 0.1598"
                ) + (0.06 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-2.928e-09 * c ** 3 + 9.675e-07 * c ** 2 - 0.0001126 * c + 0.005583"
                ) + (0.002 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-5.563e-08 * c ** 3 + 1.838e-05 * c ** 2 - 0.002139 * c + 0.1061"
                ) + (0.03778 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.006 + (0.0055 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.001 + (0.001 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-9.779e-10 * c ** 3 + 3.232e-07 * c ** 2 - 3.76e-05 * c + 0.001865"
                ) + (0.06 * (start_per_day * 365 / 12000))
            if euro_class == 3:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-6.771e-09 * c ** 3 + 4.125e-06 * c ** 2 - 0.0006655 * c + 0.04737"
                ) + (0.0237 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "5.235e-08 * c ** 3 - 2.378e-06 * c ** 2 - 0.001704 * c + 0.1599"
                ) + (0.069 * (start_per_day * 365 / 12000))
                # NOx + NO2
                em_arr[:, 2] = ne.evaluate(
                    "(2.692e-06 * c ** 3 - 0.0004576 * c ** 2 + 0.02411 * c + 0.3662)"
                ) + (1.25 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "3.662e-08 * c ** 3 - 3.196e-06 * c ** 2 - 0.0001999 * c + 0.04362"
                ) + (0.028 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-6.771e-10 * c ** 3 + 4.125e-07 * c ** 2 - 6.655e-05 * c + 0.004737"
                ) + (0.0024 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-6.094e-09 * c ** 3 + 3.713e-06 * c ** 2 - 0.000599 * c + 0.04263"
                ) + (0.0213 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.004 + (0.005092 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.001 + (0.001 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-1.131e-10 * c ** 3 + 6.889e-08 * c ** 2 - 1.111e-05 * c + 0.000791"
                ) + (0.000396 * (start_per_day * 365 / 12000))
            if euro_class == 4:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.138e-08 * c ** 3 + 3.998e-06 * c ** 2 - 0.0004721 * c + 0.02618"
                ) + (0.011575 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.238e-08 * c ** 3 + 2.399e-06 * c ** 2 - 0.00114 * c + 0.09208"
                ) + (0.0426 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.938e-06 * c ** 3 - 0.0002904 * c ** 2 + 0.01012 * c + 0.6372)"
                ) + (0.936 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-1.386e-08 * c ** 3 + 5.912e-06 * c ** 2 - 0.0006624 * c + 0.03824"
                ) + (0.002197 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-1.707e-09 * c ** 3 + 5.996e-07 * c ** 2 - 7.082e-05 * c + 0.003928"
                ) + (0.001736 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-9.675e-09 * c ** 3 + 3.398e-06 * c ** 2 - 0.0004013 * c + 0.02226"
                ) + (0.0098 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.0035 + (0.005 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.001 + (0.001 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-9.106e-11 * c ** 3 + 3.198e-08 * c ** 2 - 3.777e-06 * c + 0.0002095"
                ) + (9.26e-5 * (start_per_day * 365 / 12000))
            if euro_class == 5:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.098e-08 * c ** 3 + 3.824e-06 * c ** 2 - 4.702e-04 * c + 0.02442"
                ) + (0.008327 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "9.384e-08 * c ** 3 - 2.133e-05 * c ** 2 + 0.001079 * c + 0.01645"
                ) + (0.031848 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.634e-06 * c ** 3 - 0.0002633 * c ** 2 + 0.01063 * c + 0.5499) "
                ) + (1.112 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-3.416e-09 * c ** 3 + 1.005e-06 * c ** 2 - 0.0001049 * c + 0.005359"
                ) + (0.0023 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-6.586e-09 * c ** 3 + 2.295e-06 * c ** 2 - 0.0002821 * c + 0.01465"
                ) + (0.005 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-4.391e-09 * c ** 3 + 1.53e-06 * c ** 2 - 0.0001881 * c + 0.009768"
                )
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.008 + (0.0074 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.002 + (1.44e-4 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-8.782e-11 * c ** 3 + 3.059e-08 * c ** 2 - 3.761e-06 * c + 0.0001954"
                ) + (6.66e-5 * (start_per_day * 365 / 12000))

            if euro_class == 6.0:
                # EURO 6-ab
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "3.156e-08 * c ** 3 - 6.092e-06 * c ** 2 + 0.0002999 * c + 0.0125"
                ) + (0.0178 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "8.375e-08 * c ** 3 - 1.282e-05 * c ** 2 - 0.0001148 * c + 0.1052"
                ) + (0.0799 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.183e-06 * c ** 3 - 0.0002066 * c ** 2 + 0.01106 * c + 0.2807)"
                ) + (0.595 * (start_per_day * 365 / 12000))
                # PM
                em_arr[:, 3] = ne.evaluate(
                    "-2.807e-10 * c ** 3 + 2.184e-07 * c ** 2 - 2.521e-05 * c + 0.00126"
                ) + (7.27e-4 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.84e-08 * c ** 3 - 5.482e-06 * c ** 2 + 0.0002699 * c + 0.01125"
                ) + (0.016 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "3.156e-09 * c ** 3 - 6.092e-07 * c ** 2 + 2.999e-05 * c + 0.00125"
                ) + (0.0018 * (start_per_day * 365 / 12000))

                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.01 + (0.012 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.006 + (0.007 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "2.525e-10 * c ** 3 - 4.873e-08 * c ** 2 + 2.399e-06 * c + 0.0001"
                ) + (1.43e-4 * (start_per_day * 365 / 12000))

            if euro_class == 6.1:
                # EURO 6-c
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "9.641e-09 * c ** 3 - 1.899e-06 * c ** 2 + 0.000122 * c + 0.008276"
                ) + (0.0135 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "3.007e-08 * c ** 3 - 3.641e-06 * c ** 2 - 0.0002374 * c + 0.05297"
                ) + (0.0372 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(5.086e-07 * c ** 3 - 8.662e-05 * c ** 2 + 0.004389 * c + 0.1175)"
                ) + (0.2376 * (start_per_day * 365 / 12000))
                # PM
                em_arr[:, 3] = ne.evaluate(
                    "1.302e-09 * c ** 3 - 2.298e-07 * c ** 2 + 1.141e-05 * c + 0.0003679"
                ) + (6.01e-4 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "8.677e-09 * c ** 3 - 1.709e-06 * c ** 2 + 0.0001098 * c + 0.007448"
                ) + (0.01217 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "9.641e-10 * c ** 3 - 1.899e-07 * c ** 2 + 1.22e-05 * c + 0.0008276"
                ) + (0.00135 * (start_per_day * 365 / 12000))
                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.01 + (0.01173 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.006 + (0.0069 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "7.713e-11 * c ** 3 - 1.519e-08 * c ** 2 + 9.757e-07 * c + 6.621e-05"
                ) + (1.08e-4 * (start_per_day * 365 / 12000))

            if euro_class == 6.2:
                # EURO 6-d-temp
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "1.876e-10 * c **3 - 1.283e-07 * c ** 2 + 5.19e-05 * c + 0.005754"
                ) + (0.0111 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "3.432e-09 * c **3 + 9.136e-07 * c ** 2  - 0.0002983 * c + 0.02698"
                ) + (0.016 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.335e-07 * c **3 - 2.156e-05 * c ** 2  + 0.0009388 * c + 0.02214)"
                ) + (0.043 * (start_per_day * 365 / 12000))

                # PM
                em_arr[:, 3] = ne.evaluate(
                    "2.216e-09 * c **3 - 4.813e-07 * c ** 2  + 3.152e-05 * c - 0.0001162"
                ) + (5.3e-4 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "1.688e-10 * c **3 - 1.155e-07 * c ** 2  + 4.671e-05 * c + 0.005178"
                ) + (0.01 * (start_per_day * 365 / 12000))

                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "1.876e-11 * c **3 - 1.283e-08 * c ** 2  + 5.19e-06 * c + 0.0005754"
                ) + (0.00112 * (start_per_day * 365 / 12000))

                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.01 + (0.01173 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.006 + (0.0069 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "1.501e-12 * c **3 - 1.026e-09 * c ** 2  + 4.152e-07 * c + 4.603e-05"
                ) + (8.93e-5 * (start_per_day * 365 / 12000))

            if euro_class == 6.3:
                # EURO 6-d
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "1.42e-09 * c ** 3 - 3.905e-07 * c ** 2 + 6.247e-05 * c + 0.004804"
                ) + (0.0098 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "2.6e-09 * c ** 3 + 9.256e-07 * c ** 2 - 0.0002761 * c + 0.02438"
                ) + (0.0142 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.159e-07 * c ** 3 - 1.891e-05 * c ** 2 + 0.0008328 * c + 0.01852)"
                ) + (0.0365 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "2.12e-09 * c ** 3 - 4.61e-07 * c ** 2 + 3.018e-05 * c - 0.0001433"
                ) + (4.7e-4 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "1.278e-09 * c ** 3 - 3.514e-07 * c ** 2  + 5.622e-05 * c + 0.004324"
                ) + (0.0088 * (start_per_day * 365 / 12000))

                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "1.42e-10 * c ** 3 - 3.905e-08 * c ** 2 + 6.247e-06 * c + 0.0004804"
                ) + (9.83e-4 * (start_per_day * 365 / 12000))

                # Pb
                # No Pb!
                # N2O
                em_arr[:, 7] = 0.014 + (0.01173 * (start_per_day * 365 / 12000))

                # NH3
                em_arr[:, 8] = 0.006 + (0.0069 * (start_per_day * 365 / 12000))

                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "1.136e-11 * c ** 3 - 3.124e-09 * c ** 2 + 4.998e-07 * c + 3.844e-05"
                ) + (7.86e-5 * (start_per_day * 365 / 12000))

        if powertrain_type == "petrol":
            if euro_class == 0:

                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-2.397e-06 * c ** 3 + 0.0007389 * c ** 2 - 0.07681 * c + 3.276"
                ) + (1.218 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "-4.753e-07 * c ** 3 + 0.001648 * c ** 2 - 0.269 * c + 16.26"
                ) + (8.176 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(-2.973e-07 * c ** 3 + 0.00018 * c ** 2 - 0.01601 * c + 1.184)"
                ) + (1.156 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "6.057e-08 * c ** 3 - 1.191e-05 * c ** 2 + 0.0007664 * c - 0.01338"
                ) + (0.0065 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-7.567e-08 * c ** 3 + 2.369e-05 * c ** 2 - 0.002499 * c + 0.1095"
                ) + (0.0414 * (start_per_day * 365 / 12000))

                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(-2.321e-06 * c ** 3 + 0.0007152 * c ** 2 - 0.07431 * c + 3.166) + ((1.177 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((4.855452538 * 365) / 12000)")

                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "7.335e-12 * c ** 3 + 5.838e-09 * c ** 2 - 8.995e-07 * c + 8.89e-05"
                ) + (7.2e-5 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 0.08 + (0.0073 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.04 + (0.002 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-9.501e-08 * c ** 3 + 2.991e-05 * c ** 2 - 0.003171 * c + 0.1403"
                ) + (0.0534 * (start_per_day * 365 / 12000))

            if euro_class == 1:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "2.844e-07 * c ** 3 - 6.305e-05 * c ** 2 + 0.004306 * c - 0.003236"
                ) + (0.1145 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "8.223e-06 * c ** 3 - 0.001798 * c ** 2 + 0.1286 * c - 2.099"
                ) + (1.199 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(-5.15e-08 * c ** 3 + 8.582e-05 * c ** 2 - 0.01024 * c + 0.5778)"
                ) + (0.377 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "3.23e-08 * c ** 3 - 5.563e-06 * c ** 2 + 0.0003178 * c - 0.004209"
                ) + (0.0049 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.389e-08 * c ** 3 - 5.296e-06 * c ** 2 + 0.0003617 * c - 0.0002718"
                ) + (0.0096 * (start_per_day * 365 / 12000))

                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(2.605e-07 * c ** 3 - 5.775e-05 * c ** 2 + 0.003944 * c - 0.002964) + ((0.105 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.117 * 365) / 12000)")

                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "3.472e-13 * c ** 3 + 8.382e-09 * c ** 2 - 1.168e-06 * c + 9.978e-05"
                ) + (7.46e-5 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 0.01 + (0.0099 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.12 + (0.0955 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "3.677e-08 * c ** 3 - 8.152e-06 * c ** 2 + 0.0005568 * c - 0.0004184"
                ) + (0.0148 * (start_per_day * 365 / 12000))

            if euro_class == 2:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-2.562e-08 * c ** 3 + 7.201e-06 * c ** 2 - 0.0004912 * c + 0.03658"
                ) + (0.0423 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.656e-05 * c ** 3 - 0.003461 * c ** 2 + 0.2291 * c - 4.258"
                ) + (0.89 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(4.472e-07 * c ** 3 - 5.03e-05 * c ** 2 - 0.001425 * c + 0.3272)"
                ) + (0.202 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "5.509e-08 * c ** 3 - 9.59e-06 * c ** 2 + 0.0005602 * c - 0.008264"
                ) + (0.0079 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-5.123e-09 * c ** 3 + 1.44e-06 * c ** 2 - 9.824e-05 * c + 0.007316"
                ) + (0.0085 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(-2.049e-08 * c ** 3 + 5.76e-06 * c ** 2 - 0.000393 * c + 0.02926) + ((0.034 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.06 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-8.2e-12 * c ** 3 + 1.092e-08 * c ** 2 - 1.428e-06 * c + 0.0001095"
                ) + (7.7e-5 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 0.004 + (0.0049 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.12 + (0.123 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-3.312e-09 * c ** 3 + 9.31e-07 * c ** 2 - 6.351e-05 * c + 0.00473"
                ) + (0.00547 * (start_per_day * 365 / 12000))

            if euro_class == 3:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "6.91e-08 * c ** 3 - 1.357e-05 * c ** 2 + 0.0009447 * c - 0.01034"
                ) + (0.01838 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.55e-05 * c ** 3 - 0.003246 * c ** 2 + 0.2153 * c - 4.007"
                ) + (0.8266 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.383e-07 * c ** 3 - 1.557e-05 * c ** 2 - 0.0004461 * c + 0.1016)"
                ) + (0.06258 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.73e-08 * c ** 3 - 2.854e-06 * c ** 2 + 0.0001556 * c - 0.001652"
                ) + (0.0031 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.073e-08 * c ** 3 - 4.072e-06 * c ** 2 + 0.0002834 * c - 0.003101"
                ) + (0.0055 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(4.837e-08 * c ** 3 - 9.502e-06 * c ** 2 + 0.0006613 * c - 0.007235) + ((0.01286 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.13 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.606e-12 * c ** 3 + 1.045e-08 * c ** 2 - 1.348e-06 * c + 0.0001014"
                ) + (7.04e-5 * (start_per_day * 365 / 12000))
                # N2O
                em_arr[:, 7] = 2.5e-3 + (4.5e-5 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.05 + (0.037 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "8.935e-09 * c ** 3 - 1.755e-06 * c ** 2 + 0.0001222 * c - 0.001336"
                ) + (0.00237 * (start_per_day * 365 / 12000))
            if euro_class == 4:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "8.594e-08 * c ** 3 - 1.688e-05 * c ** 2 + 0.001022 * c - 0.01428"
                ) + (0.0071 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.519e-05 * c ** 3 - 0.003153 * c ** 2 + 0.2025 * c - 3.702"
                ) + (0.4911 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.427e-07 * c ** 3 - 2.47e-05 * c ** 2 + 0.001137 * c + 0.03733)"
                ) + (0.059 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "5.603e-09 * c ** 3 - 6.869e-07 * c ** 2 + 3.006e-05 * c - 2.313e-05"
                ) + (0.001679 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "3.437e-08 * c ** 3 - 6.753e-06 * c ** 2 + 0.0004088 * c - 0.005713"
                ) + (0.00285 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(5.156e-08 * c ** 3 - 1.013e-05 * c ** 2 + 0.0006132 * c - 0.008569) + ((0.00427 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.074 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.227e-11 * c ** 3 + 1.322e-08 * c ** 2 - 1.55e-06 * c + 0.0001023"
                ) + (6.58e-5 * (start_per_day * 365 / 12000))
                # N2O
                em_arr[:, 7] = 1.0e-3 + (5.72e-4 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.04 + (0.037 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "8.594e-10 * c ** 3 - 1.688e-07 * c ** 2 + 1.022e-05 * c - 0.0001428"
                ) + (7.12e-5 * (start_per_day * 365 / 12000))

            if euro_class == 5:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-6.38e-08 * c ** 3 - 1.259e-05 * c ** 2 + 0.0007827 * c - 0.009119"
                ) + (0.008958 * (start_per_day * 365 / 12000))

                # CO
                em_arr[:, 1] = ne.evaluate(
                    "9.941e-06 * c ** 3 - 0.00207 * c ** 2 + 0.137 * c - 2.51"
                ) + (0.608 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(2.053e-08 * c ** 3 - 2.671e-06 * c ** 2 - 0.0002311 * c + 0.04269)"
                ) + (0.02396 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.273e-08 * c ** 3 - 2.491e-06 * c ** 2 + 0.0001522 * c - 0.001577"
                ) + (0.00179 * (start_per_day * 365 / 12000))

                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.552e-08 * c ** 3 - 5.037e-06 * c ** 2 + 0.0003131 * c - 0.003647"
                ) + (0.00358 * (start_per_day * 365 / 12000))

                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(3.828e-08 * c ** 3 - 7.756e-06 * c ** 2 + 0.0004696 * c - 0.005471) + ((0.00537 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.038 * 365) / 12000)")

                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.028e-11 * c ** 3 + 1.193e-08 * c ** 2 - 1.411e-06 * c + 9.467e-05"
                ) + (5.98e-5 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 8.5e-4 + (5.72e-4 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (0.016 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-6.38e-10 * c ** 3 + 1.259e-07 * c ** 2 - 7.827e-06 * c + 9.119e-05"
                ) + (8.95e-5 * (start_per_day * 365 / 12000))

            if euro_class == 6.0:
                # EURO 6-ab
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "1.861e-08 * c ** 3 - 3.288e-06 * c ** 2 + 0.0001836 * c + 0.001005"
                ) + (6.04e-03 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "5.22e-06 * c ** 3 - 0.001098 * c ** 2 + 0.07276 * c - 1.225"
                ) + (3.90e-01 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(3.746e-09 * c ** 3 - 3.01e-07 * c ** 2 - 0.0001707 * c + 0.04678)"
                ) + (3.90e-02 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "5.402e-10 * c ** 3 + 7.947e-08 * c ** 2 - 2.936e-05 * c + 0.002263"
                ) + (1.09e-03 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "7.443e-09 * c ** 3 - 1.315e-06 * c ** 2 + 7.344e-05 * c + 0.0004021"
                ) + (2.42e-03 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(1.116e-08 * c ** 3 - 1.973e-06 * c ** 2 + 0.0001102 * c + 0.0006032) + ((3.62e-03 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((1.0 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.745e-11 * c ** 3 + 1.268e-08 * c ** 2 - 1.436e-06 * c + 9.009e-05"
                ) + (5.33e-05 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 2.5e-4 + (5.73e-04 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (1.88e-02 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "1.861e-10 * c ** 3 - 3.288e-08 * c ** 2 + 1.836e-06 * c + 1.005e-05"
                ) + (6.04e-05 * (start_per_day * 365 / 12000))

            if euro_class == 6.1:
                # EURO 6-c
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "7.335e-09 * c ** 3 - 1.343e-06 * c ** 2 + 8.808e-05 * c - 0.0008577"
                ) + (1.98e-03 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.848e-06 * c ** 3 - 0.0003671 * c ** 2 + 0.02352 * c - 0.2435"
                ) + (3.17e-01 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(-7.317e-09 * c ** 3 + 2.929e-06 * c ** 2 - 0.0003728 * c + 0.02627)"
                ) + (1.52e-02 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-3.9e-10 * c ** 3 + 1.517e-07 * c ** 2 - 1.904e-05 * c + 0.001032"
                ) + (4.11e-04 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "2.934e-09 * c ** 3 - 5.373e-07 * c ** 2 + 3.523e-05 * c - 0.0003431"
                ) + (7.91e-04 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(4.401e-09 * c ** 3 - 8.06e-07 * c ** 2 + 5.285e-05 * c - 0.0005146) + ((1.19e-03 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((.956 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.561e-11 * c ** 3 + 1.161e-08 * c ** 2 - 1.309e-06 * c + 8.137e-05"
                ) + (4.77e-05 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 2.5e-4 + (5.73e-04 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (1.88e-02 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "7.335e-11 * c ** 3 - 1.343e-08 * c ** 2 + 8.808e-07 * c - 8.577e-06"
                ) + (1.98e-05 * (start_per_day * 365 / 12000))

            if euro_class == 6.2:
                # EURO 6-d-temp
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-7.731e-09 * c ** 3 + 2.908e-06 * c ** 2 - 0.0002415 * c + 0.008536"
                ) + (5.42e-03 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "9.503e-07 * c ** 3 - 0.0001961 * c ** 2 + 0.0125 * c + 0.06203"
                ) + (3.57e-01 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.334e-07 * c ** 3 - 2.541e-05 * c ** 2 + 0.001513 * c - 0.01215)"
                ) + (2.28e-02 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-9.717e-10 * c ** 3 + 4.452e-07 * c ** 2 - 4.106e-05 * c + 0.001644"
                ) + (1.07e-03 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-3.093e-09 * c ** 3 + 1.163e-06 * c ** 2 - 9.659e-05 * c + 0.003414"
                ) + (2.17e-03 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(-4.639e-09 * c ** 3 + 1.745e-06 * c ** 2 - 0.0001449 * c + 0.005121) + ((3.25e-03 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((.93 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-2.405e-11 * c ** 3 + 1.101e-08 * c ** 2 - 1.246e-06 * c + 7.88e-05"
                ) + (4.67e-05 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 2.5e-4 + (5.73e-04 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.017 + (1.88e-02 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-7.731e-11 * c ** 3 + 2.908e-08 * c ** 2 - 2.415e-06 * c + 8.536e-05"
                ) + (5.42e-05 * (start_per_day * 365 / 12000))

            if euro_class == 6.3:
                # EURO 6-d
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.092e-08 * c ** 3 + 3.665e-06 * c ** 2 - 0.0002973 * c + 0.009515"
                ) + (5.05e-03 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.742e-07 * c ** 3 - 3.012e-05 * c ** 2 + 0.001088 * c + 0.2879"
                ) + (3.16e-01 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.437e-07 * c ** 3 - 2.75e-05 * c ** 2 + 0.001649 * c - 0.01736)"
                ) + (2.03e-02 * (start_per_day * 365 / 12000))
                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "-7.457e-10 * c ** 3 + 4.05e-07 * c ** 2 - 3.638e-05 * c + 0.001481"
                ) + (1.15e-03 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-4.369e-09 * c ** 3 + 1.466e-06 * c ** 2 - 0.0001189 * c + 0.003806"
                ) + (2.02e-03 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "(-6.554e-09 * c ** 3 + 2.199e-06 * c ** 2 - 0.0001784 * c + 0.005709) + ((3.03e-03 * start_per_day * 365) / 12000) + 0.002332 + ((0.04411241 * stops_per_day * 365) / 12000) + ((.91 * 365) / 12000)")
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-1.843e-11 * c ** 3 + 8.969e-09 * c ** 2 - 1.026e-06 * c + 6.87e-05"
                ) + (4.26e-05 * (start_per_day * 365 / 12000))

                # N2O
                em_arr[:, 7] = 2.5e-4 + (5.73e-04 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (1.88e-02 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-1.092e-10 * c ** 3 + 3.665e-08 * c ** 2 - 2.973e-06 * c + 9.515e-05"
                ) + (5.05e-05 * (start_per_day * 365 / 12000))

        if powertrain_type == "CNG":

            if euro_class == 2:
                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-4.483e-08 * c ** 3 + 1.26e-05 * c ** 2 - 0.0008596 * c + 0.06402"
                ) + (0.105707355 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.573e-05 * c ** 3 - 0.003288 * c ** 2 + 0.2176 * c - 4.045"
                ) + (0.801136017 * (start_per_day * 365 / 12000))
                # NO
                em_arr[:, 2] = ne.evaluate(
                    "(3.802e-07 * c ** 3 - 4.275e-05 * c ** 2 - 0.001211 * c + 0.2782)"
                ) + (0.141743928 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "5.509e-08 * c ** 3 - 9.59e-06 * c ** 2 + 0.0005602 * c - 0.008264"
                ) + (0.007912168 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-3.202e-08 * c ** 3 + 9.001e-06 * c ** 2 - 0.000614 * c + 0.04573"
                ) + (0.097250767 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "-1.281e-08 * c ** 3 + 3.6e-06 * c ** 2 - 0.0002456 * c + 0.01829"
                ) + (0.008456588 * (start_per_day * 365 / 12000))
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-3.679e-12 * c ** 3 + 4.9e-09 * c ** 2 - 6.404e-07 * c + 4.914e-05"
                )

                # N2O
                em_arr[:, 7] = 0.002 + (0.000572766 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.075 + (0.037029848 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "-1.656e-09 * c ** 3 + 4.655e-07 * c ** 2 - 3.176e-05 * c + 0.002365"
                )

            if euro_class == 3:

                # HC
                em_arr[:, 0] = ne.evaluate(
                    "6.91e-08 * c ** 3 - 1.357e-05 * c ** 2 + 0.0009447 * c - 0.01034"
                ) + (0.018378355 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "1.472e-05 * c ** 3 - 0.003084 * c ** 2 + 0.2045 * c - 3.806"
                ) + (0.74397397 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(2.075e-07 * c ** 3 - 2.335e-05 * c ** 2 - 0.0006691 * c + 0.1524)"
                ) + (0.125161722 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.73e-08 * c ** 3 - 2.854e-06 * c ** 2 + 0.0001556 * c - 0.001652"
                ) + (0.003107077 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "4.215e-08 * c ** 3 - 8.28e-06 * c ** 2 + 0.0005763 * c - 0.006305"
                ) + (0.016908085 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "2.695e-08 * c ** 3 - 5.294e-06 * c ** 2 + 0.0003684 * c - 0.004031"
                ) + (0.001470268 * (start_per_day * 365 / 12000))
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-4.803e-12 * c ** 3 + 5.227e-09 * c ** 2 - 6.738e-07 * c + 5.072e-05"
                )

                # N2O
                em_arr[:, 7] = 2.5e-4 + (0.000572766 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.03 + (0.037029848 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "4.468e-09 * c ** 3 - 8.776e-07 * c ** 2 + 6.108e-05 * c - 0.0006682"
                )
            if euro_class == 4:

                # HC
                em_arr[:, 0] = ne.evaluate(
                    "6.557e-09 * c ** 3 + 3.168e-06 * c ** 2 - 0.0004758 * c + 0.03329"
                ) + (0.048561502 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "6.843e-06 * c ** 3 - 0.001373 * c ** 2 + 0.08654 * c - 1.442"
                ) + (0.354609609 * (start_per_day * 365 / 12000))

                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(1.775e-06 * c ** 3 - 0.0003727 * c ** 2 + 0.02381 * c - 0.3833)"
                ) + (0.136211395 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "3.087e-08 * c ** 3 - 5.349e-06 * c ** 2 + 0.0002932 * c - 0.002924"
                ) + (0.007032036 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-1.646e-08 * c ** 3 + 7.339e-06 * c ** 2 - 0.000706 * c + 0.03439"
                ) + (0.044676583 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "2.302e-08 * c ** 3 - 4.17e-06 * c ** 2 + 0.0002302 * c - 0.001097"
                ) + (0.00388492 * (start_per_day * 365 / 12000))
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.752e-12 * c ** 3 + 5.948e-09 * c ** 2 - 7.025e-07 * c + 4.688e-05"
                )

                # N2O
                em_arr[:, 7] = 2.5e-4 + (0.000572766 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.03 + (0.037029848 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "4.326e-10 * c ** 3 - 8.507e-08 * c ** 2 + 5.158e-06 * c - 7.231e-05"
                )

            if euro_class == 5:

                # HC
                em_arr[:, 0] = ne.evaluate(
                    "-1.029e-08 * c ** 3 + 8.099e-06 * c ** 2 - 0.000892 * c + 0.04839"
                ) + (0.06067254 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "4.389e-06 * c ** 3 - 0.0008576 * c ** 2 + 0.0548 * c - 0.8267"
                ) + (0.431372195 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(2.835e-07 * c ** 3 - 5.778e-05 * c ** 2 + 0.003404 * c - 0.02752)"
                ) + (0.042887595 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "3.568e-08 * c ** 3 - 6.356e-06 * c ** 2 + 0.000347 * c - 0.002842"
                ) + (0.008190296 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "-2.619e-08 * c ** 3 + 1.075e-05 * c ** 2 - 0.001026 * c + 0.0469"
                ) + (0.055818737 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "1.59e-08 * c ** 3 - 2.653e-06 * c ** 2 + 0.000134 * c + 0.001487"
                ) + (0.004853803 * (start_per_day * 365 / 12000))
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-9.722e-12 * c ** 3 + 5.618e-09 * c ** 2 - 6.622e-07 * c + 4.414e-05"
                )

                # N2O
                em_arr[:, 7] = 8.5e-4 + (0.000572766 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (0.016007606 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "3.216e-10 * c ** 3 - 6.349e-08 * c ** 2 + 3.948e-06 * c - 4.585e-05"
                )

            if euro_class == 6:

                # HC
                em_arr[:, 0] = ne.evaluate(
                    "5.502e-08 * c ** 3 - 1.239e-05 * c ** 2 + 0.000846 * c - 0.004044"
                ) + (2.17e-02 * (start_per_day * 365 / 12000))
                # CO
                em_arr[:, 1] = ne.evaluate(
                    "2.678e-06 * c ** 3 - 0.0005585 * c ** 2 + 0.03632 * c - 0.5159"
                ) + (1.59e-01 * (start_per_day * 365 / 12000))
                # NOx
                em_arr[:, 2] = ne.evaluate(
                    "(2.163e-08 * c ** 3 - 3.758e-06 * c ** 2 + 5.485e-05 * c + 0.02848)"
                ) + (1.07e-02 * (start_per_day * 365 / 12000))

                # PM <= 2.5 um
                em_arr[:, 3] = ne.evaluate(
                    "1.997e-09 * c ** 3 - 2.116e-07 * c ** 2 - 1.277e-05 * c + 0.002135"
                ) + (5.05e-05 * (start_per_day * 365 / 12000))
                # CH4
                em_arr[:, 4] = ne.evaluate(
                    "4.568e-08 * c ** 3 - 1.052e-05 * c ** 2 + 0.0007284 * c - 0.003946"
                ) + (1.55e-03 * (start_per_day * 365 / 12000))
                # NMHC: hot emissions + running losses + soak + diurnal emissions
                em_arr[:, 5] = ne.evaluate(
                    "9.338e-09 * c ** 3 - 1.872e-06 * c ** 2 + 0.0001175 * c - 9.78e-05"
                ) + (1.74e-03 * (start_per_day * 365 / 12000))
                # Pb
                em_arr[:, 6] = ne.evaluate(
                    "-1.38e-11 * c ** 3 + 6.464e-09 * c ** 2 - 7.366e-07 * c + 4.665e-05"
                )

                # N2O
                em_arr[:, 7] = 8.5e-4 + (5.73e-04 * (start_per_day * 365 / 12000))
                # NH3
                em_arr[:, 8] = 0.019 + (1.88e-02 * (start_per_day * 365 / 12000))
                # Benzene
                em_arr[:, 9] = ne.evaluate(
                    "9.493e-11 * c ** 3 - 1.694e-08 * c ** 2 + 9.587e-07 * c + 4.34e-06"
                )

        if powertrain_type not in ("petrol", "diesel", "CNG"):
            raise TypeError("The powertrain type is not valid.")

        # In case the fit produces negative numbers
        em_arr[em_arr < 0] = 0

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
                urban = np.zeros((10))

            if "suburban start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["suburban start"]
                stop = self.cycle_environment[self.cycle_name]["suburban stop"]
                dist_suburban = self.cycle[start:stop].sum() / 3600
                suburban = np.mean(em_arr[start:stop, :], axis=0) * (
                    dist_suburban / distance
                )
                suburban /= 1000  # going from grams to kg

            else:
                suburban = np.zeros((10))

            if "rural start" in self.cycle_environment[self.cycle_name]:
                start = self.cycle_environment[self.cycle_name]["rural start"]
                stop = self.cycle_environment[self.cycle_name]["rural stop"]
                dist_rural = self.cycle[start:stop].sum() / 3600
                rural = np.mean(em_arr[start:stop, :], axis=0) * (dist_rural / distance)
                rural /= 1000  # going from grams to kg

            else:
                rural = np.zeros((10))

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
                (urban.reshape(10, -1), suburban.reshape(10, -1), rural.reshape(10, -1))
            )
        ).reshape(1, 30, 1, 1)
