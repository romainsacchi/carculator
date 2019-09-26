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
    Calculate hot pollutants emissions based on HBEFA 3.3 data, function of speed (given by the driving cycle)
    for vehicles with a combustion engine.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(self, cycle):

        self.cycle = cycle

    def get_emissions_per_powertrain(self, powertrain_type):
        """
        Calculate hot pollutants emissions given a powertrain type (i.e., diesel, petrol, CNG), per air sub-compartment
        (i.e., urban, suburban and rural).


        :return: Pollutants emission per km driven, per air compartment.
        :rtype: numpy.array
        """

        em_arr = np.zeros((len(self.cycle), 11))

        if powertrain_type == "diesel":
            # HC
            em_arr[:, 0] = (
                -7.58e-09 * self.cycle ** 3
                + 2.817e-06 * self.cycle ** 2
                - 0.0003648 * self.cycle
                + 0.02031
            )
            # CO
            em_arr[:, 1] = (
                -1.058e-07 * self.cycle ** 3
                + 3.487e-05 * self.cycle ** 2
                - 0.003944 * self.cycle
                + 0.1602
            )
            # NOx + NO2
            em_arr[:, 2] = (
                8.348e-07 * self.cycle ** 3
                - 0.0001636 * self.cycle ** 2
                + 0.009568 * self.cycle
                - 0.09659
            ) + (
                2.504e-07 * self.cycle ** 3
                - 4.909e-05 * self.cycle ** 2
                + 0.00287 * self.cycle
                - 0.02898
            )
            # PM
            em_arr[:, 3] = (
                -3.079e-09 * self.cycle ** 3
                + 9.476e-07 * self.cycle ** 2
                - 9.87e-05 * self.cycle
                + 0.00504
            )
            # CH4
            em_arr[:, 4] = (
                -1.819e-10 * self.cycle ** 3
                + 6.76e-08 * self.cycle ** 2
                - 8.756e-06 * self.cycle
                + 0.0004875
            )
            # NMVOC
            em_arr[:, 5] = (
                -7.398e-09 * self.cycle ** 3
                + 2.749e-06 * self.cycle ** 2
                - 0.0003561 * self.cycle
                + 0.01983
            )
            # Pb
            # No Pb!
            # SO2
            em_arr[:, 7] = (
                -9.686e-10 * self.cycle ** 3
                + 3.552e-07 * self.cycle ** 2
                - 3.764e-05 * self.cycle
                + 0.00192
            )
            # N2O
            em_arr[:, 8] = (
                -1.468e-08 * self.cycle ** 3
                + 4.082e-06 * self.cycle ** 2
                - 0.0003563 * self.cycle
                + 0.01366
            )
            # NH3
            em_arr[:, 9] = -3.638e-23 * self.cycle ** 2 + 5.807e-21 * self.cycle + 0.001
            # Benzene
            em_arr[:, 10] = (
                -1.266e-10 * self.cycle ** 3
                + 4.704e-08 * self.cycle ** 2
                - 6.093e-06 * self.cycle
                + 0.0003392
            )

        if powertrain_type == "petrol":
            # HC
            em_arr[:, 0] = (
                3.377e-8 * self.cycle ** 3
                - 5.513e-6 * self.cycle ** 2
                + 2.463e-4 * self.cycle
                + 3.095e-3
            )
            # CO
            em_arr[:, 1] = (
                3.787e-6 * self.cycle ** 3
                - 6.65e-4 * self.cycle ** 2
                + 3.736e-2 * self.cycle
                + 0.4161
            )
            # NOx + NO2
            em_arr[:, 2] = (
                2.345e-09 * self.cycle ** 3
                + 1.445e-06 * self.cycle ** 2
                - 0.0004946 * self.cycle
                + 0.04578
            ) + (
                1.173e-10 * self.cycle ** 3
                + 7.226e-08 * self.cycle ** 2
                - 2.473e-05 * self.cycle
                + 0.002289
            )
            # PM
            em_arr[:, 3] = (
                1.46e-08 * self.cycle ** 3
                - 2.792e-06 * self.cycle ** 2
                + 0.0001596 * self.cycle
                - 0.001391
            )
            # CH4
            em_arr[:, 4] = (
                2.837e-09 * self.cycle ** 3
                - 4.631e-07 * self.cycle ** 2
                + 2.069e-05 * self.cycle
                + 0.00026
            )
            # NMVOC
            em_arr[:, 5] = (
                3.093e-08 * self.cycle ** 3
                - 5.05e-06 * self.cycle ** 2
                + 0.0002256 * self.cycle
                + 0.002835
            )
            # Pb
            em_arr[:, 6] = (
                -4.746e-11 * self.cycle ** 3
                + 2.063e-08 * self.cycle ** 2
                - 2.29e-06 * self.cycle
                + 0.000122
            )
            # SO2
            em_arr[:, 7] = (
                -7.512e-10 * self.cycle ** 3
                + 3.266e-07 * self.cycle ** 2
                - 3.624e-05 * self.cycle
                + 0.001932
            )
            # N2O
            em_arr[:, 8] = (
                -4.136e-09 * self.cycle ** 3
                + 1.163e-06 * self.cycle ** 2
                - 0.0001038 * self.cycle
                + 0.003129
            )
            # NH3
            em_arr[:, 9] = (
                -2.04e-06 * self.cycle ** 2 + 0.0009044 * self.cycle - 0.01251
            )
            # Benzene
            em_arr[:, 10] = (
                4.367e-09 * self.cycle ** 3
                - 7.128e-07 * self.cycle ** 2
                + 3.185e-05 * self.cycle
                + 0.0004002
            )

        if powertrain_type == "CNG":
            # HC
            em_arr[:, 0] = (
                5.064e-08 * self.cycle ** 3
                - 8.279e-06 * self.cycle ** 2
                + 0.0003724 * self.cycle
                + 0.004453
            )
            # CO
            em_arr[:, 1] = (
                3.755e-06 * self.cycle ** 3
                - 0.0006603 * self.cycle ** 2
                + 0.03719 * self.cycle
                - 0.4194
            )
            # NOx + NO2
            em_arr[:, 2] = (
                3.518e-09 * self.cycle ** 3
                + 2.168e-06 * self.cycle ** 2
                - 0.0007419 * self.cycle
                + 0.06867
            ) + (
                2.932e-10 * self.cycle ** 3
                + 1.807e-07 * self.cycle ** 2
                - 6.182e-05 * self.cycle
                + 0.005723
            )
            # PM
            em_arr[:, 3] = (
                1.46e-08 * self.cycle ** 3
                - 2.792e-06 * self.cycle ** 2
                + 0.0001596 * self.cycle
                - 0.001391
            )
            # CH4
            em_arr[:, 4] = (
                3.247e-08 * self.cycle ** 3
                - 5.312e-06 * self.cycle ** 2
                + 0.0002397 * self.cycle
                + 0.002803
            )
            # NMVOC
            em_arr[:, 5] = (
                1.817e-08 * self.cycle ** 3
                - 2.967e-06 * self.cycle ** 2
                + 0.0001328 * self.cycle
                + 0.00165
            )
            # Pb
            em_arr[:, 6] = (
                -4.509e-11 * self.cycle ** 3
                + 1.96e-08 * self.cycle ** 2
                - 2.175e-06 * self.cycle
                + 0.0001159
            )
            # SO2
            em_arr[:, 7] = (
                -7.138e-10 * self.cycle ** 3
                + 3.103e-07 * self.cycle ** 2
                - 3.444e-05 * self.cycle
                + 0.001835
            )
            # N2O
            em_arr[:, 8] = (
                -4.118e-09 * self.cycle ** 3
                + 1.159e-06 * self.cycle ** 2
                - 0.0001035 * self.cycle
                + 0.003123
            )
            # NH3
            em_arr[:, 9] = (
                -2.095e-06 * self.cycle ** 2 + 0.0009158 * self.cycle - 0.01293
            )
            # Benzene
            em_arr[:, 10] = (
                4.367e-09 * self.cycle ** 3
                - 7.128e-07 * self.cycle ** 2
                + 3.185e-05 * self.cycle
                + 0.0004002
            )

        distance = self.cycle.sum() / 3600
        dist_urban = self.cycle[self.cycle <= 50].sum() / 3600
        dist_suburban = self.cycle[(self.cycle > 50) & (self.cycle <= 80)].sum() / 3600
        dist_rural = self.cycle[self.cycle > 80].sum() / 3600

        urban = np.mean(em_arr[self.cycle <= 50], axis=0) * (dist_urban / distance)
        urban /= 1000  # going from grams to kg
        suburban = np.mean(em_arr[(self.cycle > 50) & (self.cycle <= 80)], axis=0) * (
            dist_suburban / distance
        )
        suburban /= 1000  # going from grams to kg
        rural = np.mean(em_arr[self.cycle > 80], axis=0) * (dist_rural / distance)
        rural /= 1000  # going from grams to kg

        suburban = np.nan_to_num(suburban)
        rural = np.nan_to_num(rural)

        return (urban, suburban, rural)
