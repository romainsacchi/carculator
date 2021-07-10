import numexpr as ne
import numpy as np
from scipy.interpolate import interp1d

from . import DATA_DIR


class InternalNoiseModel:
    """
    Calculate internal noise function ofpowertrain and size, based on http://www.auto-decibel-db.com/index.html.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(self, cycle):

        self.cycle = cycle

        self.noise_coeff = self.get_noise_coefficients()

    @staticmethod
    def get_noise_coefficients():

        filename = "internal_noise_coefficients.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of noise coefficients could not be found."
            )

        with open(filepath) as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
            (_, _, *header), *data = csv_list

            csv_dict = {}
            for row in data:
                key, sub_key, *values = row
                values = [float(v) for v in values]
                if key not in csv_dict:
                    csv_dict[key] = {sub_key: values}
                else:
                    csv_dict[key][sub_key] = values

        list_pt = [
            "ICEV-p",
            "ICEV-d",
            "ICEV-g",
            "PHEV-p",
            "PHEV-d",
            "FCEV",
            "BEV",
            "HEV-p",
            "HEV-d",
        ]
        list_size = ["Large", "Lower medium", "Medium", "Mini", "SUV", "Small", "Van"]

        arr = np.zeros((len(list_size), len(list_pt), 6))

        for pt in csv_dict:
            for size in csv_dict[pt]:
                arr[list_size.index(size), list_pt.index(pt), :] = csv_dict[pt][size]

        return arr

    def calculate_noise(self):
        # Instantiate the interpolation function
        # x = speed points (i.e., idle, 50 km/h, 80 km/h, 100 km/h, 120 km/h and 140 km/h)
        # y = dBs
        f = interp1d([0.0, 50.0, 80.0, 100.0, 120.0, 140.0], self.noise_coeff)

        # get dBs for every second of the driving cycle
        noise = f(self.cycle)

        # convert dBs to Watts (or joule.s^-1)
        noise = (10 ** -12) * (10 ** (noise / 10))

        # sum dBs along driving cycle to get joules
        noise = noise.sum(axis=2)

        # return a 9 by 7 arrays
        noise[noise < 1e-8] = np.nan
        return noise.reshape(7, 9, 1, 1)
