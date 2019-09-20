"""
.. module: inventory.py

"""
import numpy as np
from pathlib import Path
from inspect import currentframe, getframeinfo


class InventoryCalculation:
    """
    Build and solve the inventory for results characterization

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series
    :param rho_air: Mass per unit volume of air. Set to (1.225 kg/m3) by default.
    :type rho_air: float

    :ivar rho_air: Mass per unit volume of air. Value of 1.204 at 23C (test temperature for WLTC).
    :vartype rho_air: float
    :ivar velocity: Time series of speed values, in meters per second.
    :vartype velocity: numpy.ndarray
    :ivar acceleration: Time series of acceleration, calculated as increment in velocity per interval of 1 second,
        in meter per second^2.
    :vartype acceleration: numpy.ndarray
    """

    def __init__(self, array, method = 'recipe', level = 'midpoint', split = 'components'):

        self.array = array
        self.method = method
        self.level = level
        self.split = split
        self.number_of_cars = len(self.array.coords['size']) * len(self.array.coords['powertrain']) \
                              * len(self.array.coords['year'])

    def get_A_matrix(self):
        filename = 'A_matrix.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The technology matrix could not be found.')

        initial_A = np.genfromtxt(filepath, delimiter=';')
        new_A = np.identity(len(initial_A) + self.number_of_cars)
        new_A[0:len(initial_A),0:len(initial_A)] = initial_A
        return new_A

    def get_B_matrix(self):
        filename = 'B_matrix' + self.method + '_' + self.level + '.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The method or impact level specified could not be found.')

        initial_B = np.genfromtxt(filepath, delimiter=';')
        new_B = np.identity(len(initial_B) + self.number_of_cars)
        new_B[0:len(initial_B),0:len(initial_B)] = initial_B
        return new_B

    def get_dict_input(self):
        filename = 'dict_inputs_A_matrix.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The dictionary of activity labels could not be found.')

        with open(filename) as f:
            input_dict = [[val.strip() for val in r.split(";")] for r in f.readlines()]

        (_, *header), *data = input_dict
        csv_dict = {}
        for row in data:
            key, *values = row
            csv_dict[key] = {key: value for key, value in zip(header, values)}
        print(input_dict)
        for pt in self.array.coords['powertrain']:
            for s in self.array.coords['size']:
                for y in self.array.coords['year']:
                    maximum = max(csv_dict, key = csv_dict.get)
                    csv_dict['Passenger car, ' + pt + ', ' + s + ', ' + y] = maximum + 1

        return csv_dict

    def get_rev_dict_input(self):
        return {v: k for k,v in self.inputs.items()}



