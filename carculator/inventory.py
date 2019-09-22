"""
.. module: inventory.py

"""
import numpy as np
from pathlib import Path
from inspect import currentframe, getframeinfo
import csv
from .hot_emissions import HotEmissionsModel
import itertools


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

    def __init__(self, array, cycle, method = 'recipe', level = 'midpoint', split = 'components'):

        self.year = array.coords['year'].values
        self.powertrain = array.coords['powertrain'].values
        self.size = array.coords['size'].values
        self.number_of_cars = len(self.year) * len(self.size) * len(self.powertrain)
        self.array = array.stack(desired=['powertrain', 'size', 'year'])
        self.cycle = cycle
        self.method = method
        self.level = level
        self.split = split

        self.inputs = self.get_dict_input()
        self.rev_inputs = self.get_rev_dict_input()
        self.A = self.get_A_matrix()
        self.B = self.get_B_matrix()

        self.index_cng = [self.inputs[i] for i in self.inputs if 'ICEV-g' in i[0]]
        self.index_combustion_wo_cng = [self.inputs[i] for i in self.inputs if any(ele in i[0] for ele in\
                                                                            ["ICEV-p", "HEV-p", "PHEV", "ICEV-d"])]
        self.index_diesel = [self.inputs[i] for i in self.inputs if 'ICEV-d' in i[0]]
        self.index_all_petrol = [self.inputs[i] for i in self.inputs if any(ele in i[0] for ele in\
                                                                            ["ICEV-p", "HEV-p", "PHEV"])]
        self.index_petrol = [self.inputs[i] for i in self.inputs if 'ICEV-p' in i[0]]
        self.index_hybrid = [self.inputs[i] for i in self.inputs if 'HEV-p' in i[0]]
        self.index_plugin_hybrid = [self.inputs[i] for i in self.inputs if 'PHEV' in i[0]]
        self.index_fuel_cell = [self.inputs[i] for i in self.inputs if 'FCEV' in i[0]]
        self.index_emissions = [self.inputs[i] for i in self.inputs if 'air' in i[1][0] and len(i[1])>1]
        self.index_noise = [self.inputs[i] for i in self.inputs if 'noise' in i[0]]

        f = np.zeros(np.shape(self.A)[0])
        f[-1] = 1

        for i in array.value:
            self.temp_array = self.array.sel(value=i.values)
            self.set_inputs_in_A_matrix()

            g = np.linalg.inv(self.A).dot(f)
            C = g.dot(self.B.T)
            print(C)


    def __getitem__(self, key):
        """
        Make class['foo'] automatically filter for the parameter 'foo'
        Makes the model code much cleaner

        :param key: Parameter name
        :type key: str
        :return: `array` filtered after the parameter selected
        """
        return self.temp_array.sel(parameter=key)

    def get_A_matrix(self):
        filename = 'A_matrix.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The technology matrix could not be found.')

        initial_A = np.genfromtxt(filepath, delimiter=';')
        new_A = np.identity(np.shape(initial_A)[0] + self.number_of_cars)
        new_A[0:np.shape(initial_A)[0],0:np.shape(initial_A)[0]] = initial_A

        return new_A

    def get_B_matrix(self):
        filename = 'B_matrix' + '_' + self.method + '_' + self.level + '.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The method or impact level specified could not be found.')

        initial_B = np.genfromtxt(filepath, delimiter=';')
        print(np.shape(initial_B))
        new_B = np.zeros((np.shape(initial_B)[0],np.shape(initial_B)[1] + self.number_of_cars))
        new_B[0:np.shape(initial_B)[0],0:np.shape(initial_B)[1]] = initial_B
        print(np.shape(new_B))

        return new_B

    def get_dict_input(self):
        filename = 'dict_inputs_A_matrix.csv'
        filepath = Path(getframeinfo(currentframe()).filename).resolve().parent.joinpath('data/'+ filename)
        if not filepath.is_file():
            raise FileNotFoundError('The dictionary of activity labels could not be found.')

        csv_dict = {}
        count = 0

        with open(filepath) as f:
            input_dict = csv.reader(f, delimiter=';')
            for row in input_dict:
                if "(" in row[1]:
                    new_str = row[1].replace("(","")
                    new_str = new_str.replace(")","")
                    new_str = [s.strip() for s in new_str.split(',') if s]
                    t = ()
                    for s in new_str:
                        if "low population" in s:
                            s = 'low population density, long-term'
                            t += (s,)
                            break
                        else:
                            t += (s.replace("'", ''),)
                    csv_dict[(row[0], t)] = count
                else:
                    csv_dict[(row[0], row[1])] = count
                count+=1

        for pt in self.powertrain:
            for s in self.size:
                for y in self.year:
                    maximum = csv_dict[max(csv_dict, key = csv_dict.get)]
                    name = 'Passenger car, ' + pt + ', ' + s + ', ' + str(y)
                    csv_dict[(name, 'GLO')] = maximum + 1


        return csv_dict

    def get_rev_dict_input(self):
        return {v: k for k,v in self.inputs.items()}

    def set_inputs_in_A_matrix(self):
        print('done')

        # Glider
        self.A[self.inputs[('market for glider, passenger car', 'GLO')], -self.number_of_cars :] =\
            self['glider base mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Glider lightweighting', 'GLO')], -self.number_of_cars :] =\
            self['lightweighting'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('maintenance, passenger car', 'RER')], -self.number_of_cars :] =\
            self['curb mass'] / 1600 / 150000 *-1

        # Glider EoL
        self.A[self.inputs[('market for manual dismantling of used electric passenger car','GLO')], -self.number_of_cars :] =\
            self['curb mass'] * (1 - self['combustion power share']) / 1180 / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for manual dismantling of used passenger car with internal combustion engine','GLO')], -self.number_of_cars :] =\
            self['curb mass'] * self['combustion power share'] / 1600 / self['lifetime kilometers'] *-1

        # Powertrain components
        self.A[self.inputs[('market for charger, electric passenger car','GLO')], -self.number_of_cars :] =\
            self['charger mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for converter, for electric passenger car','GLO')], -self.number_of_cars :] =\
            self['converter mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for electric motor, electric passenger car','GLO')], -self.number_of_cars :] =\
            self['electric engine mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for inverter, for electric passenger car','GLO')], -self.number_of_cars :] =\
            self['inverter mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for power distribution unit, for electric passenger car','GLO')], -self.number_of_cars :] =\
            self['power distribution unit mass'] / self['lifetime kilometers'] *-1

        l_elec_pt = ['charger mass','converter mass','inverter mass','power distribution unit mass',
            'electric engine mass', 'fuel cell stack mass', 'fuel cell ancillary BoP mass', 'fuel cell essential BoP mass',
                     'battery cell mass','battery BoP mass']

        self.A[self.inputs[('market for used powertrain from electric passenger car, manual dismantling','GLO')], -self.number_of_cars :] =\
            self[l_elec_pt].sum(axis=0) / self['lifetime kilometers'] *-1

        self.A[self.inputs[('market for internal combustion engine, passenger car','GLO')], -self.number_of_cars :] =\
            (self[['combustion engine mass','electric engine mass', 'powertrain mass']].sum(axis=0)) / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Ancillary BoP','GLO')], -self.number_of_cars :] =\
            self['fuel cell ancillary BoP mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Essential BoP','GLO')], -self.number_of_cars :] =\
            self['fuel cell essential BoP mass'] / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Stack 2020','GLO')], -self.number_of_cars :] =\
            self['fuel cell stack mass'] / self['lifetime kilometers'] *-1

        # Energy storage
        self.A[self.inputs[('Battery BoP','GLO')], -self.number_of_cars :] =\
            (self['battery BoP mass'] * (1 + self['battery lifetime replacements'])) / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Battery cell','GLO')], -self.number_of_cars :] =\
            (self['battery cell mass'] * (1 + self['fuel cell lifetime replacements'])) / self['lifetime kilometers'] *-1

        self.A[self.inputs[('Battery cell','GLO')], -self.number_of_cars :] =\
            (self['battery cell mass'] * (1 + self['fuel cell lifetime replacements'])) / self['lifetime kilometers'] *-1

        self.A[self.inputs[('polyethylene production, high density, granulate','RER')], self.index_petrol] =\
            self.temp_array.sel(parameter = 'fuel tank mass', powertrain = 'ICEV-p') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'ICEV-p') *-1

        self.A[self.inputs[('polyethylene production, high density, granulate','RER')], self.index_hybrid] =\
            self.temp_array.sel(parameter = 'fuel tank mass', powertrain = 'HEV-p') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'HEV-p') *-1

        self.A[self.inputs[('polyethylene production, high density, granulate','RER')], self.index_plugin_hybrid] =\
            self.temp_array.sel(parameter = 'fuel tank mass', powertrain = 'PHEV') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'PHEV') *-1

        self.A[self.inputs[('glass fibre reinforced plastic production, polyamide, injection moulded','RER')], self.index_cng] =\
           self.temp_array.sel(parameter = 'fuel tank mass', powertrain = 'ICEV-g') / self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'ICEV-g') *-1

        self.A[self.inputs[('Fuel tank, compressed hydrogen gas, 700bar','GLO')], self.index_fuel_cell] =\
            self.temp_array.sel(parameter = 'fuel tank mass', powertrain = 'FCEV') / self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'FCEV') *-1

        # Energy chain
        self.A[self.inputs[('market group for electricity, low voltage','ENTSO-E')], -self.number_of_cars :] =\
            self['electricity consumption'] *-1

        self.A[self.inputs[("Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station",'RER')], self.index_fuel_cell] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'FCEV') / self.temp_array.sel(parameter = 'range', powertrain = 'FCEV') *-1

        self.A[self.inputs[("market for natural gas, from high pressure network (1-5 bar), at service station",'GLO')], self.index_cng] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'ICEV-g') / self.temp_array.sel(parameter = 'range', powertrain = 'ICEV-g') *-1

        self.A[self.inputs[("market for diesel",'Europe without Switzerland')], self.index_diesel] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'ICEV-d') / self.temp_array.sel(parameter = 'range', powertrain = 'ICEV-d') *-1

        self.A[self.inputs[('market for petrol, low-sulfur','Europe without Switzerland')], self.index_petrol] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'ICEV-p') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'ICEV-p') *-1

        self.A[self.inputs[('market for petrol, low-sulfur','Europe without Switzerland')], self.index_hybrid] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'HEV-p') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'HEV-p') *-1

        self.A[self.inputs[('market for petrol, low-sulfur','Europe without Switzerland')], self.index_plugin_hybrid] =\
            self.temp_array.sel(parameter = 'fuel mass', powertrain = 'PHEV') /self.temp_array.sel(parameter = 'lifetime kilometers', powertrain = 'PHEV') *-1

        # Non-exhaust emissions
        self.A[self.inputs[("market for road wear emissions, passenger car",'GLO')], -self.number_of_cars :] =\
            self['driving mass'] * -1 * 1E-08 *-1
        self.A[self.inputs[("market for tyre wear emissions, passenger car",'GLO')], -self.number_of_cars :] =\
            self['driving mass'] * -1 * 6E-08 *-1
        self.A[self.inputs[("market for brake wear emissions, passenger car",'GLO')], -self.number_of_cars :] =\
            self['driving mass'] * -1 * 5E-09 *-1

        # Infrastructure
        self.A[self.inputs[("market for road",'GLO')], -self.number_of_cars :] =\
            5.37E-7 * self['driving mass'] *-1

        # Exhaust emissions
        # Fuel-based emissions
        self.A[self.inputs[("Carbon dioxide, fossil", ('air',))], -self.number_of_cars :] =\
            (self['CO2 per kg fuel'] * self['fuel mass'])/ self['range'] *-1

        # Non-fuel based emissions
        list_direct_emissions = ["Benzene direct emissions, urban",
            "Benzene direct emissions, suburban",
            "Benzene direct emissions, rural",
            "Sulfur dioxide direct emissions, urban",
            "Sulfur dioxide direct emissions, suburban",
            "Sulfur dioxide direct emissions, rural",
            "Methane direct emissions, urban",
            "Methane direct emissions, suburban",
            "Methane direct emissions, rural",
            "Carbon monoxide direct emissions, urban",
            "Carbon monoxide direct emissions, suburban",
            "Carbon monoxide direct emissions, rural",
            "Hydrocarbons direct emissions, urban",
            "Hydrocarbons direct emissions, suburban",
            "Hydrocarbons direct emissions, rural",
            "Dinitrogen oxide direct emissions, urban",
            "Dinitrogen oxide direct emissions, suburban",
            "Dinitrogen oxide direct emissions, rural",
            "Ammonia direct emissions, urban",
            "Ammonia direct emissions, suburban",
            "Ammonia direct emissions, rural",
            "NMVOC direct emissions, urban",
            "NMVOC direct emissions, suburban",
            "NMVOC direct emissions, rural",
            "Nitrogen oxides direct emissions, urban",
            "Nitrogen oxides direct emissions, suburban",
            "Nitrogen oxides direct emissions, rural",
            "Particulate matters direct emissions, urban",
            "Particulate matters direct emissions, suburban",
            "Particulate matters direct emissions, rural",
            "Lead direct emissions, urban",
            "Lead direct emissions, suburban",
            "Lead direct emissions, rural"]

        self.A[self.index_emissions, -self.number_of_cars :] = self[list_direct_emissions] *-1

        # Noise emissions
        self.A[self.index_noise, -self.number_of_cars :] = self[[self.rev_inputs[e][0] for e in self.index_noise]] *-1




