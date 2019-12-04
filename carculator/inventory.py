"""
.. module: inventory.py

"""
import numpy as np
from pathlib import Path
from inspect import currentframe, getframeinfo
import csv
import xarray as xr
from . import DATA_DIR
from .background_systems import BackgroundSystemModel
import itertools

class InventoryCalculation:
    """
    Build and solve the inventory for results characterization

    """

    def __init__(self, array):

        self.year = array.coords["year"].values
        self.powertrain = array.coords["powertrain"].values
        self.size = array.coords["size"].values
        self.number_of_cars = len(self.year) * len(self.size) * len(self.powertrain)
        self.array = array.stack(desired=["powertrain", "size", "year"])
        self.iterations = len(array.value.values)

        self.inputs = self.get_dict_input()
        self.rev_inputs = self.get_rev_dict_input()
        self.array_inputs = {
            x: i for i, x in enumerate(list(self.array.parameter.values), 0)
        }
        self.array_powertrains = {
            x: i for i, x in enumerate(list(self.array.powertrain.values), 0)
        }
        self.A = self.get_A_matrix()

        self.index_cng = [self.inputs[i] for i in self.inputs if "ICEV-g" in i[0]]
        self.index_combustion_wo_cng = [
            self.inputs[i]
            for i in self.inputs
            if any(ele in i[0] for ele in ["ICEV-p", "HEV-p", "PHEV", "ICEV-d"])
        ]
        self.index_diesel = [self.inputs[i] for i in self.inputs if "ICEV-d" in i[0]]
        self.index_all_petrol = [
            self.inputs[i]
            for i in self.inputs
            if any(ele in i[0] for ele in ["ICEV-p", "HEV-p", "PHEV"])
        ]
        self.index_petrol = [self.inputs[i] for i in self.inputs if "ICEV-p" in i[0]]
        self.index_hybrid = [self.inputs[i] for i in self.inputs if "HEV-p" in i[0]]
        self.index_plugin_hybrid = [
            self.inputs[i] for i in self.inputs if "PHEV" in i[0]
        ]
        self.index_fuel_cell = [self.inputs[i] for i in self.inputs if "FCEV" in i[0]]
        self.index_emissions = [
            self.inputs[i] for i in self.inputs if "air" in i[1][0] and len(i[1]) > 1
                                                   and i[0] not in ['Carbon dioxide, fossil',
                                                                    'Carbon monoxide, non-fossil',
                                                                    'Methane, non-fossil',
                                                                    'Particulates, > 10 um']
        ]
        self.index_noise = [self.inputs[i] for i in self.inputs if "noise" in i[0]]
        self.list_cat, self.split_indices = self.get_split_indices()
        self.bs = BackgroundSystemModel()

    def __getitem__(self, key):
        """
        Make class['foo'] automatically filter for the parameter 'foo'
        Makes the model code much cleaner

        :param key: Parameter name
        :type key: str
        :return: `array` filtered after the parameter selected
        """
        return self.temp_array.sel(parameter=key)

    def get_results_table(self, FU, method, level, split):

        if split == "components":
            cat = [
                "direct",
                "energy chain",
                "maintenance",
                "glider",
                "powertrain",
                "energy storage",
                "road",
                "other"
            ]

        dict_impact_cat = self.get_dict_impact_categories()

        response = xr.DataArray(
            np.zeros(
                (
                    self.B.shape[0],
                    len(FU['size']),
                    len(FU['powertrain']),
                    len(FU['year']),
                    len(cat),
                    self.iterations,
                )
            ),
            coords=[
                dict_impact_cat[method][level],
                FU['size'],
                FU['powertrain'],
                FU['year'],
                cat,
                np.arange(0, self.iterations),
            ],
            dims=[
                "impact_category",
                "size",
                "powertrain",
                "year",
                "impact",
                "value",
            ],
        )

        return response

    def get_split_indices(self):
        filename = "dict_split.csv"
        filepath = (DATA_DIR / filename)
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of splits could not be found."
            )

        with open(filepath) as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        (_, _, *header), *data = csv_list

        csv_dict = {}
        for row in data:
            key, sub_key, *values = row

            if key in csv_dict:
                if sub_key in csv_dict[key]:
                    csv_dict[key][sub_key].append({'search by': values[0], 'search for': values[1]})
                else:
                    csv_dict[key][sub_key] = [{'search by': values[0], 'search for': values[1]}]
            else:
                csv_dict[key] = {sub_key: [{'search by': values[0], 'search for': values[1]}]}

        flatten = itertools.chain.from_iterable

        d = {}
        l=[]
        for cat in csv_dict['components']:
            d[cat] = list(flatten(
                [self.get_index_of_flows([l['search for']], l['search by']) for l in csv_dict['components'][cat]]))
            l.append(d[cat])

        d['other'] = [x for x in self.rev_inputs
                      if not x in list(flatten(l))]

        list_ind = [d[x] for x in d]
        maxLen = max(map(len, list_ind))
        for row in list_ind:
            while len(row) < maxLen:
                row.extend([len(self.inputs)-1])
        return list(d.keys()), list_ind

    def calculate_impacts(
        self, scope=None, background_configuration = None, method="recipe", level="midpoint", split="components"
    ):
        if scope is None:
            scope={}
            scope['size'] = self.size.tolist()
            scope['powertrain'] = self.powertrain.tolist()
            scope['year'] = self.year.tolist()
        else:
            scope['size'] = scope.get('size', self.size.tolist())
            scope['powertrain'] = scope.get('powertrain', self.powertrain.tolist())
            scope['year'] = scope.get('year', self.year.tolist())

        # Load the B matrix
        self.B = self.get_B_matrix(method, level)

        # Prepare an array to store the results
        self.results = self.get_results_table(scope, method, level, split)

        # Iterate through the number of iterations
        for i in range(self.iterations):
            self.temp_array = self.array.sel(value=i).values
            self.set_inputs_in_A_matrix(background_configuration, i)

            for pt in scope['powertrain']:
                for y in scope['year']:
                    for s in scope['size']:
                        # Retrieve the index of a given car in the matrix A
                        car = self.inputs[('Passenger car, '+ pt + ', ' + s + ', ' + str(y),"GLO")]
                        # Set the demand vector with zeros and a 1 corresponding to the car position in the vector
                        f = np.zeros(np.shape(self.A)[0])
                        f[car] = 1

                        # Solve inventory
                        C = np.linalg.solve(self.A, f) * self.B

                        if pt=='BEV' and y==2018 and s=='Large':
                            np.savetxt('C.csv', C)
                        # Iterate through the results array to fill it
                        # TODO: optimize this section
                        self.results.loc[
                            dict(impact=self.list_cat, year=y, size=s, powertrain=pt, value=i)] = \
                            C[:, self.split_indices].sum(axis=2)

        return self.results

    def get_A_matrix(self):
        filename = "A_matrix.csv"
        filepath = (
            Path(getframeinfo(currentframe()).filename)
            .resolve()
            .parent.joinpath("data/" + filename)
        )
        if not filepath.is_file():
            raise FileNotFoundError("The technology matrix could not be found.")

        initial_A = np.genfromtxt(filepath, delimiter=";")
        new_A = np.identity(np.shape(initial_A)[0] + self.number_of_cars)
        new_A[0 : np.shape(initial_A)[0], 0 : np.shape(initial_A)[0]] = initial_A

        return new_A

    def get_B_matrix(self, method, level):
        filename = "B_matrix" + "_" + method + "_" + level + ".csv"
        filepath = (DATA_DIR / filename)

        if not filepath.is_file():
            raise FileNotFoundError(
                "The method or impact level specified could not be found."
            )

        initial_B = np.genfromtxt(filepath, delimiter=";")

        new_B = np.zeros(
            (np.shape(initial_B)[0], np.shape(initial_B)[1] + self.number_of_cars)
        )
        new_B[0 : np.shape(initial_B)[0], 0 : np.shape(initial_B)[1]] = initial_B

        return new_B

    def get_dict_input(self):
        filename = "dict_inputs_A_matrix.csv"
        filepath = (DATA_DIR / filename)
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of activity labels could not be found."
            )
        csv_dict = {}
        count = 0
        with open(filepath) as f:
            input_dict = csv.reader(f, delimiter=";")
            for row in input_dict:
                if "(" in row[1]:
                    new_str = row[1].replace("(", "")
                    new_str = new_str.replace(")", "")
                    new_str = [s.strip() for s in new_str.split(",") if s]
                    t = ()
                    for s in new_str:
                        if "low population" in s:
                            s = "low population density, long-term"
                            t += (s,)
                            break
                        else:
                            t += (s.replace("'", ""),)
                    csv_dict[(row[0], t, row[2])] = count
                else:
                    csv_dict[(row[0], row[1], row[2], row[3])] = count
                count += 1

        for pt in self.powertrain:
            for s in self.size:
                for y in self.year:
                    maximum = csv_dict[max(csv_dict, key=csv_dict.get)]
                    name = "Passenger car, " + pt + ", " + s + ", " + str(y)
                    csv_dict[(name, "GLO")] = maximum + 1

        return csv_dict

    def get_dict_impact_categories(self):
        filename = "dict_impact_categories.csv"
        filepath = (DATA_DIR / filename)
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of impact categories could not be found."
            )

        csv_dict = {}

        with open(filepath) as f:
            input_dict = csv.reader(f, delimiter=";")
            for row in input_dict:
                csv_dict[row[0]] = {row[1]: [x.strip() for x in row[2:]]}

        return csv_dict

    def get_rev_dict_input(self):
        return {v: k for k, v in self.inputs.items()}

    def get_index_from_array(self, items_to_look_for):
        return [
            self.inputs[c] - (len(self.inputs) - self.number_of_cars)
            for c in self.inputs
            if any(ele in c[0] for ele in items_to_look_for)
        ]

    def get_index_of_flows(self, items_to_look_for, search_by="name"):
        if search_by == "name":
            return [
                int(self.inputs[c])
                for c in self.inputs
                if any(ele in c[0] for ele in items_to_look_for)
            ]
        if search_by == "compartment":
            return [
                int(self.inputs[c])
                for c in self.inputs
                if any(ele in c[1] for ele in items_to_look_for)
            ]

    def set_inputs_in_A_matrix(self, background_configuration, iteration):

        # Glider
        self.A[
            self.inputs[('market for glider, passenger car', 'GLO','kilogram','glider, passenger car')],
            -self.number_of_cars :,
        ] = (
            (self.temp_array[self.array_inputs["glider base mass"], :] * (1-self.temp_array[self.array_inputs["lightweighting"], :]))
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[('Glider lightweighting', 'GLO', 'kilogram', 'Glider lightweighting')], -self.number_of_cars :
        ] = (
            (self.temp_array[self.array_inputs["lightweighting"], :] * self.temp_array[self.array_inputs["glider base mass"], :])
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[('maintenance, passenger car', 'RER', 'unit', 'passenger car maintenance')], -self.number_of_cars :
        ] = (self.temp_array[self.array_inputs["curb mass"], :] / 1600 / 150000 * -1)

        # Glider EoL
        self.A[
            self.inputs[
                ('market for manual dismantling of used electric passenger car',
                 'GLO',
                 'unit',
                 'manual dismantling of used electric passenger car')
            ],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["curb mass"], :]
            * (1 - self.temp_array[self.array_inputs["combustion power share"], :])
            / 1180
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[
                ('market for manual dismantling of used passenger car with internal combustion engine',
                 'GLO',
                 'unit',
                 'manual dismantling of used passenger car with internal combustion engine')
            ],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["curb mass"], :]
            * self.temp_array[self.array_inputs["combustion power share"], :]
            / 1600
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Powertrain components
        self.A[
            self.inputs[('market for charger, electric passenger car',
              'GLO',
              'kilogram',
              'charger, electric passenger car')],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["charger mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[('market for converter, for electric passenger car',
              'GLO',
              'kilogram',
              'converter, for electric passenger car')],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["converter mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[('market for electric motor, electric passenger car',
              'GLO',
              'kilogram',
              'electric motor, electric passenger car')],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["electric engine mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[('market for inverter, for electric passenger car',
              'GLO',
              'kilogram',
              'inverter, for electric passenger car')],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["inverter mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[
                ('market for power distribution unit, for electric passenger car',
                 'GLO',
                 'kilogram',
                 'power distribution unit, for electric passenger car')
            ],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["power distribution unit mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        l_elec_pt = [
            "charger mass",
            "converter mass",
            "inverter mass",
            "power distribution unit mass",
            "electric engine mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
            "battery cell mass",
            "battery BoP mass",
        ]

        self.A[
            self.inputs[
                ('market for used powertrain from electric passenger car, manual dismantling',
                 'GLO',
                 'kilogram',
                 'used powertrain from electric passenger car, manual dismantling')
            ],
            -self.number_of_cars :,
        ] = (
            self.temp_array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
        )

        self.A[
            self.inputs[
                ('market for internal combustion engine, passenger car',
                 'GLO',
                 'kilogram',
                 'internal combustion engine, for passenger car')
            ],
            -self.number_of_cars :,
        ] = (
            (
                self.temp_array[
                    [
                        self.array_inputs[l]
                        for l in [
                            "combustion engine mass",
                            "electric engine mass",
                            "powertrain mass",
                        ]
                    ],
                    :,
                ].sum(axis=0)
            )
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[('Ancillary BoP', 'GLO', 'kilogram', 'Ancillary BoP')], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell ancillary BoP mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[('Essential BoP', 'GLO', 'kilogram', 'Essential BoP')], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell essential BoP mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[('Stack 2020', 'GLO', 'kilowatt', 'Stack 2020')], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell stack mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Energy storage
        self.A[self.inputs[('Battery BoP', 'GLO', 'kilogram', 'Battery BoP')], -self.number_of_cars :] = (
            (
                self.temp_array[self.array_inputs["battery BoP mass"], :]
                * (
                    1
                    + self.temp_array[
                        self.array_inputs["battery lifetime replacements"], :
                    ]
                )
            )
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[('Battery cell', 'GLO', 'kilogram', 'Battery cell')], -self.number_of_cars :] = (
            (
                self.temp_array[self.array_inputs["battery cell mass"], :]
                * (
                    1
                    + self.temp_array[
                        self.array_inputs["fuel cell lifetime replacements"], :
                    ]
                )
            )
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if any(ele in c[0] for ele in ["ICEV-d", "ICEV-p", "HEV-p", "PHEV"])
        ]
        index = self.get_index_from_array(["ICEV-d", "ICEV-p", "HEV-p", "PHEV"])

        self.A[
            self.inputs[('polyethylene production, high density, granulate',
              'RER',
              'kilogram',
              'polyethylene, high density, granulate')],
            index_A,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        index = self.get_index_from_array(["ICEV-g"])
        self.A[
            self.inputs[
                ('glass fibre reinforced plastic production, polyamide, injection moulded',
                 'RER',
                 'kilogram',
                 'glass fibre reinforced plastic, polyamide, injection moulded')
            ],
            self.index_cng,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        index = self.get_index_from_array(["FCEV"])
        self.A[
            self.inputs[('Fuel tank, compressed hydrogen gas, 700bar',
              'GLO',
              'kilogram',
              'Fuel tank, compressed hydrogen gas, 700bar')],
            self.index_fuel_cell,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        # Energy chain
        dict_map = {
            'Hydro': ('electricity production, hydro, run-of-river',
                      'DE',
                      'kilowatt hour',
                      'electricity, high voltage',
                      ),
            'Nuclear': ('electricity production, nuclear, pressure water reactor',
                        'DE',
                        'kilowatt hour',
                        'electricity, high voltage'),
            'Gas': ('electricity production, natural gas, conventional power plant',
                    'DE',
                    'kilowatt hour',
                    'electricity, high voltage'),
            'Solar': ('electricity production, photovoltaic, 3kWp slanted-roof installation, multi-Si, panel, mounted',
                      'DE',
                      'kilowatt hour',
                      'electricity, low voltage'),
            'Wind': ('electricity production, wind, 1-3MW turbine, onshore',
                     'DE',
                     'kilowatt hour',
                     'electricity, high voltage'),
            'Biomass': ('heat and power co-generation, wood chips, 6667 kW, state-of-the-art 2014',
                        'DE',
                        'kilowatt hour',
                        'electricity, high voltage'),
            'Coal': ('electricity production, hard coal',
                     'DE',
                     'kilowatt hour',
                     'electricity, high voltage'),
            'Oil': ('electricity production, oil',
                    'DE',
                    'kilowatt hour',
                    'electricity, high voltage'),
            'Geo': ('electricity production, deep geothermal',
                    'DE',
                    'kilowatt hour',
                    'electricity, high voltage'),
            'Waste': (
                    'electricity, from municipal waste incineration to generic market for electricity, medium voltage',
                    'DE',
                    'kilowatt hour',
                    'electricity, medium voltage'),
        }
        if background_configuration:
            # If a customization dict is passed
            if 'country' in background_configuration:
                # If a country is specified
                country = background_configuration['country']
                losses_to_low = float(self.bs.losses[country]['LV'])

                if 'custom electricity mix' in background_configuration:
                    # If a special electricity mix is specified, we use it
                    mix = background_configuration['custom electricity mix']
                else:
                    mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.year).values
            else:
                country = 'RER'
                losses_to_low = 1.07
                mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.year).values
        else:
            country = 'RER'
            losses_to_low = 1.07
            mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.year).values


        for y in self.year:
            self.A[np.ix_([self.inputs[dict_map[t]] for t in dict_map],
                          [self.inputs[i] for i in self.inputs if str(y) in i[0]])] = \
                np.outer(mix[self.year.tolist().index(y)],
                         self.temp_array[
                             self.array_inputs["electricity consumption"],
                             self.get_index_from_array([str(y)])
                         ]) * -1 * losses_to_low


        index = self.get_index_from_array(["FCEV"])

        if background_configuration:
            if 'hydrogen technology' in background_configuration:
                # If a customization dict is passed
                hydro_technology = background_configuration['hydrogen technology']
            else:
                hydro_technology = 'Electrolysis'
        else:
            hydro_technology = 'Electrolysis'

        dict_h_map = {
            'Electrolysis': ('Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station',
                          'RER',
                          'kilogram',
                          'Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station'),
            'Electrolysis - solar': ('Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - solar PV',
                          'RER',
                          'kilogram',
                          'Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - solar PV'),
            'Electrolysis - hydro': ('Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - hydro reservoir',
                                  'RER',
                                  'kilogram',
                                  'Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - hydro reservoir'),
            'Electrolysis - nuclear': ('Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - nuclear PWR',
                                      'RER',
                                      'kilogram',
                                      'Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station - nuclear PWR'),
            'SMR': ('Hydrogen, gaseous, 700 bar, from SMR NG w/o CCS, at H2 fuelling station',
                  'RER',
                  'kilogram',
                  'Hydrogen, gaseous, 700 bar, from SMR NG w/o CCS, at H2 fuelling station'),
        }


        self.A[self.inputs[dict_h_map[hydro_technology]],
                self.index_fuel_cell] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
            )

        # If hydrolysis is chosen, adjust the electricity mix

        if hydro_technology == 'hydrolysis' and iteration == 0:

            # Zero out initial electricity provider
            old_amount = self.A[

                self.inputs[('market group for electricity, medium voltage',
                                              'Europe without Switzerland',
                                              'kilowatt hour',
                                              'electricity, medium voltage')],
                            self.inputs[dict_h_map[hydro_technology]]
            ]

            self.A[self.inputs[('market group for electricity, medium voltage',
              'Europe without Switzerland',
              'kilowatt hour',
              'electricity, medium voltage')],
               self.inputs[dict_h_map[hydro_technology]]] = 0


            # TODO: differentiate hydrogen production in time

            self.A[[self.inputs[dict_map[t]] for t in dict_map],
                    self.inputs[dict_h_map[hydro_technology]]
                    ] = \
                (np.outer(mix[0], old_amount) * losses_to_low).reshape(10,)

        index = self.get_index_from_array(["ICEV-g"])

        if background_configuration:
            if 'cng technology' in background_configuration:
                # If a customization dict is passed
                cng_technology = background_configuration['cng technology']
            else:
                cng_technology = 'cng'
        else:
            cng_technology = 'cng'

        dict_cng_map = {
            'cng': ('market for natural gas, from high pressure network (1-5 bar), at service station',
                  'GLO',
                  'kilogram',
                  'natural gas, from high pressure network (1-5 bar), at service station'),
            'biogas': ('biogas upgrading - sewage sludge - amine scrubbing - best',
                      'CH',
                      'cubic meter',
                      'biogas upgrading - sewage sludge - amine scrubbing - best')
        }

        if cng_technology == 'cng':
            self.A[self.inputs[dict_cng_map[cng_technology]],
                    self.index_cng] = (
                self.temp_array[self.array_inputs["fuel mass"], index]
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )
        else:
            # biogas
            self.A[self.inputs[dict_cng_map[cng_technology]],
                    self.index_cng] = (
                (self.temp_array[self.array_inputs["fuel mass"], index] / 180 ) #kg/m3 @ 200 bar
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )

        index = self.get_index_from_array(["ICEV-d"])
        if background_configuration:
            if 'diesel technology' in background_configuration:
                # If a customization dict is passed
                diesel_technology = background_configuration['diesel technology']
            else:
                diesel_technology = 'diesel'
        else:
            diesel_technology = 'diesel'

        dict_diesel_map = {
            'diesel': ('market for diesel', 'Europe without Switzerland', 'kilogram', 'diesel'),
            'biodiesel - algae': ('biodiesel production from algae',
                                  'RER',
                                  'megajoule',
                                  'biodiesel production from algae'),
            'biodiesel - cooking oil': ('Biodiesel ( Fatty Acid Methy Esters) production from Waste (Used) cooking Oil  {RER} | transesterification of WVO Europe | Alloc Rec, U',
                                      'RER',
                                      'megajoule',
                                      'Biodiesel ( Fatty Acid Methy Esters) production from Waste (Used) cooking Oil  {RER} | transesterification of WVO Europe | Alloc Rec, U')
        }

        if diesel_technology == 'diesel':
            self.A[self.inputs[dict_diesel_map[diesel_technology]],
                    self.index_diesel] = (
                self.temp_array[self.array_inputs["fuel mass"], index]
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )
        else:
            # biodiesel
            self.A[self.inputs[dict_diesel_map[diesel_technology]],
                    self.index_diesel] = (
                (self.temp_array[self.array_inputs["fuel mass"], index] / 48) # LHV biodiesel 40 MJ/kg
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )

        index = self.get_index_from_array(["ICEV-p", 'HEV-p', 'PHEV'])

        if background_configuration:
            if 'petrol technology' in background_configuration:
                # If a customization dict is passed
                petrol_technology = background_configuration['petrol technology']
            else:
                petrol_technology = 'petrol'
        else:
            petrol_technology = 'petrol'

        dict_petrol_map = {
            'petrol': ('market for petrol, low-sulfur',
                      'Europe without Switzerland',
                      'kilogram',
                      'petrol, low-sulfur'),
            'bioethanol - wheat straw': ('ethanol from wheat straw pellets',
                                          'RER',
                                          'megajoule',
                                          'ethanol from wheat straw pellets'),
            'bioethanol - forest residues': ('Ethanol from forest residues',
                                            'RER',
                                              'megajoule',
                                              'Ethanol from forest residues'),
            'bioethanol - sugarbeet' : ('Ethanol from sugarbeet (purity >95%) production & distribution {RER} | fermentation of sugar beet juice with anaerobic digestion of wastes | Alloc Rec,  U',
                                      'RER',
                                      'megajoule',
                                      'Ethanol from sugarbeet (purity >95%) production & distribution {RER} | fermentation of sugar beet juice with anaerobic digestion of wastes | Alloc Rec,  U') ,
            'bioethanol - maize starch' : ('Ethanol from maize (purity >95%) production & distribution {RER} | fermentation of maize starch | Alloc Rec,  U',
                                          'RER',
                                          'megajoule',
                                          'Ethanol from maize (purity >95%) production & distribution {RER} | fermentation of maize starch | Alloc Rec,  U')
        }

        if petrol_technology == 'petrol':
            self.A[self.inputs[dict_petrol_map[petrol_technology]],
                    self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = (
                self.temp_array[self.array_inputs["fuel mass"], index]
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )
        else:
            # bioethanol
            self.A[self.inputs[dict_petrol_map[petrol_technology]],
                    self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = (
                (self.temp_array[self.array_inputs["fuel mass"], index] / 42.4) # LHV petrol 42.4 MJ/kg
                / self.temp_array[self.array_inputs["range"], index]
                * -1
                )

        # Non-exhaust emissions
        self.A[
            self.inputs[('market for road wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'road wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 1e-08)
        self.A[
            self.inputs[('market for tyre wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'tyre wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 6e-08)
        self.A[
            self.inputs[('market for brake wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'brake wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 5e-09)

        # Infrastructure
        self.A[self.inputs[('market for road', 'GLO', 'meter-year', 'road')], -self.number_of_cars :] = (
            5.37e-7 * self.temp_array[self.array_inputs["driving mass"], :] * -1
        )

        # Exhaust emissions
        # Fuel-based emissions
        self.A[
            self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')], -self.number_of_cars :
        ] = (
            (
                self.temp_array[self.array_inputs["CO2 per kg fuel"], :]
                * self.temp_array[self.array_inputs["fuel mass"], :]
            )
            / self.temp_array[self.array_inputs["range"], :]
            * -1
        )

        # Non-fuel based emissions
        list_direct_emissions = [
            "Benzene direct emissions, urban",
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
            "Lead direct emissions, rural",
        ]

        self.A[self.index_emissions, -self.number_of_cars :] = (
            self.temp_array[[self.array_inputs[l] for l in list_direct_emissions], :]
            * -1
        )

        if cng_technology=='biogas':
            # change fossil emissions to non-fossil, if first iteration
            old_co2_amounts = self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_cng]
            self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_cng] = 0
            self.A[
                self.inputs[('Carbon dioxide, non-fossil', ('air',), 'kilogram')],
                self.index_cng] = old_co2_amounts


        if diesel_technology in ('biodiesel - algae', 'biodiesel - cooking oil'):
            # change fossil emissions to non-fossil, if first iteration
            old_co2_amounts = self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_diesel]
            self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_diesel] = 0
            self.A[
                self.inputs[('Carbon dioxide, non-fossil', ('air',), 'kilogram')],
                self.index_diesel] = old_co2_amounts

        if petrol_technology in ('bioethanol - wheat straw', 'bioethanol - forest residues',
                                 'bioethanol - sugarbeet', 'bioethanol - maize starch'):
            # change fossil emissions to non-fossil, if first iteration
            old_co2_amounts = self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_petrol + self.index_hybrid + self.index_plugin_hybrid]
            self.A[
                self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = 0
            self.A[
                self.inputs[('Carbon dioxide, non-fossil', ('air',), 'kilogram')],
                self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = old_co2_amounts

        # Noise emissions
        self.A[self.index_noise, -self.number_of_cars :] = (
            self.temp_array[
                [self.array_inputs[self.rev_inputs[e][0]] for e in self.index_noise], :
            ]
            * -1
        )
