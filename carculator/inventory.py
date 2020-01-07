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
from .export import ExportInventory

class InventoryCalculation:
    """
    Build and solve the inventory for results characterization and inventory export

    Vehicles to be analyzed can be filtered by passing a `scope` dictionary.
    Some assumptions in the background system can also be adjusted by passing a `background_configuration` dictionary.

    .. code-block:: python

        scope = {
                        'powertrain':['BEV', 'FCEV', 'ICEV-p'],
                    }
        background_configuration = {
                                            'country' : 'DE', # will use the network electricity losses of Germany
                                            'custom electricity mix' : [[1,0,0,0,0,0,0,0,0,0], # in this case, 100% hydropower for the first year
                                                                        [0.5,0.5,0,0,0,0,0,0,0,0]], # in this case, 50% hydro, 50% nuclear for the second year
                                            'hydrogen technology' : 'Electrolysis',
                                            'petrol technology': 'bioethanol - wheat straw',
                                            'battery technology': 'LFP',
                                            'battery origin': 'NO'
                                        }

    The `custom electricity mix` key in the background_configuration dictionary defines an electricity mix to apply,
    under the form of one or several array(s), depending on teh number of years to analyze,
    that should total 1, of which the indices correspond to:

        - [0]: hydro-power
        - [1]: nuclear
        - [2]: natural gas
        - [3]: solar power
        - [4]: wind power
        - [5]: biomass
        - [6]: coal
        - [7]: oil
        - [8]: geothermal
        - [9]: waste incineration

    If none is given, the electricity mix corresponding tot he country specified in `country` will be selected.
    If no country is specified, Europe applies.

    :ivar array: array from the CarModel class
    :vartype array: CarModel.array
    :ivar scope: dictionary that contains filters for narrowing the analysis
    :ivar background_configuration: dictionary that contains choices for background system


    """

    def __init__(self, array, scope = None, background_configuration = None):

        if scope is None:
            scope = {}
            scope['size'] = array.coords['size'].values.tolist()
            scope['powertrain'] = array.coords['powertrain'].values.tolist()
            scope['year'] = array.coords['year'].values.tolist()
        else:
            scope['size'] = scope.get('size', array.coords['size'].values.tolist())
            scope['powertrain'] = scope.get('powertrain', array.coords['powertrain'].values.tolist())
            scope['year'] = scope.get('year', array.coords['year'].values.tolist())

        self.scope = scope

        array = array.sel(powertrain=self.scope["powertrain"], year=self.scope["year"], size=self.scope["size"])

        self.number_of_cars = len(self.scope["year"]) * len(self.scope["size"]) * len(self.scope["powertrain"])
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
        self.background_configuration = background_configuration

    def __getitem__(self, key):
        """
        Make class['foo'] automatically filter for the parameter 'foo'
        Makes the model code much cleaner

        :param key: Parameter name
        :type key: str
        :return: `array` filtered after the parameter selected
        """
        return self.temp_array.sel(parameter=key)

    def get_results_table(self, method, level, split):
        """
        Format an xarray.DataArray array to receive the results.

        :param method: impact assessment method. Only "ReCiPe" method available at the moment.
        :param level: "midpoint" or "endpoint" impact assessment level. Only "midpoint" available at the moment.
        :param split: "components" or "impact categories". Split by impact categories only applicable when "endpoint" level is applied.
        :return: xarrray.DataArray
        """

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
                    len(self.scope['size']),
                    len(self.scope['powertrain']),
                    len(self.scope['year']),
                    len(cat),
                    self.iterations,
                )
            ),
            coords=[
                dict_impact_cat[method][level],
                self.scope['size'],
                self.scope['powertrain'],
                self.scope['year'],
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
        """
        Return list of indices to split the results into categories.

        :return: list of indices
        :rtype: list
        """
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

    def calculate_impacts(self, method="recipe", level="midpoint", split="components"):
        """
        Solve the inventory, fill in the results array, return the results array.

        :param method: Impact assessment method. Only "recipe" available at the moment.
        :param level: Impact assessment level ("midpoint" or "endpoint"). Only "midpoint" available at the moment.
        :param split: Splitting mode ("components" or "impact categories"). Only "components" available at the moment.
        :return: an array with characterized results.
        :rtype: xarray.DataArray

        """


        # Load the B matrix
        self.B = self.get_B_matrix(method, level)

        # Prepare an array to store the results
        results = self.get_results_table(method, level, split)

        # Fill in the A matrix with car parameters
        self.set_inputs_in_A_matrix(self.array.values)

        for pt in self.scope['powertrain']:
            for y in self.scope['year']:
                for s in self.scope['size']:
                    # Retrieve the index of a given car in the matrix A
                    car = self.inputs[('Passenger car, '+ pt + ', ' + s + ', ' + str(y),"GLO", "kilometer", "transport, passenger car, EURO6")]
                    # Set the demand vector with zeros and a 1 corresponding to the car position in the vector
                    f = np.zeros((np.shape(self.A)[0],np.shape(self.A)[1] ))
                    f[:,car] = 1

                    X = np.linalg.solve(self.A, f)
                    X_reshaped = np.zeros((X.shape[0],X.shape[1] , X.shape[1]))
                    indices = np.arange(X.shape[1])
                    X_reshaped[:, indices, indices] = X

                    C = np.dot(X_reshaped, self.B.T).T

                    # Iterate through the results array to fill it
                    results.loc[
                        dict(impact=self.list_cat, year=y, size=s, powertrain=pt)] = \
                        C[:, self.split_indices].sum(axis=2)

        return results

    def get_A_matrix(self):
        """
        Load the A matrix. The A matrix contains exchanges of products (rows) between activities (columns).

        :return: A matrix with three dimensions of shape (number of values, number of products, number of activities).
        :rtype: numpy.ndarray

        """
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
        # Resize the matrix to fit the number of iterations in `array`
        new_A = np.resize(new_A,
                           (self.array.shape[1],
                               new_A.shape[0],
                            new_A.shape[1]
                            )
                           )
        return new_A

    def get_B_matrix(self, method, level):
        """
        Load the B matrix. The B matrix contains impact assessment figures for a give impact assessment method,
        per unit of activity. Its length column-wise equals the length of the A matrix row-wise.
        Its length row-wise equals the number of impact assessment methods.

        :param method: only "recipe" available at the moment.
        :param level: only "midpoint" available at the moment.
        :return: an array with impact values per unit of activity for each method.
        :rtype: numpy.ndarray

        """
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
        """
        Load a dictionary with tuple ("name of activity", "location", "unit", "reference product") as key, row/column
        indices as values.

        :return: dictionary with `label:index` pairs.
        :rtype: dict

        """
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

        for pt in self.scope["powertrain"]:
            for s in self.scope["size"]:
                for y in self.scope["year"]:
                    maximum = csv_dict[max(csv_dict, key=csv_dict.get)]
                    name = "Passenger car, " + pt + ", " + s + ", " + str(y)
                    csv_dict[(name, "GLO", "kilometer", "transport, passenger car, EURO6")] = maximum + 1

        return csv_dict

    def get_dict_impact_categories(self):
        """
        Load a dictionary with available impact assessment methods as keys, and assessment level and categories as values.

        ..code-block:: python

            {'recipe': {'midpoint': ['freshwater ecotoxicity',
                                   'human toxicity',
                                   'marine ecotoxicity',
                                   'terrestrial ecotoxicity',
                                   'metal depletion',
                                   'agricultural land occupation',
                                   'climate change',
                                   'fossil depletion',
                                   'freshwater eutrophication',
                                   'ionising radiation',
                                   'marine eutrophication',
                                   'natural land transformation',
                                   'ozone depletion',
                                   'particulate matter formation',
                                   'photochemical oxidant formation',
                                   'terrestrial acidification',
                                   'urban land occupation',
                                   'water depletion',
                                   'human noise']
                       }
           }

        :return: dictionary
        :rtype: dict
        """
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
        """
        Reverse the self.inputs dictionary.

        :return: reversed dictionary
        :rtype: dict
        """
        return {v: k for k, v in self.inputs.items()}

    def get_index_from_array(self, items_to_look_for):
        """
        Return list of row/column indices of self.array of labels that contain the string defined in `items_to_look_for`.

        :param items_to_look_for: string to search for
        :return: list
        """
        return [
            self.inputs[c] - (len(self.inputs) - self.number_of_cars)
            for c in self.inputs
            if any(ele in c[0] for ele in items_to_look_for)
            and (self.inputs[c] - (len(self.inputs) - self.number_of_cars))>=0
        ]

    def get_index_of_flows(self, items_to_look_for, search_by="name"):
        """
        Return list of row/column indices of self.A of labels that contain the string defined in `items_to_look_for`.

        :param items_to_look_for: string
        :param search_by: "name" or "compartment" (for elementary flows)
        :return: list of row/column indices
        :rtype: list
        """
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

    def export_lci(self, presamples = True):
        """
        Export the inventory as a dictionary. Also return a list of arrays that contain pre-sampled random values if
        :meth:`stochastic` of :class:`CarModel` class has been called.

        :param presamples: boolean.
        :return: inventory, and optionally, list of arrays containing pre-sampled values.
        :rtype: list
        """
        self.set_inputs_in_A_matrix(self.array.values)
        if presamples == True:
            lci, array = ExportInventory(self.A, self.rev_inputs).write_lci(presamples)
            return (lci, array)
        else:
            lci = ExportInventory(self.A, self.rev_inputs).write_lci(presamples)
            return lci

    def export_lci_to_bw(self, presamples = True):
        """
        Export the inventory as a `brightway2` bw2io.importers.base_lci.LCIImporter object
        with the inventory in the `data` attribute.

        .. code-block:: python

            # get the invenotry
            i, _ = ic.export_lci_to_bw()

            # import it in a Brightway2 project
            i.match_database('ecoinvent 3.6 cutoff', fields=('name', 'unit', 'location', 'reference product'))
            i.match_database("biosphere3", fields=('name', 'unit', 'categories'))
            i.match_database(fields=('name', 'unit', 'location', 'reference product'))
            i.match_database(fields=('name', 'unit', 'categories'))

            # Create an additional biosphere database for the few flows that do not
            # exist in "biosphere3"
            i.create_new_biosphere("additional_biosphere", relink=True)

            # Check if all exchanges link
            i.statistics()

            # Register the database
            i.write_database()

        :return: LCIImport object that can be directly registered in a `brightway2` project.
        :rtype: bw2io.importers.base_lci.LCIImporter
        """

        self.set_inputs_in_A_matrix(self.array.values)
        if presamples == True:
            lci, array = ExportInventory(self.A, self.rev_inputs).write_lci_to_bw(presamples)
            return (lci, array)
        else:
            lci = ExportInventory(self.A, self.rev_inputs).write_lci_to_bw(presamples)
            return lci

    def export_lci_to_excel(self, directory=None):
        """
        Export the inventory as an Excel file. Also return the file path where the file is stored.

        :param directory: directory where to save the file.
        :type directory: str
        :return: file path where the file is stored.
        :rtype: str
        """

        self.set_inputs_in_A_matrix(self.array.values)
        fp = ExportInventory(self.A, self.rev_inputs).write_lci_to_excel(directory=directory)
        return fp

    def set_inputs_in_A_matrix(self, array):
        """
        Fill-in the A matrix. Does not return anything. Modifies in place.
        Shape of the A matrix (values, products, activities).

        :param array: :attr:`array` from :class:`CarModel` class
        """

        # Glider
        self.A[:,
            self.inputs[('market for glider, passenger car', 'GLO','kilogram','glider, passenger car')],
            -self.number_of_cars :
        ] = (
            (array[self.array_inputs["glider base mass"], :] * (1-array[self.array_inputs["lightweighting"], :]))
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[('Glider lightweighting', 'GLO', 'kilogram', 'Glider lightweighting')], -self.number_of_cars :
        ] = (
            (array[self.array_inputs["lightweighting"], :] * array[self.array_inputs["glider base mass"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[('maintenance, passenger car', 'RER', 'unit', 'passenger car maintenance')], -self.number_of_cars :
        ] = (array[self.array_inputs["curb mass"], :] / 1600 / 150000 * -1)

        # Glider EoL
        self.A[:,
            self.inputs[
                ('market for manual dismantling of used electric passenger car',
                 'GLO',
                 'unit',
                 'manual dismantling of used electric passenger car')
            ],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["curb mass"], :]
            * (1 - array[self.array_inputs["combustion power share"], :])
            / 1180
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[
                ('market for manual dismantling of used passenger car with internal combustion engine',
                 'GLO',
                 'unit',
                 'manual dismantling of used passenger car with internal combustion engine')
            ],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["curb mass"], :]
            * array[self.array_inputs["combustion power share"], :]
            / 1600
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Powertrain components
        self.A[:,
            self.inputs[('market for charger, electric passenger car',
              'GLO',
              'kilogram',
              'charger, electric passenger car')],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["charger mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[('market for converter, for electric passenger car',
              'GLO',
              'kilogram',
              'converter, for electric passenger car')],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["converter mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[('market for electric motor, electric passenger car',
              'GLO',
              'kilogram',
              'electric motor, electric passenger car')],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["electric engine mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[('market for inverter, for electric passenger car',
              'GLO',
              'kilogram',
              'inverter, for electric passenger car')],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["inverter mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,
            self.inputs[
                ('market for power distribution unit, for electric passenger car',
                 'GLO',
                 'kilogram',
                 'power distribution unit, for electric passenger car')
            ],
            -self.number_of_cars :
        ] = (
            array[self.array_inputs["power distribution unit mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
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

        self.A[:,
            self.inputs[
                ('market for used powertrain from electric passenger car, manual dismantling',
                 'GLO',
                 'kilogram',
                 'used powertrain from electric passenger car, manual dismantling')
            ],
            -self.number_of_cars :
        ] = (
            array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)
            / array[self.array_inputs["lifetime kilometers"], :]
        )

        self.A[:,
            self.inputs[
                ('market for internal combustion engine, passenger car',
                 'GLO',
                 'kilogram',
                 'internal combustion engine, for passenger car')
            ],
            -self.number_of_cars :
        ] = (
            (
                array[
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
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,self.inputs[('Ancillary BoP', 'GLO', 'kilogram', 'Ancillary BoP')], -self.number_of_cars :] = (
            array[self.array_inputs["fuel cell ancillary BoP mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,self.inputs[('Essential BoP', 'GLO', 'kilogram', 'Essential BoP')], -self.number_of_cars :] = (
            array[self.array_inputs["fuel cell essential BoP mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[:,self.inputs[('Stack 2020', 'GLO', 'kilowatt', 'Stack 2020')], -self.number_of_cars :] = (
            array[self.array_inputs["fuel cell stack mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
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

        # Energy storage
        if self.background_configuration:
            if "battery technology" in self.background_configuration:
                battery_tech = self.background_configuration["battery technology"]
            else:
                battery_tech = "NMC"

            if "battery origin" in self.background_configuration:
                battery_origin = self.background_configuration['battery origin']
                losses_to_medium = float(self.bs.losses[battery_origin]['MV'])
                mix = self.bs.electricity_mix.sel(country=battery_origin, value=0).interp(year=self.scope["year"]).values
            else:
                # If not specified, origin set to China
                battery_origin = "CN"
                losses_to_medium = float(self.bs.losses[battery_origin]['MV'])
                mix = self.bs.electricity_mix.sel(country=battery_origin, value=0).interp(year=self.scope["year"]).values
        else:
            battery_tech = "NMC"
            battery_origin = "CN"
            losses_to_medium = float(self.bs.losses[battery_origin]['MV'])
            mix = self.bs.electricity_mix.sel(country=battery_origin, value=0).interp(year=self.scope["year"]).values

        if battery_tech == "NMC":
            # Use the NMC inventory of Schmidt et al. 2019
            self.A[:,self.inputs[('Battery BoP', 'GLO', 'kilogram', 'Battery BoP')], -self.number_of_cars :] = (
                (
                    array[self.array_inputs["battery BoP mass"], :]
                    * (
                        1
                        + array[
                            self.array_inputs["battery lifetime replacements"], :
                        ]
                    )
                )
                / array[self.array_inputs["lifetime kilometers"], :]
                * -1
            )

            self.A[:,self.inputs[('Battery cell', 'GLO', 'kilogram', 'Battery cell')], -self.number_of_cars :] = (
                (
                    array[self.array_inputs["battery cell mass"], :]
                    * (
                        1
                        + array[
                            self.array_inputs["fuel cell lifetime replacements"], :
                        ]
                    )
                )
                / array[self.array_inputs["lifetime kilometers"], :]
                * -1
            )

            # Set an input of electricity, given the country of manufacture
            self.A[:,self.inputs[('market group for electricity, medium voltage', 'GLO',
                                'kilowatt hour', 'electricity, medium voltage')],
                   self.inputs[('Battery cell', 'GLO', 'kilogram', 'Battery cell')]] = 0

            self.A[:,[self.inputs[dict_map[t]] for t in dict_map],
                   self.inputs[('Battery cell', 'GLO', 'kilogram', 'Battery cell')]] =\
                    (np.outer(mix[0],
                              (array[self.array_inputs["battery cell production heat"], :].max(axis=1)
                               * -1)) * losses_to_medium).T

        if battery_tech == "LFP":
            self.A[:,self.inputs[('Li-ion (LFP)', 'JP', 'kilowatt hour', 'Li-ion (LFP)')], -self.number_of_cars :] = (
                (
                    array[self.array_inputs["electric energy stored"], :]
                    * (
                        1
                        + array[
                            self.array_inputs["battery lifetime replacements"], :
                        ]
                    )
                )
                / array[self.array_inputs["lifetime kilometers"], :]
                * -1
            )

            # Set an input of electricity, given the country of manufacture
            old_amount = self.A[:,self.inputs[('market for electricity, medium voltage', 'JP',
                                'kilowatt hour', 'electricity, medium voltage')],
                   self.inputs[('Li-ion (LFP)', 'JP', 'kilowatt hour', 'Li-ion (LFP)')]]
            self.A[:,self.inputs[('market for electricity, medium voltage', 'JP',
                                'kilowatt hour', 'electricity, medium voltage')],
                   self.inputs[('Li-ion (LFP)', 'JP', 'kilowatt hour', 'Li-ion (LFP)')]] = 0

            self.A[:,[self.inputs[dict_map[t]] for t in dict_map],
                   self.inputs[('Li-ion (LFP)', 'JP', 'kilowatt hour', 'Li-ion (LFP)')]] =\
                    (np.outer(mix[0], old_amount) * losses_to_medium).T

        if battery_tech == "NCA":
            self.A[:,self.inputs[('Li-ion (NCA)', 'JP', 'kilowatt hour', 'Li-ion (NCA)')], -self.number_of_cars :] = (
                (
                    array[self.array_inputs["electric energy stored"], :]
                    * (
                        1
                        + array[
                            self.array_inputs["battery lifetime replacements"], :
                        ]
                    )
                )
                / array[self.array_inputs["lifetime kilometers"], :]
                * -1
            )

            # Set an input of electricity, given the country of manufacture
            old_amount = self.A[:,self.inputs[('market for electricity, medium voltage', 'JP',
                                'kilowatt hour', 'electricity, medium voltage')],
                   self.inputs[('Li-ion (NCA)', 'JP', 'kilowatt hour', 'Li-ion (NCA)')]]
            self.A[:,self.inputs[('market for electricity, medium voltage', 'JP',
                                'kilowatt hour', 'electricity, medium voltage')],
                   self.inputs[('Li-ion (NCA)', 'JP', 'kilowatt hour', 'Li-ion (NCA)')]] = 0

            self.A[:,[self.inputs[dict_map[t]] for t in dict_map],
                   self.inputs[('Li-ion (NCA)', 'JP', 'kilowatt hour', 'Li-ion (NCA)')]] =\
                    (np.outer(mix[0], old_amount) * losses_to_medium).T


        index_A = [
            self.inputs[c]
            for c in self.inputs
            if any(ele in c[0] for ele in ["ICEV-d", "ICEV-p", "HEV-p", "PHEV"])
        ]
        index = self.get_index_from_array(["ICEV-d", "ICEV-p", "HEV-p", "PHEV"])

        self.A[:,
            self.inputs[('polyethylene production, high density, granulate',
              'RER',
              'kilogram',
              'polyethylene, high density, granulate')],
            index_A
        ] = (
            array[self.array_inputs["fuel tank mass"],:, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        index = self.get_index_from_array(["ICEV-g"])
        self.A[:,
            self.inputs[
                ('glass fibre reinforced plastic production, polyamide, injection moulded',
                 'RER',
                 'kilogram',
                 'glass fibre reinforced plastic, polyamide, injection moulded')
            ],
            self.index_cng,
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        index = self.get_index_from_array(["FCEV"])
        self.A[:,
            self.inputs[('Fuel tank, compressed hydrogen gas, 700bar',
              'GLO',
              'kilogram',
              'Fuel tank, compressed hydrogen gas, 700bar')],
            self.index_fuel_cell,
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T


        if self.background_configuration:
            # If a customization dict is passed
            if 'country' in self.background_configuration:
                # If a country is specified
                country = self.background_configuration['country']
                try:
                    losses_to_low = float(self.bs.losses[country]['LV'])
                except KeyError:
                    # If losses for the country are not found, assume 15%
                    losses_to_low = 1.15

                if 'custom electricity mix' in self.background_configuration:
                    # If a special electricity mix is specified, we use it
                    mix = self.background_configuration['custom electricity mix']
                else:
                    mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.scope["year"]).values
            else:
                country = 'RER'
                losses_to_low = float(self.bs.losses["RER"]['LV'])
                mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.scope["year"]).values
        else:
            country = 'RER'
            losses_to_low = float(self.bs.losses["RER"]['LV'])
            mix = self.bs.electricity_mix.sel(country=country, value=0).interp(year=self.scope["year"]).values

        for y in self.scope["year"]:
            index = self.get_index_from_array([str(y)])

            self.A[np.ix_(np.arange(self.iterations),[self.inputs[dict_map[t]] for t in dict_map],
                          [self.inputs[i] for i in self.inputs
                           if str(y) in i[0]
                           and "Passenger" in i[0]
                           ])] = \
                (np.outer(mix[self.scope["year"].index(y)],
                         array[self.array_inputs["electricity consumption"], :, index]
                         ) * -1 * losses_to_low).reshape(self.iterations, len(mix[self.scope["year"].index(y)]), len(index))

        if "FCEV" in self.scope['powertrain']:

            index = self.get_index_from_array(["FCEV"])

            if self.background_configuration:
                if 'hydrogen technology' in self.background_configuration:
                    # If a customization dict is passed
                    hydro_technology = self.background_configuration['hydrogen technology']
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


            self.A[:,self.inputs[dict_h_map[hydro_technology]],
                    self.index_fuel_cell] = (
                array[self.array_inputs["fuel mass"], :, index]
                / array[self.array_inputs["range"], :, index]
                * -1
                ).T

            # If hydrolysis is chosen, adjust the electricity mix

            if hydro_technology == 'Electrolysis':

                # Zero out initial electricity provider
                old_amount = self.A[:,self.inputs[('market group for electricity, medium voltage',
                  'Europe without Switzerland',
                  'kilowatt hour',
                  'electricity, medium voltage')],
                   self.inputs[dict_h_map[hydro_technology]]]

                self.A[:, self.inputs[('market group for electricity, medium voltage',
                                       'Europe without Switzerland',
                                       'kilowatt hour',
                                       'electricity, medium voltage')],
                self.inputs[dict_h_map[hydro_technology]]] = 0

                # TODO: differentiate hydrogen production in time

                for y in self.scope["year"]:
                    index = self.get_index_from_array([str(y)])

                    self.A[np.ix_(np.arange(self.iterations), [self.inputs[dict_map[t]] for t in dict_map],
                                  [self.inputs[i] for i in self.inputs
                                   if str(y) in i[0]
                                   and "Passenger" in i[0]
                                   and "FCEV" in i[0]
                                   ])] = \
                        (np.outer(mix[self.scope["year"].index(y)], old_amount) * losses_to_medium)\
                            .reshape(self.iterations, 10, -1)

        if 'ICEV-g' in self.scope['powertrain']:
            index = self.get_index_from_array(["ICEV-g"])

            if self.background_configuration:
                if 'cng technology' in self.background_configuration:
                    # If a customization dict is passed
                    cng_technology = self.background_configuration['cng technology']
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
                          'kilogram',
                          'biogas upgrading - sewage sludge - amine scrubbing - best')
            }

            if cng_technology == 'cng':
                self.A[:,self.inputs[dict_cng_map[cng_technology]],
                        self.index_cng] = (
                    array[self.array_inputs["fuel mass"], :, index]
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T
            else:
                # biogas
                self.A[:,self.inputs[dict_cng_map[cng_technology]],
                        self.index_cng] = (
                    (array[self.array_inputs["fuel mass"], :, index])
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T

        if "ICEV-d" in self.scope['powertrain']:

            index = self.get_index_from_array(["ICEV-d"])
            if self.background_configuration:
                if 'diesel technology' in self.background_configuration:
                    # If a customization dict is passed
                    diesel_technology = self.background_configuration['diesel technology']
                else:
                    diesel_technology = 'diesel'
            else:
                diesel_technology = 'diesel'

            dict_diesel_map = {
                'diesel': ('market for diesel', 'Europe without Switzerland', 'kilogram', 'diesel'),
                'biodiesel - algae': ('Biodiesel from algae','RER','kilogram','Biodiesel from algae'),
                'biodiesel - cooking oil': ('Biodiesel from cooking oil','RER','kilogram','Biodiesel from cooking oil')
            }

            if diesel_technology == 'diesel':
                self.A[:,self.inputs[dict_diesel_map[diesel_technology]],
                        self.index_diesel] = (
                    array[self.array_inputs["fuel mass"], :, index]
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T
            else:
                # biodiesel
                self.A[:,self.inputs[dict_diesel_map[diesel_technology]],
                        self.index_diesel] = (
                    (array[self.array_inputs["fuel mass"], :, index] / 48) # LHV biodiesel 40 MJ/kg
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T

        if [i for i in self.scope['powertrain'] if i in ["ICEV-p", 'HEV-p', 'PHEV']]:
            index = self.get_index_from_array(["ICEV-p", 'HEV-p', 'PHEV'])

            if self.background_configuration:
                if 'petrol technology' in self.background_configuration:
                    # If a customization dict is passed
                    petrol_technology = self.background_configuration['petrol technology']
                else:
                    petrol_technology = 'petrol'
            else:
                petrol_technology = 'petrol'

            dict_petrol_map = {
                'petrol': ('market for petrol, low-sulfur',
                          'Europe without Switzerland',
                          'kilogram',
                          'petrol, low-sulfur'),
                'bioethanol - wheat straw': ('Ethanol from wheat straw pellets','RER','kilogram','Ethanol from wheat straw pellets'),
                'bioethanol - forest residues': ('Ethanol from forest residues','RER','kilogram','Ethanol from forest residues'),
                'bioethanol - sugarbeet' : ('Ethanol from sugarbeet','RER','kilogram','Ethanol from sugarbeet') ,
                'bioethanol - maize starch' : ('Ethanol from maize starch','RER','kilogram','Ethanol from maize starch')
            }

            if petrol_technology == 'petrol':
                self.A[:,self.inputs[dict_petrol_map[petrol_technology]],
                        self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = (
                    array[self.array_inputs["fuel mass"], :, index]
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T
            else:
                # bioethanol
                self.A[:,self.inputs[dict_petrol_map[petrol_technology]],
                        self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = (
                    (array[self.array_inputs["fuel mass"], :, index] / 42.4) # LHV petrol 42.4 MJ/kg
                    / array[self.array_inputs["range"], :, index]
                    * -1
                    ).T

        # Non-exhaust emissions
        self.A[:,
            self.inputs[('market for road wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'road wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 1e-08)
        self.A[:,
            self.inputs[('market for tyre wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'tyre wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 6e-08)
        self.A[:,
            self.inputs[('market for brake wear emissions, passenger car',
                          'GLO',
                          'kilogram',
                          'brake wear emissions, passenger car')],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 5e-09)

        # Infrastructure
        self.A[:,self.inputs[('market for road', 'GLO', 'meter-year', 'road')], -self.number_of_cars :] = (
            5.37e-7 * array[self.array_inputs["driving mass"], :] * -1
        )

        # Exhaust emissions
        # Fuel-based emissions
        self.A[:,
            self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')], -self.number_of_cars :
        ] = (
            (
                array[self.array_inputs["CO2 per kg fuel"], :]
                * array[self.array_inputs["fuel mass"], :]
            )
            / array[self.array_inputs["range"], :]
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

        self.A[:,self.index_emissions, -self.number_of_cars :] = (
            array[[self.array_inputs[l] for l in list_direct_emissions]]
            * -1
        ).transpose([1,0,2])

        if 'cng_technology' in locals():
            if cng_technology == 'biogas':
                # change fossil emissions to non-fossil, if first iteration
                self.A[:,
                    self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                    self.index_cng] = 0

        if 'diesel_technology' in locals():
            if diesel_technology in ('biodiesel - algae', 'biodiesel - cooking oil'):
                # change fossil emissions to non-fossil, if first iteration

                self.A[:,
                    self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                    self.index_diesel] = 0

        if 'petrol_technology' in locals():
            if petrol_technology in ('bioethanol - wheat straw', 'bioethanol - forest residues',
                                 'bioethanol - sugarbeet', 'bioethanol - maize starch'):
            # change fossil emissions to non-fossil, if first iteration
                self.A[:,
                    self.inputs[('Carbon dioxide, fossil', ('air',), 'kilogram')],
                    self.index_petrol + self.index_hybrid + self.index_plugin_hybrid] = 0

        # Noise emissions
        self.A[:,self.index_noise, -self.number_of_cars :] = (
            array[
                [self.array_inputs[self.rev_inputs[e][0]] for e in self.index_noise], :
            ]
            * -1
        ).transpose([1,0,2])

