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
from .export import ExportInventory
from scipy import sparse
import numexpr as ne
import itertools


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

        InventoryCalculation(CarModel.array,
                            background_configuration=background_configuration,
                            scope=scope,
                            scenario="RCP26")

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

    If none is given, the electricity mix corresponding to the country specified in `country` will be selected.
    If no country is specified, Europe applies.

    :ivar array: array from the CarModel class
    :vartype array: CarModel.array
    :ivar scope: dictionary that contains filters for narrowing the analysis
    :ivar background_configuration: dictionary that contains choices for background system
    :ivar scenario: REMIND energy scenario to use ("BAU": business-as-usual or "RCP26": limits radiative forcing to 2.6 W/m^2.).
                    "BAU" selected by default.

    .. code-block:: python




    """

    def __init__(
        self, array, scope=None, background_configuration=None, scenario="BAU"
    ):

        if scope is None:
            scope = {}
            scope["size"] = array.coords["size"].values.tolist()
            scope["powertrain"] = array.coords["powertrain"].values.tolist()
            scope["year"] = array.coords["year"].values.tolist()
        else:
            scope["size"] = scope.get("size", array.coords["size"].values.tolist())
            scope["powertrain"] = scope.get(
                "powertrain", array.coords["powertrain"].values.tolist()
            )
            scope["year"] = scope.get("year", array.coords["year"].values.tolist())

        self.scope = scope
        self.scenario = scenario

        if background_configuration is None:
            self.background_configuration = {"country": "RER"}
        else:
            self.background_configuration = background_configuration

        if "country" not in self.background_configuration:
            self.background_configuration["country"] = "RER"

        if "battery technology" not in self.background_configuration:
            self.background_configuration["battery technology"] = "NMC"

        if "battery origin" not in self.background_configuration:
            self.background_configuration["battery origin"] = "CN"

        array = array.sel(
            powertrain=self.scope["powertrain"],
            year=self.scope["year"],
            size=self.scope["size"],
        )
        self.array = array.stack(desired=["size", "powertrain", "year"])
        self.iterations = len(array.value.values)

        self.number_of_cars = (
            len(self.scope["size"])
            * len(self.scope["powertrain"])
            * len(self.scope["year"])
        )

        self.array_inputs = {
            x: i for i, x in enumerate(list(self.array.parameter.values), 0)
        }
        self.array_powertrains = {
            x: i for i, x in enumerate(list(self.array.powertrain.values), 0)
        }

        self.A = self.get_A_matrix()
        self.inputs = self.get_dict_input()
        self.add_additional_activities()

        self.rev_inputs = self.get_rev_dict_input()

        self.index_cng = [self.inputs[i] for i in self.inputs if "ICEV-g" in i[0]]
        self.index_combustion_wo_cng = [
            self.inputs[i]
            for i in self.inputs
            if any(
                ele in i[0]
                for ele in ["ICEV-p", "HEV-p", "PHEV-p", "ICEV-d", "PHEV-d", "HEV-d"]
            )
        ]
        self.index_diesel = [self.inputs[i] for i in self.inputs if "ICEV-d" in i[0]]
        self.index_all_petrol = [
            self.inputs[i]
            for i in self.inputs
            if any(ele in i[0] for ele in ["ICEV-p", "HEV-p", "PHEV-p"])
        ]
        self.index_petrol = [self.inputs[i] for i in self.inputs if "ICEV-p" in i[0]]
        self.index_hybrid = [
            self.inputs[i]
            for i in self.inputs
            if any(ele in i[0] for ele in ["HEV-p", "HEV-d"])
        ]
        self.index_plugin_hybrid = [
            self.inputs[i] for i in self.inputs if "PHEV" in i[0]
        ]
        self.index_fuel_cell = [self.inputs[i] for i in self.inputs if "FCEV" in i[0]]
        self.index_emissions = [
            self.inputs[i]
            for i in self.inputs
            if "air" in i[1][0]
            and len(i[1]) > 1
            and i[0]
            not in [
                "Carbon dioxide, fossil",
                "Carbon monoxide, non-fossil",
                "Methane, non-fossil",
                "Particulates, > 10 um",
            ]
        ]

        self.map_non_fuel_emissions = {
            (
                "Methane, fossil",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Methane direct emissions, suburban",
            (
                "Methane, fossil",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Methane direct emissions, rural",
            (
                "Lead",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Lead direct emissions, suburban",
            (
                "Ammonia",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Ammonia direct emissions, suburban",
            (
                "NMVOC, non-methane volatile organic compounds, unspecified origin",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "NMVOC direct emissions, urban",
            (
                "PAH, polycyclic aromatic hydrocarbons",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Hydrocarbons direct emissions, urban",
            (
                "Dinitrogen monoxide",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Dinitrogen oxide direct emissions, rural",
            (
                "Nitrogen oxides",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Nitrogen oxides direct emissions, urban",
            (
                "Ammonia",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Ammonia direct emissions, urban",
            (
                "Particulates, < 2.5 um",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Particulate matters direct emissions, suburban",
            (
                "Carbon monoxide, fossil",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Carbon monoxide direct emissions, urban",
            (
                "Nitrogen oxides",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Nitrogen oxides direct emissions, rural",
            (
                "NMVOC, non-methane volatile organic compounds, unspecified origin",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "NMVOC direct emissions, suburban",
            (
                "Benzene",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Benzene direct emissions, suburban",
            (
                "Ammonia",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Ammonia direct emissions, rural",
            (
                "Sulfur dioxide",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Sulfur dioxide direct emissions, rural",
            (
                "NMVOC, non-methane volatile organic compounds, unspecified origin",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "NMVOC direct emissions, rural",
            (
                "Particulates, < 2.5 um",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Particulate matters direct emissions, urban",
            (
                "Sulfur dioxide",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Sulfur dioxide direct emissions, urban",
            (
                "Dinitrogen monoxide",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Dinitrogen oxide direct emissions, suburban",
            (
                "Carbon monoxide, fossil",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Carbon monoxide direct emissions, rural",
            (
                "Methane, fossil",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Methane direct emissions, urban",
            (
                "Carbon monoxide, fossil",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Carbon monoxide direct emissions, suburban",
            (
                "Lead",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Lead direct emissions, urban",
            (
                "Particulates, < 2.5 um",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Particulate matters direct emissions, rural",
            (
                "Sulfur dioxide",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Sulfur dioxide direct emissions, suburban",
            (
                "Benzene",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Benzene direct emissions, rural",
            (
                "Nitrogen oxides",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Nitrogen oxides direct emissions, suburban",
            (
                "Lead",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Lead direct emissions, rural",
            (
                "Benzene",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Benzene direct emissions, urban",
            (
                "PAH, polycyclic aromatic hydrocarbons",
                ("air", "low population density, long-term"),
                "kilogram",
            ): "Hydrocarbons direct emissions, rural",
            (
                "PAH, polycyclic aromatic hydrocarbons",
                ("air", "non-urban air or from high stacks"),
                "kilogram",
            ): "Hydrocarbons direct emissions, suburban",
            (
                "Dinitrogen monoxide",
                ("air", "urban air close to ground"),
                "kilogram",
            ): "Dinitrogen oxide direct emissions, urban",
        }

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

    def get_results_table(self, method, level, split, sensitivity=False):
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
                "EoL",
                "powertrain",
                "energy storage",
                "road",
            ]

        dict_impact_cat = self.get_dict_impact_categories()

        if sensitivity == False:

            response = xr.DataArray(
                np.zeros(
                    (
                        self.B.shape[1],
                        len(self.scope["size"]),
                        len(self.scope["powertrain"]),
                        len(self.scope["year"]),
                        len(cat),
                        self.iterations,
                    )
                ),
                coords=[
                    dict_impact_cat[method][level],
                    self.scope["size"],
                    self.scope["powertrain"],
                    self.scope["year"],
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

        else:
            params = ["reference"]
            params.extend([a for a in self.array_inputs])
            response = xr.DataArray(
                np.zeros(
                    (
                        self.B.shape[1],
                        len(self.scope["size"]),
                        len(self.scope["powertrain"]),
                        len(self.scope["year"]),
                        self.iterations,
                    )
                ),
                coords=[
                    dict_impact_cat[method][level],
                    self.scope["size"],
                    self.scope["powertrain"],
                    self.scope["year"],
                    params,
                ],
                dims=["impact_category", "size", "powertrain", "year", "parameter"],
            )

        return response

    def get_split_indices(self):
        """
        Return list of indices to split the results into categories.

        :return: list of indices
        :rtype: list
        """
        filename = "dict_split.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError("The dictionary of splits could not be found.")

        with open(filepath) as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        (_, _, *header), *data = csv_list

        csv_dict = {}
        for row in data:
            key, sub_key, *values = row

            if key in csv_dict:
                if sub_key in csv_dict[key]:
                    csv_dict[key][sub_key].append(
                        {"search by": values[0], "search for": values[1]}
                    )
                else:
                    csv_dict[key][sub_key] = [
                        {"search by": values[0], "search for": values[1]}
                    ]
            else:
                csv_dict[key] = {
                    sub_key: [{"search by": values[0], "search for": values[1]}]
                }

        flatten = itertools.chain.from_iterable

        d = {}
        l = []

        for cat in csv_dict["components"]:
            d[cat] = list(
                flatten(
                    [
                        self.get_index_of_flows([l["search for"]], l["search by"])
                        for l in csv_dict["components"][cat]
                    ]
                )
            )
            l.append(d[cat])

        list_ind = [d[x] for x in d]
        maxLen = max(map(len, list_ind))
        for row in list_ind:
            while len(row) < maxLen:
                row.extend([len(self.inputs) - 1])
        return list(d.keys()), list_ind

    def calculate_impacts(
        self, method="recipe", level="midpoint", split="components", sensitivity=False
    ):

        # Load the B matrix
        self.B = self.get_B_matrix()

        # Prepare an array to store the results
        results = self.get_results_table(method, level, split, sensitivity=sensitivity)

        # Fill in the A matrix with car parameters
        self.set_inputs_in_A_matrix(self.array.values)

        # Collect indices of activities contributing to the first level
        arr = self.A[0, : -self.number_of_cars, -self.number_of_cars :].sum(axis=1)
        ind = np.nonzero(arr)[0]

        new_arr = np.float16(
            np.zeros((self.A.shape[1], self.B.shape[1], len(self.scope["year"])))
        )

        f = np.float16(np.zeros((np.shape(self.A)[1])))

        for y in self.scope["year"]:
            B = self.B.interp(year=y, kwargs={"fill_value": "extrapolate"}).values

            for a in ind:
                f[:] = 0
                f[a] = 1
                X = np.float16(sparse.linalg.spsolve(self.A[0], f.T))
                C = X * B
                new_arr[a, :, self.scope["year"].index(y)] = C.sum(axis=1)

        new_arr = new_arr.T.reshape(
            len(self.scope["year"]), B.shape[0], 1, 1, self.A.shape[-1]
        )

        a = np.float16(self.A[:, :, -self.number_of_cars :].transpose(0, 2, 1))

        arr = np.float16(ne.evaluate("a * new_arr * -1"))

        arr = arr.transpose(1, 3, 0, 4, 2)
        arr = arr[:, :, :, self.split_indices, :].sum(axis=4)

        for y in range(0, len(self.scope["year"])):
            results[:, :, :, y, :, :] = arr[
                :, y :: len(self.scope["year"]), y, :, :
            ].reshape(
                (
                    B.shape[0],
                    len(self.scope["size"]),
                    len(self.scope["powertrain"]),
                    len(results.impact.values),
                    self.iterations,
                )
            )

        return results.astype("float16")

    def add_additional_activities(self):
        # Add as many rows and columns as cars to consider
        # Also add additional columns and rows for electricity markets
        # for fuel preparation and energy battery production

        maximum = max(self.inputs.values())

        for y in self.scope["year"]:
            maximum += 1
            self.inputs[
                (
                    "electricity market for fuel preparation, " + str(y),
                    self.background_configuration["country"],
                    "kilowatt hour",
                    "electricity, low voltage",
                )
            ] = maximum
            maximum += 1
            self.inputs[
                (
                    "electricity market for energy storage production, " + str(y),
                    self.background_configuration["battery origin"],
                    "kilowatt hour",
                    "electricity, low voltage",
                )
            ] = maximum

        for s in self.scope["size"]:
            for pt in self.scope["powertrain"]:
                for y in self.scope["year"]:
                    maximum += 1
                    name = "Passenger car, " + pt + ", " + s + ", " + str(y)
                    self.inputs[
                        (name, "GLO", "kilometer", "transport, passenger car, EURO6")
                    ] = maximum

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

        new_A = np.identity(
            np.shape(initial_A)[0] + (len(self.scope["year"]) * 2) + self.number_of_cars
        )

        new_A[0 : np.shape(initial_A)[0], 0 : np.shape(initial_A)[0]] = initial_A
        # Resize the matrix to fit the number of iterations in `array`
        new_A = np.resize(new_A, (self.array.shape[1], new_A.shape[0], new_A.shape[1]))
        return new_A

    def get_B_matrix(self):
        """
        Load the B matrix. The B matrix contains impact assessment figures for a give impact assessment method,
        per unit of activity. Its length column-wise equals the length of the A matrix row-wise.
        Its length row-wise equals the number of impact assessment methods.

        :param method: only "recipe" available at the moment.
        :param level: only "midpoint" available at the moment.
        :return: an array with impact values per unit of activity for each method.
        :rtype: numpy.ndarray

        """

        if self.scenario == "BAU":
            list_file_names = [
                "B_matrix_recipe_midpoint_BAU_2020.csv",
                "B_matrix_recipe_midpoint_BAU_2030.csv",
                "B_matrix_recipe_midpoint_BAU_2040.csv",
                "B_matrix_recipe_midpoint_BAU_2050.csv",
            ]

        if self.scenario == "RCP26":
            list_file_names = [
                "B_matrix_recipe_midpoint_RCP26_2020.csv",
                "B_matrix_recipe_midpoint_RCP26_2030.csv",
                "B_matrix_recipe_midpoint_RCP26_2040.csv",
                "B_matrix_recipe_midpoint_RCP26_2050.csv",
            ]

        B = np.zeros((len(list_file_names), 19, self.A.shape[1]))

        for f in list_file_names:
            filepath = DATA_DIR / "IAM" / f
            initial_B = np.genfromtxt(filepath, delimiter=";")

            new_B = np.zeros(
                (
                    np.shape(initial_B)[0],
                    np.shape(initial_B)[1]
                    + (self.number_of_cars + (len(self.scope["year"]) * 2)),
                )
            )

            new_B[0 : np.shape(initial_B)[0], 0 : np.shape(initial_B)[1]] = initial_B

            B[list_file_names.index(f), :, :] = new_B

        response = xr.DataArray(
            B,
            coords=[
                [2020, 2030, 2040, 2050],
                self.get_dict_impact_categories()["recipe"]["midpoint"],
                list(self.inputs.keys()),
            ],
            dims=["year", "category", "activity",],
        )
        return response

    def get_dict_input(self):
        """
        Load a dictionary with tuple ("name of activity", "location", "unit", "reference product") as key, row/column
        indices as values.

        :return: dictionary with `label:index` pairs.
        :rtype: dict

        """
        filename = "dict_inputs_A_matrix.csv"
        filepath = DATA_DIR / filename
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
        filepath = DATA_DIR / filename
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

    def get_index_from_array(
        self, items_to_look_for, items_to_look_for_also=None, method="or"
    ):
        """
        Return list of row/column indices of self.array of labels that contain the string defined in `items_to_look_for`.

        :param items_to_look_for: string to search for
        :return: list
        """
        # return [
        #    self.inputs[c] - (len(self.inputs) + self.number_of_cars)
        #    for c in self.inputs
        #    if any(ele in c[0] for ele in items_to_look_for)
        #    and (self.inputs[c] - (len(self.inputs) + self.number_of_cars)) >= 0
        # ]
        if method == "or":
            return [
                self.inputs[c] - (len(self.inputs) - self.number_of_cars)
                for c in self.inputs
                if any(ele in c[0] for ele in items_to_look_for)
                and (self.inputs[c] - (len(self.inputs) - self.number_of_cars)) >= 0
            ]
        if method == "and":
            return [
                self.inputs[c] - (len(self.inputs) - self.number_of_cars)
                for c in self.inputs
                if any(ele in c[0] for ele in items_to_look_for)
                and any(ele in c[0] for ele in items_to_look_for_also)
                and (self.inputs[c] - (len(self.inputs) - self.number_of_cars)) >= 0
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
                if all(ele in c[0].lower() for ele in items_to_look_for)
            ]
        if search_by == "compartment":
            return [
                int(self.inputs[c])
                for c in self.inputs
                if all(ele in c[1] for ele in items_to_look_for)
            ]

    def export_lci(self, presamples=True, ecoinvent_compatibility=True):
        """
        Export the inventory as a dictionary. Also return a list of arrays that contain pre-sampled random values if
        :meth:`stochastic` of :class:`CarModel` class has been called.

        :param presamples: boolean.
        :return: inventory, and optionally, list of arrays containing pre-sampled values.
        :rtype: list
        """
        self.set_inputs_in_A_matrix(self.array.values)
        if presamples == True:
            lci, array = ExportInventory(self.A, self.rev_inputs).write_lci(
                presamples, ecoinvent_compatibility
            )
            return (lci, array)
        else:
            lci = ExportInventory(self.A, self.rev_inputs).write_lci(
                presamples, ecoinvent_compatibility
            )
            return lci

    def export_lci_to_bw(self, presamples=True, ecoinvent_compatibility=True):
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
            lci, array = ExportInventory(self.A, self.rev_inputs).write_lci_to_bw(
                presamples, ecoinvent_compatibility
            )
            return (lci, array)
        else:
            lci = ExportInventory(self.A, self.rev_inputs).write_lci_to_bw(
                presamples, ecoinvent_compatibility
            )
            return lci

    def export_lci_to_excel(self, directory=None, ecoinvent_compatibility=True):
        """
        Export the inventory as an Excel file. Also return the file path where the file is stored.

        :param directory: directory where to save the file.
        :type directory: str
        :return: file path where the file is stored.
        :rtype: str
        """

        self.set_inputs_in_A_matrix(self.array.values)
        fp = ExportInventory(self.A, self.rev_inputs).write_lci_to_excel(
            directory, ecoinvent_compatibility
        )
        return fp

    def set_inputs_in_A_matrix(self, array):
        """
        Fill-in the A matrix. Does not return anything. Modifies in place.
        Shape of the A matrix (values, products, activities).

        :param array: :attr:`array` from :class:`CarModel` class
        """

        # Glider
        self.A[
            :,
            self.inputs[
                (
                    "market for glider, passenger car",
                    "GLO",
                    "kilogram",
                    "glider, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            (array[self.array_inputs["glider base mass"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                ("Glider lightweighting", "GLO", "kilogram", "Glider lightweighting")
            ],
            -self.number_of_cars :,
        ] = (
            (array[self.array_inputs["lightweighting"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "maintenance, passenger car",
                    "RER",
                    "unit",
                    "passenger car maintenance",
                )
            ],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["curb mass"], :] / 1600 / 150000 * -1)

        # Glider EoL
        self.A[
            :,
            self.inputs[
                (
                    "market for manual dismantling of used electric passenger car",
                    "GLO",
                    "unit",
                    "manual dismantling of used electric passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["curb mass"], :]
            * (1 - array[self.array_inputs["combustion power share"], :])
            / 1180
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for manual dismantling of used passenger car with internal combustion engine",
                    "GLO",
                    "unit",
                    "manual dismantling of used passenger car with internal combustion engine",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["curb mass"], :]
            * array[self.array_inputs["combustion power share"], :]
            / 1600
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Powertrain components
        self.A[
            :,
            self.inputs[
                (
                    "market for charger, electric passenger car",
                    "GLO",
                    "kilogram",
                    "charger, electric passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["charger mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for converter, for electric passenger car",
                    "GLO",
                    "kilogram",
                    "converter, for electric passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["converter mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for electric motor, electric passenger car",
                    "GLO",
                    "kilogram",
                    "electric motor, electric passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["electric engine mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for inverter, for electric passenger car",
                    "GLO",
                    "kilogram",
                    "inverter, for electric passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["inverter mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for power distribution unit, for electric passenger car",
                    "GLO",
                    "kilogram",
                    "power distribution unit, for electric passenger car",
                )
            ],
            -self.number_of_cars :,
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

        self.A[
            :,
            self.inputs[
                (
                    "market for used powertrain from electric passenger car, manual dismantling",
                    "GLO",
                    "kilogram",
                    "used powertrain from electric passenger car, manual dismantling",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)
            / array[self.array_inputs["lifetime kilometers"], :]
        )

        self.A[
            :,
            self.inputs[
                (
                    "market for internal combustion engine, passenger car",
                    "GLO",
                    "kilogram",
                    "internal combustion engine, for passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            (
                array[
                    [
                        self.array_inputs[l]
                        for l in ["combustion engine mass", "powertrain mass"]
                    ],
                    :,
                ].sum(axis=0)
            )
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[("Ancillary BoP", "GLO", "kilogram", "Ancillary BoP")],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["fuel cell ancillary BoP mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[("Essential BoP", "GLO", "kilogram", "Essential BoP")],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["fuel cell essential BoP mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[("Stack", "GLO", "kilowatt", "Stack")],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["fuel cell stack mass"], :]
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Start of printout

        print(
            "****************** IMPORTANT BACKGROUND PARAMETERS ******************",
            end="\n * ",
        )

        # Energy chain
        dict_map = {
            "Hydro": (
                "electricity production, hydro, run-of-river",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Nuclear": (
                "electricity production, nuclear, pressure water reactor",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Gas": (
                "electricity production, natural gas, conventional power plant",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Solar": (
                "electricity production, photovoltaic, 3kWp slanted-roof installation, multi-Si, panel, mounted",
                "DE",
                "kilowatt hour",
                "electricity, low voltage",
            ),
            "Wind": (
                "electricity production, wind, 1-3MW turbine, onshore",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Biomass": (
                "heat and power co-generation, wood chips, 6667 kW, state-of-the-art 2014",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Coal": (
                "electricity production, hard coal",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Oil": (
                "electricity production, oil",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Geo": (
                "electricity production, deep geothermal",
                "DE",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            "Waste": (
                "treatment of municipal solid waste, incineration",
                "DE",
                "kilowatt hour",
                "electricity, for reuse in municipal waste incineration only",
            ),
        }

        # Energy storage

        print(
            "The country of use is " + self.background_configuration["country"],
            end="\n * ",
        )

        battery_tech = self.background_configuration["battery technology"]
        battery_origin = self.background_configuration["battery origin"]

        print(
            "Power and energy batteries produced in "
            + battery_origin
            + " using "
            + battery_tech
            + " chemistry.",
            end="\n * ",
        )

        losses_to_low = float(self.bs.losses[battery_origin]["LV"])
        mix_battery_manufacturing = (
            self.bs.electricity_mix.sel(country=battery_origin, value=0)
            .interp(year=self.scope["year"])
            .values
        )

        # Fill the electricity markets for battery production
        for y in self.scope["year"]:
            m = np.array(
                mix_battery_manufacturing[self.scope["year"].index(y)]
            ).reshape(-1, 10, 1)

            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    [self.inputs[dict_map[t]] for t in dict_map],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(y) in i[0]
                        and "electricity market for energy storage production" in i[0]
                    ],
                )
            ] = (m * losses_to_low * -1)

            # Add transmission network for high and medium voltage
            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, electricity, high voltage",
                        "CH",
                        "kilometer",
                        "transmission network, electricity, high voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (6.58e-9 * -1 * losses_to_low)

            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, electricity, medium voltage",
                        "CH",
                        "kilometer",
                        "transmission network, electricity, medium voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (1.86e-8 * -1 * losses_to_low)

            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, long-distance",
                        "UCTE",
                        "kilometer",
                        "transmission network, long-distance",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (3.17e-10 * -1 * losses_to_low)

            # Add distribution network, low voltage
            self.A[
                :,
                self.inputs[
                    (
                        "distribution network construction, electricity, low voltage",
                        "CH",
                        "kilometer",
                        "distribution network, electricity, low voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (8.74e-8 * -1 * losses_to_low)

            # Add supply of sulfur hexafluoride for transformers
            self.A[
                :,
                self.inputs[
                    (
                        "market for sulfur hexafluoride, liquid",
                        "RER",
                        "kilogram",
                        "sulfur hexafluoride, liquid",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = ((5.4e-8 + 2.99e-9) * -1 * losses_to_low)

            # Add SF_6 leakage

            self.A[
                :,
                self.inputs[("Sulfur hexafluoride", ("air",), "kilogram")],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = ((5.4e-8 + 2.99e-9) * -1 * losses_to_low)

        # Use the NMC inventory of Schmidt et al. 2019
        self.A[
            :,
            self.inputs[("Battery BoP", "GLO", "kilogram", "Battery BoP")],
            -self.number_of_cars :,
        ] = (
            (
                array[self.array_inputs["battery BoP mass"], :]
                * (1 + array[self.array_inputs["battery lifetime replacements"], :])
            )
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        battery_cell_label = (
            "Battery cell, " + battery_tech,
            "GLO",
            "kilogram",
            "Battery cell",
        )

        self.A[:, self.inputs[battery_cell_label], -self.number_of_cars :,] = (
            (
                array[self.array_inputs["battery cell mass"], :]
                * (1 + array[self.array_inputs["fuel cell lifetime replacements"], :])
            )
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Set an input of electricity, given the country of manufacture
        self.A[
            :,
            self.inputs[
                (
                    "market group for electricity, medium voltage",
                    "World",
                    "kilowatt hour",
                    "electricity, medium voltage",
                )
            ],
            self.inputs[battery_cell_label],
        ] = 0

        for y in self.scope["year"]:

            index = self.get_index_from_array([str(y)])

            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(y) in i[0]
                        and "electricity market for energy storage production" in i[0]
                    ],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(y) in i[0] and "Passenger" in i[0]
                    ],
                )
            ] = (
                array[
                    self.array_inputs["battery cell production electricity"], :, index
                ].T
                * self.A[
                    :,
                    self.inputs[battery_cell_label],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(y) in i[0] and "Passenger" in i[0]
                    ],
                ]
            ).reshape(
                self.iterations, 1, -1
            )

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if any(
                ele in c[0]
                for ele in ["ICEV-d", "ICEV-p", "HEV-p", "PHEV-p", "PHEV-d", "HEV-d"]
            )
        ]
        index = self.get_index_from_array(
            ["ICEV-d", "ICEV-p", "HEV-p", "PHEV-p", "PHEV-d", "HEV-d"]
        )

        self.A[
            :,
            self.inputs[
                (
                    "polyethylene production, high density, granulate",
                    "RER",
                    "kilogram",
                    "polyethylene, high density, granulate",
                )
            ],
            index_A,
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        index = self.get_index_from_array(["ICEV-g"])
        self.A[
            :,
            self.inputs[
                (
                    "glass fibre reinforced plastic production, polyamide, injection moulded",
                    "RER",
                    "kilogram",
                    "glass fibre reinforced plastic, polyamide, injection moulded",
                )
            ],
            self.index_cng,
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        index = self.get_index_from_array(["FCEV"])
        self.A[
            :,
            self.inputs[
                (
                    "Fuel tank, compressed hydrogen gas, 700bar",
                    "GLO",
                    "kilogram",
                    "Fuel tank, compressed hydrogen gas, 700bar",
                )
            ],
            self.index_fuel_cell,
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        if not any(
            True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in self.scope["powertrain"]
        ):
            self.background_configuration["custom electricity mix"] = [
                self.bs.electricity_mix.sel(
                    country=self.background_configuration["country"], value=0
                )
                .interp(year=y)
                .values
                for y in self.scope["year"]
            ]

        country = self.background_configuration["country"]
        try:
            losses_to_low = float(self.bs.losses[country]["LV"])
        except KeyError:
            # If losses for the country are not found, assume 15%
            losses_to_low = 1.15

        if "custom electricity mix" in self.background_configuration or not any(
            True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in self.scope["powertrain"]
        ):
            # If a special electricity mix is specified, we use it
            mix = self.background_configuration["custom electricity mix"]
            print("The electricity mix is provided by the user.", end="\n \t * ")
        else:
            print("The electricity mix is provided by carculator.", end="\n \t * ")
            use_year = [
                int(i)
                for i in (
                    array[
                        self.array_inputs["lifetime kilometers"],
                        :,
                        self.get_index_from_array(["BEV", "PHEV-p", "PHEV-d"]),
                    ]
                    / array[
                        self.array_inputs["kilometers per year"],
                        :,
                        self.get_index_from_array(["BEV", "PHEV-p", "PHEV-d"]),
                    ]
                )
                .mean(axis=1)
                .reshape(-1, len(self.scope["year"]))
                .mean(axis=0)
            ]

            mix = [
                self.bs.electricity_mix.sel(country=country, value=0)
                .interp(year=np.arange(y, y + use_year[self.scope["year"].index(y)]))
                .mean(axis=0)
                .values
                for y in self.scope["year"]
            ]

        for y in self.scope["year"]:
            sum_renew = (
                mix[self.scope["year"].index(y)][0]
                + mix[self.scope["year"].index(y)][4]
                + mix[self.scope["year"].index(y)][5]
                + mix[self.scope["year"].index(y)][8]
            )

            if self.scope["year"].index(y) + 1 == len(self.scope["year"]):
                end_str = "\n * "
            else:
                end_str = "\n \t * "

            print(
                "in "
                + str(y)
                + ", % of renewable _________________________ "
                + str(np.round(sum_renew * 100, 0))
                + "%",
                end=end_str,
            )

        # Fill the electricity markets for battery charging and hydrogen production
        for y in self.scope["year"]:
            m = np.array(mix[self.scope["year"].index(y)]).reshape(-1, 10, 1)
            # Add electricity technology shares
            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    [self.inputs[dict_map[t]] for t in dict_map],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(y) in i[0]
                        and "electricity market for fuel preparation" in i[0]
                    ],
                )
            ] = (m * -1 * losses_to_low)

            # Add transmission network for high and medium voltage
            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, electricity, high voltage",
                        "CH",
                        "kilometer",
                        "transmission network, electricity, high voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (6.58e-9 * -1 * losses_to_low)

            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, electricity, medium voltage",
                        "CH",
                        "kilometer",
                        "transmission network, electricity, medium voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (1.86e-8 * -1 * losses_to_low)

            self.A[
                :,
                self.inputs[
                    (
                        "transmission network construction, long-distance",
                        "UCTE",
                        "kilometer",
                        "transmission network, long-distance",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (3.17e-10 * -1 * losses_to_low)

            # Add distribution network, low voltage
            self.A[
                :,
                self.inputs[
                    (
                        "distribution network construction, electricity, low voltage",
                        "CH",
                        "kilometer",
                        "distribution network, electricity, low voltage",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (8.74e-8 * -1 * losses_to_low)

            # Add supply of sulfur hexafluoride for transformers
            self.A[
                :,
                self.inputs[
                    (
                        "market for sulfur hexafluoride, liquid",
                        "RER",
                        "kilogram",
                        "sulfur hexafluoride, liquid",
                    )
                ],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = ((5.4e-8 + 2.99e-9) * -1 * losses_to_low)

            # Add SF_6 leakage

            self.A[
                :,
                self.inputs[("Sulfur hexafluoride", ("air",), "kilogram")],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = ((5.4e-8 + 2.99e-9) * -1 * losses_to_low)

        if any(
            True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in self.scope["powertrain"]
        ):
            for y in self.scope["year"]:
                index = self.get_index_from_array(
                    ["BEV", "PHEV-p", "PHEV-d"], [str(y)], method="and"
                )

                self.A[
                    np.ix_(
                        np.arange(self.iterations),
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "electricity market for fuel preparation" in i[0]
                        ],
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "Passenger" in i[0]
                            and any(
                                True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in i[0]
                            )
                        ],
                    )
                ] = (
                    array[self.array_inputs["electricity consumption"], :, index] * -1
                ).T.reshape(
                    self.iterations, 1, -1
                )

        if "FCEV" in self.scope["powertrain"]:

            index = self.get_index_from_array(["FCEV"])

            if "hydrogen technology" in self.background_configuration:
                # If a customization dict is passed
                hydro_technology = self.background_configuration["hydrogen technology"]
            else:
                hydro_technology = "Electrolysis"

            dict_h_map = {
                "Electrolysis": (
                    "Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station",
                    "RER",
                    "kilogram",
                    "Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station",
                ),
                "SMR": (
                    "Hydrogen, gaseous, 700 bar, from SMR NG w/o CCS, at H2 fuelling station",
                    "RER",
                    "kilogram",
                    "Hydrogen, gaseous, 700 bar, from SMR NG w/o CCS, at H2 fuelling station",
                ),
            }

            self.A[
                :, self.inputs[dict_h_map[hydro_technology]], self.index_fuel_cell
            ] = (
                array[self.array_inputs["fuel mass"], :, index]
                / array[self.array_inputs["range"], :, index]
                * -1
            ).T

            print("Hydrogen produced by " + hydro_technology + ".", end="\n * ")

            # If hydrolysis is chosen, adjust the electricity mix

            if hydro_technology == "Electrolysis":
                index_FCEV = self.get_index_from_array(["FCEV"])

                for y in self.scope["year"]:
                    index = [
                        x
                        for x in self.get_index_from_array([str(y)])
                        if x in index_FCEV
                    ]
                    new_mix = np.zeros((len(index), 10, self.iterations))
                    m = np.array(mix[self.scope["year"].index(y)]).reshape(-1, 10)
                    new_mix[:, np.arange(10), :] = m.T
                    b = (
                        array[self.array_inputs["fuel mass"], :, index]
                        / array[self.array_inputs["range"], :, index]
                    )

                    self.A[
                        :,
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "electricity market for fuel preparation" in i[0]
                        ],
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0] and "Passenger" in i[0] and "FCEV" in i[0]
                        ],
                    ] = (
                        (
                            array[self.array_inputs["fuel mass"], :, index]
                            / array[self.array_inputs["range"], :, index]
                        )
                        * -58
                    ).T

        if "ICEV-g" in self.scope["powertrain"]:
            index = self.get_index_from_array(["ICEV-g"])

            if "cng technology" in self.background_configuration:
                # If a customization dict is passed
                cng_technology = self.background_configuration["cng technology"]
            else:
                cng_technology = "biogas"

            if "alternative CNG share" in self.background_configuration:
                # If a share of alternative fuel to CNG is provided
                share = np.array(self.background_configuration["alternative CNG share"])
                try:
                    assert len(share) == len(self.scope["year"])

                except AssertionError:
                    print(
                        "The values for biofuel shares do not match the number of years."
                    )
                    exit(1)

            else:
                if "biogas" in cng_technology:
                    fuel_type = "Biomass fuel"
                else:
                    fuel_type = "Synthetic fuel"

                region = self.bs.region_map[self.background_configuration["country"]][
                    "RegionCode"
                ]
                share = (
                    self.bs.biofuel.sel(
                        region=region,
                        value=0,
                        fuel_type=fuel_type,
                        scenario=self.scenario,
                    )
                    .interp(year=self.scope["year"])
                    .values
                )
                self.background_configuration["alternative CNG share"] = share

            print("CNG is completed by " + cng_technology + ".", end="\n \t * ")

            for y in self.scope["year"]:
                if self.scope["year"].index(y) + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "
                print(
                    "in "
                    + str(y)
                    + " _________________________________________ "
                    + str(np.round(share[self.scope["year"].index(y)] * 100, 0))
                    + "%",
                    end=end_str,
                )

            dict_cng_map = {
                "biogas": (
                    "biogas upgrading - sewage sludge - amine scrubbing - best",
                    "CH",
                    "kilogram",
                    "biogas upgrading - sewage sludge - amine scrubbing - best",
                ),
                "syngas": (
                    "methane, from electrochemical methanation, with carbon from atmospheric CO2 capture",
                    "RER",
                    "kilogram",
                    "methane, from electrochemical methanation",
                ),
            }

            # Conventional CNG minus biogas or syngas share
            for y in self.scope["year"]:
                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0] and "Passenger" in i[0] and "ICEV-g" in i[0]
                ]
                ind_array = [
                    x for x in self.get_index_from_array([str(y)]) if x in index
                ]

                self.A[
                    :,
                    self.inputs[
                        (
                            "market for natural gas, from high pressure network (1-5 bar), at service station",
                            "GLO",
                            "kilogram",
                            "natural gas, from high pressure network (1-5 bar), at service station",
                        )
                    ],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        * (1 - share[self.scope["year"].index(y)])
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # biogas share
                self.A[:, self.inputs[dict_cng_map[cng_technology]], ind_A] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        * share[self.scope["year"].index(y)]
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # Fuel-based emissions from CNG, CO2
                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["CO2 per kg fuel"], :, ind_array]
                        * (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * (1 - share[self.scope["year"].index(y)])
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if [i for i in self.scope["powertrain"] if i in ["ICEV-d", "PHEV-d", "HEV-d"]]:
            index = self.get_index_from_array(["ICEV-d", "PHEV-d", "HEV-d"])

            if "diesel technology" in self.background_configuration:
                # If a customization dict is passed
                diesel_technology = self.background_configuration["diesel technology"]
            else:
                diesel_technology = "biodiesel - algae"

            if "alternative diesel share" in self.background_configuration:
                # If a share of alternative fuel to diesel is provided
                share = np.array(
                    self.background_configuration["alternative diesel share"]
                )
                try:
                    assert len(share) == len(self.scope["year"])

                except AssertionError:
                    print(
                        "The values for biofuel shares do not match the number of years."
                    )
                    exit(1)

            else:

                if "biodiesel" in diesel_technology:
                    fuel_type = "Biomass fuel"
                else:
                    fuel_type = "Synthetic fuel"

                region = self.bs.region_map[self.background_configuration["country"]][
                    "RegionCode"
                ]
                share = (
                    self.bs.biofuel.sel(
                        region=region,
                        value=0,
                        fuel_type=fuel_type,
                        scenario=self.scenario,
                    )
                    .interp(year=self.scope["year"])
                    .values
                )
                self.background_configuration["alternative diesel share"] = share

            dict_diesel_map = {
                "diesel": (
                    "market for diesel",
                    "Europe without Switzerland",
                    "kilogram",
                    "diesel",
                ),
                "biodiesel - algae": (
                    "Biodiesel from algae",
                    "RER",
                    "kilogram",
                    "Biodiesel from algae",
                ),
                "biodiesel - cooking oil": (
                    "Biodiesel from cooking oil",
                    "RER",
                    "kilogram",
                    "Biodiesel from cooking oil",
                ),
                "synthetic diesel": (
                    "Diesel production, Fischer Tropsch process",
                    "RER",
                    "kilogram",
                    "Diesel",
                ),
            }

            print("Diesel is completed by " + diesel_technology + ".", end="\n \t * ")

            for y in self.scope["year"]:
                if self.scope["year"].index(y) + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "
                print(
                    "in "
                    + str(y)
                    + " _________________________________________ "
                    + str(np.round(share[self.scope["year"].index(y)] * 100, 0))
                    + "%",
                    end=end_str,
                )

            for y in self.scope["year"]:
                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "Passenger" in i[0]
                    and any(x in i[0] for x in ["ICEV-d", "PHEV-d", "HEV-d"])
                ]

                ind_array = [
                    x for x in self.get_index_from_array([str(y)]) if x in index
                ]

                # conventional diesel minus biodiesel share
                self.A[
                    :,
                    self.inputs[
                        (
                            "market for diesel",
                            "Europe without Switzerland",
                            "kilogram",
                            "diesel",
                        )
                    ],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        * (1 - share[self.scope["year"].index(y)])
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # biodiesel
                self.A[:, self.inputs[dict_diesel_map[diesel_technology]], ind_A,] = (
                    (
                        (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * share[self.scope["year"].index(y)]
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # Fuel-based emissions from conventional diesel, CO2
                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["CO2 per kg fuel"], :, ind_array]
                        * (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * (1 - share[self.scope["year"].index(y)])
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if [i for i in self.scope["powertrain"] if i in ["ICEV-p", "HEV-p", "PHEV-p"]]:
            index = self.get_index_from_array(["ICEV-p", "HEV-p", "PHEV-p"])

            if "petrol technology" in self.background_configuration:
                # If a customization dict is passed
                petrol_technology = self.background_configuration["petrol technology"]
            else:
                petrol_technology = "bioethanol - wheat straw"

            if "alternative petrol share" in self.background_configuration:
                # If a share of alternative fuel to gasoline is provided
                share = np.array(
                    self.background_configuration["alternative petrol share"]
                )
                try:
                    assert len(share) == len(self.scope["year"])

                except AssertionError:
                    print(
                        "The values for biofuel shares do not match the number of years."
                    )
                    exit(1)

            else:
                if "bioethanol" in petrol_technology:
                    fuel_type = "Biomass fuel"
                else:
                    fuel_type = "Synthetic fuel"

                region = self.bs.region_map[self.background_configuration["country"]][
                    "RegionCode"
                ]

                share = (
                    self.bs.biofuel.sel(
                        region=region,
                        value=0,
                        fuel_type=fuel_type,
                        scenario=self.scenario,
                    )
                    .interp(year=self.scope["year"])
                    .values
                )
                self.background_configuration["alternative petrol share"] = share

            dict_petrol_map = {
                "bioethanol - wheat straw": (
                    "Ethanol from wheat straw pellets",
                    "RER",
                    "kilogram",
                    "Ethanol from wheat straw pellets",
                ),
                "bioethanol - forest residues": (
                    "Ethanol from forest residues",
                    "RER",
                    "kilogram",
                    "Ethanol from forest residues",
                ),
                "bioethanol - sugarbeet": (
                    "Ethanol from sugarbeet",
                    "RER",
                    "kilogram",
                    "Ethanol from sugarbeet",
                ),
                "bioethanol - maize starch": (
                    "Ethanol from maize starch",
                    "RER",
                    "kilogram",
                    "Ethanol from maize starch",
                ),
                "synthetic gasoline": (
                    "Gasoline production",
                    "RER",
                    "kilogram",
                    "Gasoline",
                ),
            }

            print("Gasoline is completed by " + petrol_technology + ".", end="\n \t * ")

            for y in self.scope["year"]:
                if self.scope["year"].index(y) + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "
                print(
                    "in "
                    + str(y)
                    + " _________________________________________ "
                    + str(np.round(share[self.scope["year"].index(y)] * 100, 0))
                    + "%",
                    end=end_str,
                )

            for y in self.scope["year"]:
                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(y) in i[0]
                    and "Passenger" in i[0]
                    and any(x in i[0] for x in ["ICEV-p", "HEV-p", "PHEV-p"])
                ]
                ind_array = [
                    x for x in self.get_index_from_array([str(y)]) if x in index
                ]
                # conventional petrol supply minus share of bioethanol
                self.A[
                    :,
                    self.inputs[
                        (
                            "market for petrol, low-sulfur",
                            "Europe without Switzerland",
                            "kilogram",
                            "petrol, low-sulfur",
                        )
                    ],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        * (1 - share[self.scope["year"].index(y)])
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # bioethanol supply
                self.A[:, self.inputs[dict_petrol_map[petrol_technology]], ind_A,] = (
                    (
                        (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * share[self.scope["year"].index(y)]
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # Fuel-based emissions from conventional petrol, CO2

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["CO2 per kg fuel"], :, ind_array]
                        * (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * (1 - share[self.scope["year"].index(y)])
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        # Non-exhaust emissions
        self.A[
            :,
            self.inputs[
                (
                    "market for road wear emissions, passenger car",
                    "GLO",
                    "kilogram",
                    "road wear emissions, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 1e-08)
        self.A[
            :,
            self.inputs[
                (
                    "market for tyre wear emissions, passenger car",
                    "GLO",
                    "kilogram",
                    "tyre wear emissions, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 6e-08)
        self.A[
            :,
            self.inputs[
                (
                    "market for brake wear emissions, passenger car",
                    "GLO",
                    "kilogram",
                    "brake wear emissions, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (array[self.array_inputs["driving mass"], :] * 5e-09)

        # Infrastructure
        self.A[
            :,
            self.inputs[("market for road", "GLO", "meter-year", "road")],
            -self.number_of_cars :,
        ] = (5.37e-7 * array[self.array_inputs["driving mass"], :] * -1)

        # Exhaust emissions

        # Non-fuel based emissions
        self.A[:, self.index_emissions, -self.number_of_cars :] = (
            array[
                [
                    self.array_inputs[self.map_non_fuel_emissions[self.rev_inputs[x]]]
                    for x in self.index_emissions
                ]
            ]
            * -1
        ).transpose([1, 0, 2])

        # Noise emissions
        self.A[:, self.index_noise, -self.number_of_cars :] = (
            array[
                [self.array_inputs[self.rev_inputs[e][0]] for e in self.index_noise], :
            ]
            * -1
        ).transpose([1, 0, 2])
        print("*********************************************************************")
