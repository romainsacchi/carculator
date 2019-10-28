"""
.. module: inventory.py

"""
import numpy as np
from pathlib import Path
from inspect import currentframe, getframeinfo
import csv
import xarray as xr
import itertools
from . import DATA_DIR


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
        ]
        self.index_noise = [self.inputs[i] for i in self.inputs if "noise" in i[0]]

        self.split_dict = self.get_split_dict()

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
                "energy storage",
                "glider",
                "powertrain",
                "road",
                "maintenance",
                "other"
            ]

        dict_impact_cat = self.get_dict_impact_categories()

        if FU == None:
            response = xr.DataArray(
                np.zeros(
                    (
                        self.B.shape[0],
                        len(self.size),
                        len(self.powertrain),
                        len(self.year),
                        len(cat),
                        self.iterations,
                    )
                ),
                coords=[
                    dict_impact_cat[method][level],
                    self.size,
                    self.powertrain,
                    self.year,
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

    def get_split_dict(self):
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

        return csv_dict

    def calculate_impacts(
        self, FU=None, method="recipe", level="midpoint", split="components"
    ):

        # Load the B matrix
        self.B = self.get_B_matrix(method, level)

        # Prepare an array to store the results
        self.results = self.get_results_table(FU, method, level, split)

        # List the splitting categories
        split_categories = [cat for cat in self.results.coords['impact'].values
                      if cat != 'other']

        # Iterate through the number of iterations
        for i in range(self.iterations):
            self.temp_array = self.array.sel(value=i).values
            self.set_inputs_in_A_matrix()

            # TODO: optimize this whole section
            for pt in self.powertrain:
                for y in self.year:
                    for s in self.size:
                        # Retrieve the index of a given car in the matrix A
                        car = self.inputs[('Passenger car, '+ pt + ', ' + s + ', ' + str(y),"GLO")]
                        # Set the demand vector with zeros and a 1 corresponding to the car position in the vector
                        f = np.zeros(np.shape(self.A)[0])
                        f[car] = 1

                        # Solve inventory
                        g = np.linalg.inv(self.A).dot(f)
                        C = g * self.B


                        # Iterate through the results array to fill it
                        for cat in split_categories:
                            # Retrieve position of certain datasets to split results into categories
                            # (direct emissions, energy chain, etc.)
                            index = [self.get_index_of_flows([l['search for']], l['search by'])
                                     for l in self.split_dict[split][cat]]
                            index = [item for sublist in index for item in sublist]

                            self.results.loc[dict(impact=cat, year=y, size=s, powertrain=pt, value=i)] = \
                                C[:, index].sum(axis=1)
                        # Fill the 'other' section by subtracting the total impact by what has already been
                        # accounted for.
                        self.results.loc[dict(impact='other', year=y, size=s, powertrain=pt, value=i)] = \
                            (C[:, :].sum(axis=1) - self.results.loc[dict(year=y, size=s, powertrain=pt, value=i)].sum(axis=1).values)

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
                    csv_dict[(row[0], t)] = count
                else:
                    csv_dict[(row[0], row[1])] = count
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
                self.inputs[c]
                for c in self.inputs
                if any(ele in c[0] for ele in items_to_look_for)
            ]
        if search_by == "compartment":
            return [
                self.inputs[c]
                for c in self.inputs
                if any(ele in c[1] for ele in items_to_look_for)
            ]

    def set_inputs_in_A_matrix(self):

        # Glider
        self.A[
            self.inputs[("market for glider, passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["glider base mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[("Glider lightweighting", "GLO")], -self.number_of_cars :
        ] = (
            self.temp_array[self.array_inputs["lightweighting"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[("maintenance, passenger car", "RER")], -self.number_of_cars :
        ] = (self.temp_array[self.array_inputs["curb mass"], :] / 1600 / 150000 * -1)

        # Glider EoL
        self.A[
            self.inputs[
                ("market for manual dismantling of used electric passenger car", "GLO")
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
                (
                    "market for manual dismantling of used passenger car with internal combustion engine",
                    "GLO",
                )
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
            self.inputs[("market for charger, electric passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["charger mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[("market for converter, for electric passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["converter mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[("market for electric motor, electric passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["electric engine mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[("market for inverter, for electric passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (
            self.temp_array[self.array_inputs["inverter mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            self.inputs[
                (
                    "market for power distribution unit, for electric passenger car",
                    "GLO",
                )
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
                (
                    "market for used powertrain from electric passenger car, manual dismantling",
                    "GLO",
                )
            ],
            -self.number_of_cars :,
        ] = (
            self.temp_array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
        )

        self.A[
            self.inputs[
                ("market for internal combustion engine, passenger car", "GLO")
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

        self.A[self.inputs[("Ancillary BoP", "GLO")], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell ancillary BoP mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[("Essential BoP", "GLO")], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell essential BoP mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[self.inputs[("Stack 2020", "GLO")], -self.number_of_cars :] = (
            self.temp_array[self.array_inputs["fuel cell stack mass"], :]
            / self.temp_array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Energy storage
        self.A[self.inputs[("Battery BoP", "GLO")], -self.number_of_cars :] = (
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

        self.A[self.inputs[("Battery cell", "GLO")], -self.number_of_cars :] = (
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

        self.A[self.inputs[("Battery cell", "GLO")], -self.number_of_cars :] = (
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
            self.inputs[("polyethylene production, high density, granulate", "RER")],
            index_A,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        index = self.get_index_from_array(["ICEV-g"])
        self.A[
            self.inputs[
                (
                    "glass fibre reinforced plastic production, polyamide, injection moulded",
                    "RER",
                )
            ],
            self.index_cng,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        index = self.get_index_from_array(["FCEV"])
        self.A[
            self.inputs[("Fuel tank, compressed hydrogen gas, 700bar", "GLO")],
            self.index_fuel_cell,
        ] = (
            self.temp_array[self.array_inputs["fuel tank mass"], index]
            / self.temp_array[self.array_inputs["lifetime kilometers"], index]
            * -1
        )

        # Energy chain
        self.A[
            self.inputs[("market group for electricity, low voltage", "ENTSO-E")],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["electricity consumption"], :] * -1)

        index = self.get_index_from_array(["FCEV"])
        self.A[
            self.inputs[
                (
                    "Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station",
                    "RER",
                )
            ],
            self.index_fuel_cell,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        index = self.get_index_from_array(["ICEV-g"])

        self.A[
            self.inputs[
                (
                    "market for natural gas, from high pressure network (1-5 bar), at service station",
                    "GLO",
                )
            ],
            self.index_cng,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        index = self.get_index_from_array(["ICEV-d"])
        self.A[
            self.inputs[("market for diesel", "Europe without Switzerland")],
            self.index_diesel,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        index = self.get_index_from_array(["ICEV-p"])
        self.A[
            self.inputs[
                ("market for petrol, low-sulfur", "Europe without Switzerland")
            ],
            self.index_petrol,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        index = self.get_index_from_array(["HEV-p"])
        self.A[
            self.inputs[
                ("market for petrol, low-sulfur", "Europe without Switzerland")
            ],
            self.index_hybrid,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        index = self.get_index_from_array(["PHEV"])
        self.A[
            self.inputs[
                ("market for petrol, low-sulfur", "Europe without Switzerland")
            ],
            self.index_plugin_hybrid,
        ] = (
            self.temp_array[self.array_inputs["fuel mass"], index]
            / self.temp_array[self.array_inputs["range"], index]
            * -1
        )

        # Non-exhaust emissions
        self.A[
            self.inputs[("market for road wear emissions, passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 1e-08)
        self.A[
            self.inputs[("market for tyre wear emissions, passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 6e-08)
        self.A[
            self.inputs[("market for brake wear emissions, passenger car", "GLO")],
            -self.number_of_cars :,
        ] = (self.temp_array[self.array_inputs["driving mass"], :] * 5e-09)

        # Infrastructure
        self.A[self.inputs[("market for road", "GLO")], -self.number_of_cars :] = (
            5.37e-7 * self.temp_array[self.array_inputs["driving mass"], :] * -1
        )

        # Exhaust emissions
        # Fuel-based emissions
        self.A[
            self.inputs[("Carbon dioxide, fossil", ("air",))], -self.number_of_cars :
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

        # Noise emissions
        self.A[self.index_noise, -self.number_of_cars :] = (
            self.temp_array[
                [self.array_inputs[self.rev_inputs[e][0]] for e in self.index_noise], :
            ]
            * -1
        )
