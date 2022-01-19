"""
inventory.py contains InventoryCalculation which provides all methods to solve inventories.
"""

import csv
import glob
import itertools
from inspect import currentframe, getframeinfo
from pathlib import Path

import numpy as np
import xarray as xr
import yaml
from scipy import sparse

from . import DATA_DIR
from .background_systems import BackgroundSystemModel
from .export import ExportInventory
from .utils import build_fleet_array

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

IAM_FILES_DIR = DATA_DIR / "IAM"


class InventoryCalculation:
    """
    Build and solve the inventory for results characterization and inventory export

    Vehicles to be analyzed can be filtered by passing a `scope` dictionary.
    Some assumptions in the background system can also be adjusted by passing a `background_configuration` dictionary.

    .. code-block:: python

        scope = {
                        'powertrain':['BEV', 'FCEV', 'ICEV-p'],
                    }
        bc = {'country':'CH', # considers electricity network losses for Switzerland
              'custom electricity mix' : [[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                                          [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                                          [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                                          [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                                         ], # in this case, 100% nuclear for the second year
              'fuel blend':{
                  'cng':{ #specify fuel bland for compressed gas
                        'primary fuel':{
                            'type':'biogas',
                            'share':[0.9, 0.8, 0.7, 0.6], # shares per year. Must total 1 for each year.
                            'origin':'FR' # can specify a location. The geographically closest dataset will be chosen.
                            },
                        'secondary fuel':{
                            'type':'syngas',
                            'share': [0.1, 0.2, 0.3, 0.4],
                            'origin': can specify a location. In case of synthetic fuel, the electricity mix is adjusted
                            # to the location. If not location is specified, it is assumed that synfuels are
                            # locally produced (i.e., same location as country of use).
                            }
                        },
                 'diesel':{
                        'primary fuel':{
                            'type':'synthetic diesel - energy allocation',
                            'share':[0.9, 0.8, 0.7, 0.6]
                            },
                        'secondary fuel':{
                            'type':'biodiesel - cooking oil',
                            'share': [0.1, 0.2, 0.3, 0.4]
                            }
                        },
                 'petrol':{
                        'primary fuel':{
                            'type':'petrol',
                            'share':[0.9, 0.8, 0.7, 0.6]
                            },
                        'secondary fuel':{
                            'type':'bioethanol - wheat straw',
                            'share': [0.1, 0.2, 0.3, 0.4]
                            }
                        },
                'hydrogen':{
                        'primary fuel':{'type':'electrolysis', 'share':[1, 0, 0, 0]},
                        'secondary fuel':{'type':'smr - natural gas', 'share':[0, 1, 1, 1]}
                        }
                    },
              'energy storage': {
                  'electric': {
                      'type':'NMC',
                      'origin': 'NO'
                  },
                  'hydrogen': {
                      'type':'carbon fiber'
                  }
              }
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
        - [10]: biomass with CCS
        - [11]: biogas with CCS
        - [12]: coal with CCS
        - [13]: natural gas with CCS
        - [14]: wood with CCS

    If none is given, the electricity mix corresponding to the country specified in `country` will be selected.
    If no country is specified, Europe applies.

    The `primary` and `secondary` fuel keys contain an array with shares of alternative petrol fuel for each year, to create a custom blend.
    If none is provided, a blend provided by the Integrated Assessment model REMIND is used, which will depend on the REMIND energy scenario selected.

    Here is a list of available fuel pathways:


    Hydrogen technologies
    --------------------
    "electrolysis"
    "smr - natural gas"
    "smr - natural gas with CCS"
    "smr - biogas"
    "smr - biogas with CCS"
    "coal gasification"
    "wood gasification"
    "wood gasification with CCS"
    "wood gasification with EF"
    "wood gasification with EF with CCS"
    "atr - natural gas"
    "atr - natural gas with CCS"
    "atr - biogas"
    "atr - biogas with CCS"


    Natural gas technologies
    ------------------------
    cng
    biogas - sewage sludge
    biogas - biowaste
    syngas

    Diesel technologies
    -------------------
    diesel
    biodiesel - algae
    biodiesel - cooking oil
    synthetic diesel - economic allocation
    synthetic diesel - energy allocation

    Petrol technologies
    -------------------
    petrol
    bioethanol - wheat straw
    bioethanol - maize starch
    bioethanol - sugarbeet
    bioethanol - forest residues
    synthetic gasoline - economic allocation
    synthetic gasoline - energy allocation

    :ivar array: array from the CarModel class
    :vartype array: CarModel.array
    :ivar scope: dictionary that contains filters for narrowing the analysis
    :ivar background_configuration: dictionary that contains choices for background system
    :ivar scenario: REMIND energy scenario to use ("SSP2-Baseline": business-as-usual,
                                                    "SSP2-PkBudg1300": limits temperature increase by 2100 to 2 degrees Celsius,
                                                    "static": no forward-looking modification of the background inventories).
                    "SSP2-Baseline" selected by default.

    .. code-block:: python

    """

    def __init__(
        self,
        array,
        scope=None,
        background_configuration=None,
        scenario="SSP2-Base",
        method="recipe",
        method_type="midpoint",
    ):

        if scope is None:
            scope = {
                "size": array.coords["size"].values.tolist(),
                "powertrain": array.coords["powertrain"].values.tolist(),
                "year": array.coords["year"].values.tolist(),
                "fu": {"unit": "vkm", "quantity": 1},
            }
        else:
            scope["size"] = scope.get("size", array.coords["size"].values.tolist())
            scope["powertrain"] = scope.get(
                "powertrain", array.coords["powertrain"].values.tolist()
            )
            scope["year"] = scope.get("year", array.coords["year"].values.tolist())
            scope["fu"] = scope.get("fu", {"unit": "vkm", "quantity": 1})

            if "unit" not in scope["fu"]:
                scope["fu"]["unit"] = "vkm"
            else:
                if scope["fu"]["unit"] not in ["vkm", "pkm"]:
                    raise NameError(
                        "Incorrect specification of functional unit. Must be 'vkm' or 'pkm'."
                    )

            if "quantity" not in scope["fu"]:
                scope["fu"]["quantity"] = 1
            else:
                try:
                    float(scope["fu"]["quantity"])
                except ValueError:
                    raise ValueError(
                        "Incorrect quantity for the functional unit defined."
                    )

        self.scope = scope
        self.scenario = scenario

        # Check if a fleet composition is specified
        if "fleet" in self.scope["fu"]:

            if isinstance(self.scope["fu"]["fleet"], xr.DataArray):
                self.fleet = self.scope["fu"]["fleet"]
            else:

                # check if a path as string is provided
                if isinstance(self.scope["fu"]["fleet"], str):
                    filepath = Path(self.scope["fu"]["fleet"])

                # check if instance of pathlib is provided instead
                elif isinstance(self.scope["fu"]["fleet"], Path):
                    filepath = self.scope["fu"]["fleet"]

                else:
                    raise TypeError(
                        "The format used to specify fleet compositions is not valid."
                        "A file path that points to a CSV file is expected. "
                        "Or an array of type xarray.DataArray."
                    )

                if not filepath.is_file():
                    raise FileNotFoundError(
                        "The CSV file that contains fleet composition could not be found."
                    )

                if filepath.suffix != ".csv":
                    raise TypeError(
                        "A CSV file is expected to build the fleet composition."
                    )

                self.fleet = build_fleet_array(filepath, self.scope)

        else:
            self.fleet = None

        array = array.sel(
            powertrain=self.scope["powertrain"],
            year=self.scope["year"],
            size=self.scope["size"],
        )

        self.array = array.stack(desired=["size", "powertrain", "year"])

        self.compliant_vehicles = 1 - array.sel(parameter="has_low_range")

        # store some important specs for inventory documentation
        self.specs = array.sel(
            parameter=[
                "combustion power",
                "electric power",
                "combustion power share",
                "lifetime kilometers",
                "kilometers per year",
                "range",
                "TtW efficiency",
                "TtW energy",
                "fuel cell system efficiency",
                "electric energy stored",
                "oxidation energy stored",
                "energy battery mass",
                "curb mass",
                "driving mass",
            ]
        )

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

        if background_configuration is not None:
            self.background_configuration = background_configuration
        else:
            self.background_configuration = {}

        if "energy storage" not in self.background_configuration:
            self.background_configuration["energy storage"] = {
                "electric": {"type": "NMC-622", "origin": "CN"}
            }
        else:
            if "electric" not in self.background_configuration["energy storage"]:
                self.background_configuration["energy storage"]["electric"] = {
                    "type": "NMC-622",
                    "origin": "CN",
                }
            else:
                if (
                    "origin"
                    not in self.background_configuration["energy storage"]["electric"]
                ):
                    self.background_configuration["energy storage"]["electric"][
                        "origin"
                    ] = "CN"
                if (
                    "type"
                    not in self.background_configuration["energy storage"]["electric"]
                ):
                    self.background_configuration["energy storage"]["electric"][
                        "type"
                    ] = "NMC-622"

        self.inputs = self.get_dict_input()
        self.bs = BackgroundSystemModel()
        self.country = self.get_country_of_use()
        self.add_additional_activities()
        self.rev_inputs = self.get_rev_dict_input()

        with open(DATA_DIR / "fuel_specs.yaml", "r", encoding="utf-8") as stream:
            self.fuel_specs = yaml.safe_load(stream)

        with open(DATA_DIR / "elec_tech_map.yaml", "r", encoding="utf-8") as stream:
            self.elec_map = yaml.safe_load(stream)
            self.elec_map = {k: tuple(v) for k, v in self.elec_map.items()}

        self.elec_tech = list(self.elec_map.keys())
        self.A = self.get_A_matrix()
        self.mix = self.define_electricity_mix_for_fuel_prep()
        self.create_fuel_dictionary()
        self.fuel_blends = {}
        self.define_fuel_blends()
        self.set_actual_range()

        if "direct air capture" in self.background_configuration:
            if "heat source" in self.background_configuration["direct air capture"]:
                heat_source = self.background_configuration["direct air capture"][
                    "heat source"
                ]
                self.select_heat_supplier(heat_source)
            else:
                self.background_configuration["direct air capture"][
                    "heat source"
                ] = "heat pump"
                self.select_heat_supplier("heat pump")
        else:
            self.background_configuration["direct air capture"] = {
                "heat source": "heat pump"
            }
            self.select_heat_supplier("heat pump")

        with open(
            DATA_DIR / "exhaust_and_noise_flows.yaml", "r", encoding="utf-8"
        ) as stream:
            flows = yaml.safe_load(stream)["exhaust"]

            d_comp = {
                "urban": "urban air close to ground",
                "suburban": "non-urban air or from high stacks",
                "rural": "low population density, long-term",
            }

            self.map_fuel_emissions = {
                (v, ("air", d_comp[comp]), "kilogram"): f"{k} direct emissions, {comp}"
                for k, v in flows.items()
                for comp in ["urban", "suburban", "rural"]
            }

        self.index_emissions = [self.inputs[i] for i in self.map_fuel_emissions.keys()]

        self.map_noise_emissions = {
            (
                f"noise, octave {i}, day time, {comp}",
                (f"octave {i}", "day time", comp),
                "joule",
            ): f"noise, octave {i}, day time, {comp}"
            for i in range(1, 9)
            for comp in ["urban", "suburban", "rural"]
        }

        self.index_noise = [self.inputs[i] for i in self.map_noise_emissions.keys()]
        self.list_cat, self.split_indices = self.get_split_indices()
        self.method = method

        if self.method == "recipe":
            self.method_type = method_type
        else:
            self.method_type = "midpoint"

        self.impact_categories = self.get_dict_impact_categories()

        # Load the B matrix
        self.B = None

    def get_results_table(self, split, sensitivity=False):
        """
        Format an xarray.DataArray array to receive the results.

        :param sensitivity:
        :param split: "components" or "impact categories". Split by impact categories only applicable when "endpoint" level is applied.
        :return: xarrray.DataArray
        """

        if split == "components":
            cat = [
                "direct - exhaust",
                "direct - non-exhaust",
                "energy chain",
                "maintenance",
                "glider",
                "EoL",
                "powertrain",
                "energy storage",
                "road",
            ]

        dict_impact_cat = list(self.impact_categories.keys())

        sizes = self.scope["size"]

        if isinstance(self.fleet, xr.core.dataarray.DataArray):
            sizes += ["fleet average"]

        if not sensitivity:

            response = xr.DataArray(
                np.zeros(
                    (
                        self.B.shape[1],
                        len(sizes),
                        len(self.scope["powertrain"]),
                        len(self.scope["year"]),
                        len(cat),
                        self.iterations,
                    )
                ),
                coords=[
                    dict_impact_cat,
                    sizes,
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
            params = [a for a in self.array.value.values]
            response = xr.DataArray(
                np.zeros(
                    (
                        self.B.shape[1],
                        len(sizes),
                        len(self.scope["powertrain"]),
                        len(self.scope["year"]),
                        self.iterations,
                    )
                ),
                coords=[
                    dict_impact_cat,
                    sizes,
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

        with open(filepath, encoding="utf-8") as f:
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

        d["direct - exhaust"] = []
        d["direct - exhaust"].append(
            self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")]
        )
        d["direct - exhaust"].append(
            self.inputs[("Carbon dioxide, non-fossil", ("air",), "kilogram")]
        )
        d["direct - exhaust"].append(
            self.inputs[("Methane, fossil", ("air",), "kilogram")]
        )

        d["direct - exhaust"].extend(self.index_emissions)
        d["direct - exhaust"].extend(self.index_noise)

        l.append(d["direct - exhaust"])

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

        # idx for an input that has no burden
        # oxygen in this case
        extra_idx = [j for i, j in self.inputs.items() if i[0].lower() == "oxygen"][0]

        list_ind = [d[x] for x in d]
        maxLen = max(map(len, list_ind))
        for row in list_ind:
            while len(row) < maxLen:
                row.append(extra_idx)
        return list(d.keys()), list_ind

    def calculate_impacts(self, split="components", sensitivity=False):

        self.B = self.get_B_matrix()

        # Prepare an array to store the results
        results = self.get_results_table(split, sensitivity=sensitivity)

        # Create electricity and fuel market datasets
        self.create_electricity_market_for_fuel_prep()

        # Create electricity market dataset for battery production
        self.create_electricity_market_for_battery_production()

        # Fill in the A matrix with car parameters
        self.set_inputs_in_A_matrix(self.array.values)

        # Add rows for fleet vehicles, if any
        if isinstance(self.fleet, xr.DataArray):
            self.build_fleet_vehicles()

            # Update number of cars
            self.number_of_cars += len(self.scope["year"]) * len(
                self.scope["powertrain"]
            )

            # Update dictionary
            self.get_rev_dict_input()

            # Update B matrix
            self.B = self.get_B_matrix()

        new_arr = np.float32(
            np.zeros((self.A.shape[1], self.B.shape[1], len(self.scope["year"])))
        )

        f_vector = np.zeros((np.shape(self.A)[1]))

        # Collect indices of activities contributing to the first level
        ind_cars = [
            self.inputs[i] for i in self.inputs if "transport, passenger" in i[0]
        ]
        arr = self.A[0, : -self.number_of_cars, ind_cars].sum(axis=0)
        ind = np.nonzero(arr)[0]

        if self.scenario != "static":
            B = self.B.interp(
                year=self.scope["year"], kwargs={"fill_value": "extrapolate"}
            ).values
        else:
            B = self.B[0].values

        for a in ind:
            f_vector[:] = 0
            f_vector[a] = 1
            X = np.float32(
                sparse.linalg.spsolve(sparse.csr_matrix(self.A[0]), f_vector.T)
            )

            if self.scenario == "static":
                new_arr[a] = np.float32(X * B).sum(axis=-1).T[..., None]
            else:
                new_arr[a] = np.float32(X * B).sum(axis=-1).T

        shape = (
            self.iterations,
            len(self.scope["size"]),
            len(self.scope["powertrain"]),
            len(self.scope["year"]),
            self.A.shape[1],
        )

        arr = (
            self.A[:, :, -self.number_of_cars :].transpose(0, 2, 1).reshape(shape)
            * new_arr.transpose(1, 2, 0)[:, None, None, None, ...]
            * -1
        )
        arr = arr[..., self.split_indices].sum(axis=-1)

        if sensitivity:
            results[...] = arr.transpose(0, 2, 3, 4, 5, 1).sum(axis=-2)
            results /= results.sel(parameter="reference")
        else:
            results[...] = arr.transpose(0, 2, 3, 4, 5, 1)

        # If the FU is in passenger-km, we normalize the results by the number of passengers
        if self.scope["fu"]["unit"] == "vkm":
            load_factor = 1
        else:

            load_factor = np.resize(
                self.array[self.array_inputs["average passengers"]].values,
                (
                    1,
                    len(self.scope["size"]),
                    len(self.scope["powertrain"]),
                    len(self.scope["year"]),
                    1,
                    1,
                ),
            )

        if sensitivity:
            return (
                results.astype("float32")
                / load_factor
                * int(self.scope["fu"]["quantity"])
                * self.compliant_vehicles.values[None, ...]
            )
        else:
            return (
                results.astype("float32")
                / load_factor
                * int(self.scope["fu"]["quantity"])
                * self.compliant_vehicles.values[None, :, :, :, None, ...]
            )

    def add_additional_activities(self):
        # Add as many rows and columns as cars to consider
        # Also add additional columns and rows for electricity markets
        # for fuel preparation and energy battery production

        maximum = max(self.inputs.values())

        for y in self.scope["year"]:

            if {"ICEV-p", "HEV-p", "PHEV-p"}.intersection(
                set(self.scope["powertrain"])
            ):
                maximum += 1
                self.inputs[
                    (
                        "fuel supply for gasoline vehicles, " + str(y),
                        self.country,
                        "kilogram",
                        "fuel",
                    )
                ] = maximum

            if {"ICEV-d", "HEV-d", "PHEV-d"}.intersection(
                set(self.scope["powertrain"])
            ):
                maximum += 1
                self.inputs[
                    (
                        "fuel supply for diesel vehicles, " + str(y),
                        self.country,
                        "kilogram",
                        "fuel",
                    )
                ] = maximum

            if {"ICEV-g"}.intersection(set(self.scope["powertrain"])):
                maximum += 1
                self.inputs[
                    (
                        "fuel supply for gas vehicles, " + str(y),
                        self.country,
                        "kilogram",
                        "fuel",
                    )
                ] = maximum

            if {"FCEV"}.intersection(set(self.scope["powertrain"])):
                maximum += 1
                self.inputs[
                    (
                        "fuel supply for hydrogen vehicles, " + str(y),
                        self.country,
                        "kilogram",
                        "fuel",
                    )
                ] = maximum

            if {"BEV", "PHEV-p", "PHEV-d"}.intersection(set(self.scope["powertrain"])):
                maximum += 1
                self.inputs[
                    (
                        "electricity supply for electric vehicles, " + str(y),
                        self.country,
                        "kilowatt hour",
                        "electricity, low voltage, for battery electric vehicles",
                    )
                ] = maximum

            maximum += 1
            self.inputs[
                (
                    "electricity market for fuel preparation, " + str(y),
                    self.country,
                    "kilowatt hour",
                    "electricity, low voltage",
                )
            ] = maximum

            maximum += 1
            self.inputs[
                (
                    "electricity market for energy storage production, " + str(y),
                    self.background_configuration["energy storage"]["electric"][
                        "origin"
                    ],
                    "kilowatt hour",
                    "electricity, low voltage, for energy storage production",
                )
            ] = maximum

        for s in self.scope["size"]:
            for pt in self.scope["powertrain"]:
                for y in self.scope["year"]:
                    maximum += 1

                    if y < 1993:
                        euro_class = "EURO-0"
                    elif 1993 <= y < 1997:
                        euro_class = "EURO-1"
                    elif 1997 <= y < 2001:
                        euro_class = "EURO-2"
                    elif 2001 <= y < 2006:
                        euro_class = "EURO-3"
                    elif 2006 <= y < 2011:
                        euro_class = "EURO-4"
                    elif 2011 <= y < 2015:
                        euro_class = "EURO-5"
                    elif 2015 <= y < 2017:
                        euro_class = "EURO-6ab"
                    elif 2017 <= y < 2019:
                        euro_class = "EURO-6c"
                    elif 2019 <= y <= 2020:
                        euro_class = "EURO-6d-TEMP"
                    else:
                        euro_class = "EURO-6d"

                    if self.scope["fu"]["unit"] == "vkm":
                        unit = "kilometer"
                    else:
                        unit = "person kilometer"

                    if pt == "BEV":
                        chemistry = self.background_configuration["energy storage"][
                            "electric"
                        ]["type"]
                        name = f"transport, passenger car, {pt}, {chemistry} battery, {s}, {y}"
                        ref = "transport, passenger car"

                    elif pt == "FCEV":
                        name = f"transport, passenger car, {pt}, {s}, {y}"
                        ref = "transport, passenger car"

                    else:
                        name = f"transport, passenger car, {pt}, {s}, {y}, {euro_class}"
                        ref = f"transport, passenger car, {euro_class}"

                    self.inputs[
                        (name, self.background_configuration["country"], unit, ref)
                    ] = maximum

    def add_additional_activities_for_export(self):
        # Add as many rows and columns as cars to consider
        # Also add additional columns and rows for electricity markets
        # for fuel preparation and energy battery production

        maximum = max(self.inputs.values())

        for s in self.scope["size"]:
            for pt in self.scope["powertrain"]:
                for y in self.scope["year"]:
                    maximum += 1

                    if y < 1993:
                        euro_class = "EURO-0"
                    elif 1993 <= y < 1997:
                        euro_class = "EURO-1"
                    elif 1997 <= y < 2001:
                        euro_class = "EURO-2"
                    elif 2001 <= y < 2006:
                        euro_class = "EURO-3"
                    elif 2006 <= y < 2011:
                        euro_class = "EURO-4"
                    elif 2011 <= y < 2015:
                        euro_class = "EURO-5"
                    elif 2015 <= y < 2017:
                        euro_class = "EURO-6ab"
                    elif 2017 <= y < 2019:
                        euro_class = "EURO-6c"
                    elif 2019 <= y < 2020:
                        euro_class = "EURO-6d-TEMP"
                    else:
                        euro_class = "EURO-6d"

                    if pt == "BEV":
                        chemistry = self.background_configuration["energy storage"][
                            "electric"
                        ]["type"]
                        name = f"Passenger car, {pt}, {chemistry} battery, {s}, {y}"
                        ref = "Passenger car"

                    elif pt == "FCEV":
                        name = f"Passenger car, {pt}, {s}, {y}"
                        ref = "Passenger car"

                    else:
                        name = f"Passenger car, {pt}, {s}, {y}, {euro_class}"
                        ref = f"Passenger car, {euro_class}"

                    self.inputs[
                        (name, self.background_configuration["country"], "unit", ref)
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

        # build matrix A from coordinates
        A_coords = np.genfromtxt(filepath, delimiter=";")
        I = A_coords[:, 0].astype(int)
        J = A_coords[:, 1].astype(int)
        initial_A = sparse.csr_matrix((A_coords[:, 2], (I, J))).toarray()

        new_A = np.identity(len(self.inputs)).astype("float32")

        new_A[0 : np.shape(initial_A)[0], 0 : np.shape(initial_A)[0]] = initial_A

        # Resize the matrix to fit the number of iterations in `array`

        new_A = np.resize(
            new_A,
            (self.array.shape[1], new_A.shape[0], new_A.shape[1]),
        )
        return new_A

    def build_fleet_vehicles(self):

        # additional cars
        n_cars = len(self.scope["year"]) * len(self.scope["powertrain"]) + len(
            self.scope["year"]
        )
        self.A = np.pad(self.A, ((0, 0), (0, n_cars), (0, n_cars)))
        maximum = max(self.inputs.values())

        for pt in self.scope["powertrain"]:

            for y in self.scope["year"]:

                # share of the powertrain that year, all sizes
                share_pt = self.fleet.sel(powertrain=pt, variable=y).sum().values

                name = "transport, passenger car, fleet average, " + pt + ", " + str(y)

                maximum += 1

                if self.scope["fu"]["unit"] == "vkm":
                    unit = "kilometer"
                else:
                    unit = "person kilometer"

                self.inputs[
                    (
                        name,
                        self.background_configuration["country"],
                        unit,
                        "transport, passenger car, fleet average",
                    )
                ] = maximum

                self.A[:, maximum, maximum] = 1

                if share_pt > 0:
                    for s in self.fleet.coords["size"].values:
                        for vin_year in range(min(self.scope["year"]), y + 1):
                            if vin_year in self.fleet.vintage_year:
                                fleet_share = (
                                    self.fleet.sel(
                                        powertrain=pt,
                                        vintage_year=vin_year,
                                        size=s,
                                        variable=y,
                                    )
                                    .sum()
                                    .values
                                    / share_pt
                                )

                                if fleet_share > 0:

                                    car_index = [
                                        self.inputs[i]
                                        for i in self.inputs
                                        if all(
                                            [
                                                item in i[0]
                                                for item in [
                                                    pt,
                                                    str(vin_year),
                                                    s,
                                                    "transport",
                                                ]
                                            ]
                                        )
                                    ][0]

                                    car_inputs = (
                                        self.A[:, : car_index - 1, car_index]
                                        * fleet_share
                                    )

                                    self.A[:, : car_index - 1, maximum] += car_inputs

                    # Fuel and electricity supply must be from the fleet year, not the vintage year

                    d_map_fuel = {
                        "ICEV-p": "gasoline",
                        "ICEV-d": "diesel",
                        "ICEV-g": "gas",
                        "HEV-p": "gasoline",
                        "HEV-d": "diesel",
                        "PHEV-p": "gasoline",
                        "PHEV-d": "diesel",
                        "BEV": "electric",
                        "FCEV": "hydrogen",
                    }

                    ind_supply = [
                        self.inputs[i]
                        for i in self.inputs
                        if "supply for " + d_map_fuel[pt] + " vehicles, " in i[0]
                    ]
                    amount_fuel = self.A[:, ind_supply, maximum].sum(axis=1)

                    # zero out initial fuel inputs
                    self.A[:, ind_supply, maximum] = 0

                    # set saved amount to current fuel supply provider
                    current_provider = [
                        self.inputs[i]
                        for i in self.inputs
                        if "supply for " + d_map_fuel[pt] + " vehicles, " + str(y)
                        in i[0]
                    ]
                    self.A[:, current_provider, maximum] = amount_fuel

                    if pt in ["PHEV-p", "PHEV-d"]:
                        ind_supply = [
                            self.inputs[i]
                            for i in self.inputs
                            if "supply for electric vehicles, " in i[0]
                        ]
                        amount_fuel = self.A[:, ind_supply, maximum].sum(axis=1)

                        # zero out initial fuel inputs
                        self.A[:, ind_supply, maximum] = 0

                        # set saved amount to current fuel supply provider
                        current_provider = [
                            self.inputs[i]
                            for i in self.inputs
                            if "supply for electric vehicles, " + str(y) in i[0]
                        ]
                        self.A[:, current_provider, maximum] = amount_fuel

        # We also want to produce a fleet average vehicle, with all powertrain types

        for y in self.scope["year"]:

            # share of that year, all sizes and powertrains
            share_pt = self.fleet.sel(variable=y).sum().values

            name = "transport, passenger car, fleet average, all powertrains, " + str(y)

            maximum += 1

            if self.scope["fu"]["unit"] == "vkm":
                unit = "kilometer"
            else:
                unit = "person kilometer"

            self.inputs[
                (
                    name,
                    self.background_configuration["country"],
                    unit,
                    "transport, passenger car, fleet average",
                )
            ] = maximum

            self.A[:, maximum, maximum] = 1

            if share_pt > 0:
                for pt in self.fleet.coords["powertrain"].values:
                    for s in self.fleet.coords["size"].values:
                        for vin_year in range(min(self.scope["year"]), y + 1):

                            if vin_year in self.fleet.vintage_year:

                                fleet_share = (
                                    self.fleet.sel(
                                        powertrain=pt,
                                        vintage_year=vin_year,
                                        size=s,
                                        variable=y,
                                    )
                                    .sum()
                                    .values
                                    / share_pt
                                )

                                if fleet_share > 0:
                                    car_index = [
                                        self.inputs[i]
                                        for i in self.inputs
                                        if all(
                                            [
                                                item in i[0]
                                                for item in [
                                                    pt,
                                                    str(vin_year),
                                                    s,
                                                    "transport",
                                                ]
                                            ]
                                        )
                                    ][0]

                                    car_inputs = (
                                        self.A[:, : car_index - 1, car_index]
                                        * fleet_share
                                    )

                                    self.A[:, : car_index - 1, maximum] += car_inputs

            # Fuel and electricity supply must be from the fleet year, not the vintage year
            d_map_fuel = {
                "ICEV-p": "gasoline",
                "ICEV-d": "diesel",
                "ICEV-g": "gas",
                "HEV-p": "gasoline",
                "HEV-d": "diesel",
                "PHEV-p": "gasoline",
                "PHEV-d": "diesel",
                "BEV": "electric",
                "FCEV": "hydrogen",
            }

            for fuel_type in set(d_map_fuel.values()):

                ind_supply = [
                    self.inputs[i]
                    for i in self.inputs
                    if "supply for " + fuel_type + " vehicles, " in i[0]
                ]
                amount_fuel = self.A[:, ind_supply, maximum].sum(axis=1)

                if amount_fuel > 0:

                    # zero out initial fuel inputs
                    self.A[:, ind_supply, maximum] = 0

                    # set saved amount to current fuel supply provider
                    current_provider = [
                        self.inputs[i]
                        for i in self.inputs
                        if "supply for " + fuel_type + " vehicles, " + str(y) in i[0]
                    ]

                    self.A[:, current_provider, maximum] = amount_fuel

    def get_B_matrix(self):
        """
        Load the B matrix. The B matrix contains impact assessment figures for a give impact assessment method,
        per unit of activity. Its length column-wise equals the length of the A matrix row-wise.
        Its length row-wise equals the number of impact assessment methods.

        :return: an array with impact values per unit of activity for each method.
        :rtype: numpy.ndarray

        """

        if self.method == "recipe":
            if self.method_type == "midpoint":
                list_file_names = glob.glob(
                    str(IAM_FILES_DIR)
                    + "/*recipe_midpoint*{}*.csv".format(self.scenario)
                )
                list_file_names = sorted(list_file_names)
                B = np.zeros((len(list_file_names), 22, len(self.inputs)))
            elif self.method_type == "endpoint":
                list_file_names = glob.glob(
                    str(IAM_FILES_DIR)
                    + "/*recipe_endpoint*{}*.csv".format(self.scenario)
                )
                list_file_names = sorted(list_file_names)
                B = np.zeros((len(list_file_names), 6, len(self.inputs)))
            else:
                raise TypeError(
                    "The LCIA method type should be either 'midpoint' or 'endpoint'."
                )

        else:
            list_file_names = glob.glob(
                str(IAM_FILES_DIR) + "/*ilcd*{}*.csv".format(self.scenario)
            )
            list_file_names = sorted(list_file_names)
            B = np.zeros((len(list_file_names), 19, len(self.inputs)))

        for f, filepath in enumerate(list_file_names):
            initial_B = np.genfromtxt(filepath, delimiter=";")
            new_B = np.zeros(
                (
                    np.shape(initial_B)[0],
                    len(self.inputs),
                )
            )
            new_B[0 : np.shape(initial_B)[0], 0 : np.shape(initial_B)[1]] = initial_B
            B[f, :, :] = new_B

        list_impact_categories = list(self.impact_categories.keys())

        if self.scenario != "static":
            response = xr.DataArray(
                B,
                coords=[
                    [2005, 2010, 2020, 2030, 2040, 2050],
                    list_impact_categories,
                    list(self.inputs.keys()),
                ],
                dims=["year", "category", "activity"],
            )
        else:
            response = xr.DataArray(
                B,
                coords=[[2020], list_impact_categories, list(self.inputs.keys())],
                dims=["year", "category", "activity"],
            )

        return response

    @staticmethod
    def get_dict_input():
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
        with open(filepath, encoding="utf-8") as f:
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
                                   'human noise',
                                   'primary energy, non-renewable',
                                   'primary energy, renewable']
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

        with open(filepath, encoding="utf-8") as f:
            input_dict = csv.reader(f, delimiter=";")
            for row in input_dict:
                if row[0] == self.method and row[3] == self.method_type:
                    csv_dict[row[2]] = {
                        "method": row[1],
                        "category": row[2],
                        "type": row[3],
                        "abbreviation": row[4],
                        "unit": row[5],
                        "source": row[6],
                    }

        return csv_dict

    def get_rev_dict_input(self):
        """
        Reverse the self.inputs dictionary.

        :return: reversed dictionary
        :rtype: dict
        """
        return {v: k for k, v in self.inputs.items()}

    def get_index_vehicle_from_array(
        self, items_to_look_for, items_to_look_for_also=None, method="or"
    ):
        """
        Return list of row/column indices of self.array of labels that contain the string defined in `items_to_look_for`.

        :param items_to_look_for_also:
        :param method:
        :param items_to_look_for: string to search for
        :return: list
        """
        if not isinstance(items_to_look_for, list):
            items_to_look_for = [items_to_look_for]

        if not items_to_look_for_also is None:
            if not isinstance(items_to_look_for_also, list):
                items_to_look_for_also = [items_to_look_for_also]

        list_vehicles = self.array.desired.values

        if method == "or":
            return [
                c
                for c, v in enumerate(list_vehicles)
                if set(items_to_look_for).intersection(v)
            ]

        if method == "and":
            return [
                c
                for c, v in enumerate(list_vehicles)
                if set(items_to_look_for).intersection(v)
                and set(items_to_look_for_also).intersection(v)
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

    def resize_A_matrix_for_export(self):
        """Removes some vehicles if they do not comply with
        requirements in terms of range"""

        indices_to_remove = []

        for i in self.inputs:
            if (
                "passenger car, " in i[0].lower()
                and "fleet average" not in i[0].lower()
                and "market" not in i[0].lower()
                and "used" not in i[0].lower()
            ):

                split_name = [x.strip() for x in i[0].split(", ")]

                if split_name[0] == "transport":
                    if len(split_name) == 5:
                        _, _, pt, size, year = split_name
                    else:
                        if split_name[2] == "BEV":
                            _, _, pt, _, size, year = split_name
                        else:
                            _, _, pt, size, year, _ = split_name
                else:
                    if len(split_name) == 5:
                        if split_name[1] == "BEV":
                            # electric vehicle
                            _, pt, _, size, year = i[0].split(", ")
                        else:
                            # combustion vehicles
                            _, pt, size, year, _ = i[0].split(", ")
                    else:
                        # FCEV vehicle
                        _, pt, size, year = i[0].split(", ")

                if (
                    self.compliant_vehicles.sel(
                        powertrain=pt, size=size, year=int(year)
                    )
                    == 0
                ):
                    indices_to_remove.append(self.inputs[i])
                    self.rev_inputs.pop(self.inputs[i])

        indices_to_preserve = [
            i for i in range(self.A.shape[1]) if i not in indices_to_remove
        ]

        self.A = self.A[
            np.ix_(range(self.A.shape[0]), indices_to_preserve, indices_to_preserve)
        ]

        self.rev_inputs = dict(enumerate(self.rev_inputs.values()))

    def export_lci(
        self,
        presamples=False,
        ecoinvent_version="3.8",
        db_name="carculator db",
        create_vehicle_datasets=True,
    ):
        """
        Export the inventory as a dictionary. Also return a list of arrays that contain pre-sampled random values if
        :param db_name:
        :meth:`stochastic` of :class:`CarModel` class has been called.

        :param presamples: boolean.
        :param ecoinvent_version: str. "3.5", "3.6" or "uvek"
        :param create_vehicle_datasets: bool. Whether vehicles datasets (as structured in ecoinvent) should be created too.
        :return: inventory, and optionally, list of arrays containing pre-sampled values.
        :rtype: list
        """

        self.inputs = self.get_dict_input()
        self.bs = BackgroundSystemModel()
        self.country = self.get_country_of_use()
        self.add_additional_activities()
        self.rev_inputs = self.get_rev_dict_input()
        self.A = self.get_A_matrix()

        if "direct air capture" in self.background_configuration:
            if "heat source" in self.background_configuration["direct air capture"]:
                heat_source = self.background_configuration["direct air capture"][
                    "heat source"
                ]
                self.select_heat_supplier(heat_source)
            else:
                self.background_configuration["direct air capture"][
                    "heat source"
                ] = "heat pump"
                self.select_heat_supplier("heat pump")
        else:
            self.background_configuration["direct air capture"] = {
                "heat source": "heat pump"
            }
            self.select_heat_supplier("heat pump")

        if create_vehicle_datasets:

            # add vehicles datasets
            self.add_additional_activities_for_export()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # resize A matrix
            self.A = self.get_A_matrix()

            if "direct air capture" in self.background_configuration:
                if "heat source" in self.background_configuration["direct air capture"]:
                    heat_source = self.background_configuration["direct air capture"][
                        "heat source"
                    ]
                    self.select_heat_supplier(heat_source)
                else:
                    self.background_configuration["direct air capture"][
                        "heat source"
                    ] = "heat pump"
                    self.select_heat_supplier("heat pump")
            else:
                self.background_configuration["direct air capture"] = {
                    "heat source": "heat pump"
                }
                self.select_heat_supplier("heat pump")

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            self.set_inputs_in_A_matrix_for_export(self.array.values)

        else:

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            self.set_inputs_in_A_matrix(self.array.values)

        # Add rows for fleet vehicles, if any
        if isinstance(self.fleet, xr.DataArray):
            print("Building fleet average vehicles...")
            self.build_fleet_vehicles()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # Update number of cars
            self.number_of_cars += len(self.scope["year"]) * len(
                self.scope["powertrain"]
            )

        # Remove vehicles not compliant or available
        self.resize_A_matrix_for_export()

        if presamples:
            lci, array = ExportInventory(
                self.A, self.rev_inputs, db_name=db_name
            ).write_lci(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=self.specs,
            )
            return lci, array
        else:
            lci = ExportInventory(self.A, self.rev_inputs, db_name=db_name).write_lci(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=self.specs,
            )
            return lci

    def export_lci_to_bw(
        self,
        presamples=False,
        ecoinvent_version="3.8",
        db_name="carculator db",
        create_vehicle_datasets=True,
    ):
        """
        Export the inventory as a `brightway2` bw2io.importers.base_lci.LCIImporter object
        with the inventory in the `data` attribute.

        .. code-block:: python

            # get the inventory
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

        self.inputs = self.get_dict_input()
        self.bs = BackgroundSystemModel()
        self.country = self.get_country_of_use()
        self.add_additional_activities()
        self.rev_inputs = self.get_rev_dict_input()
        self.A = self.get_A_matrix()

        if "direct air capture" in self.background_configuration:
            if "heat source" in self.background_configuration["direct air capture"]:
                heat_source = self.background_configuration["direct air capture"][
                    "heat source"
                ]
                self.select_heat_supplier(heat_source)
            else:
                self.background_configuration["direct air capture"][
                    "heat source"
                ] = "heat pump"
                self.select_heat_supplier("heat pump")
        else:
            self.background_configuration["direct air capture"] = {
                "heat source": "heat pump"
            }
            self.select_heat_supplier("heat pump")

        if create_vehicle_datasets:

            # add vehicles datasets
            self.add_additional_activities_for_export()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # resize A matrix
            self.A = self.get_A_matrix()

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            if "direct air capture" in self.background_configuration:
                if "heat source" in self.background_configuration["direct air capture"]:
                    heat_source = self.background_configuration["direct air capture"][
                        "heat source"
                    ]
                    self.select_heat_supplier(heat_source)
                else:
                    self.background_configuration["direct air capture"][
                        "heat source"
                    ] = "heat pump"
                    self.select_heat_supplier("heat pump")
            else:
                self.background_configuration["direct air capture"] = {
                    "heat source": "heat pump"
                }
                self.select_heat_supplier("heat pump")

            self.set_inputs_in_A_matrix_for_export(self.array.values)

        else:

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            self.set_inputs_in_A_matrix(self.array.values)

        # Add rows for fleet vehicles, if any
        if isinstance(self.fleet, xr.core.dataarray.DataArray):
            print("Building fleet average vehicles...")
            self.build_fleet_vehicles()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # Update number of cars
            self.number_of_cars += len(self.scope["year"]) * len(
                self.scope["powertrain"]
            )

        # Remove vehicles not compliant or available
        self.resize_A_matrix_for_export()

        if presamples:
            lci, array = ExportInventory(
                self.A, self.rev_inputs, db_name=db_name
            ).write_lci_to_bw(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=self.specs,
            )
            return lci, array
        else:

            lci = ExportInventory(
                self.A, self.rev_inputs, db_name=db_name
            ).write_lci_to_bw(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=self.specs,
            )

            return lci

    def export_lci_to_excel(
        self,
        directory=None,
        ecoinvent_version="3.8",
        software_compatibility="brightway2",
        filename=None,
        create_vehicle_datasets=True,
        export_format="file",
    ):
        """
        Export the inventory as an Excel file (if the destination software is Brightway2) or a CSV file (if the destination software is Simapro) file.
        Also return the file path where the file is stored.

        :param filename:
        :param create_vehicle_datasets:
        :param export_format:
        :param directory: directory where to save the file.
        :type directory: str
        :param ecoinvent_version: "3.8", "3.7", "3.6", "3.5" or "uvek"
        :param software_compatibility: "brightway2" or "simapro"
        :return: file path where the file is stored.
        :rtype: str
        """

        if software_compatibility not in ("brightway2", "simapro"):
            raise NameError(
                "The destination software argument is not valid. Choose between 'brightway2' or 'simapro'."
            )

        # Simapro inventory only for ecoinvent 3.5 or UVEK
        if software_compatibility == "simapro":
            if ecoinvent_version == "3.8":
                print(
                    "Simapro-compatible inventory export is only available for ecoinvent 3.5, 3.6, 3.7 or UVEK."
                )
                return

        self.inputs = self.get_dict_input()
        self.bs = BackgroundSystemModel()
        self.country = self.get_country_of_use()
        self.add_additional_activities()
        self.rev_inputs = self.get_rev_dict_input()
        self.A = self.get_A_matrix()

        if "direct air capture" in self.background_configuration:
            if "heat source" in self.background_configuration["direct air capture"]:
                heat_source = self.background_configuration["direct air capture"][
                    "heat source"
                ]
                self.select_heat_supplier(heat_source)
            else:
                self.background_configuration["direct air capture"][
                    "heat source"
                ] = "heat pump"
                self.select_heat_supplier("heat pump")
        else:
            self.background_configuration["direct air capture"] = {
                "heat source": "heat pump"
            }
            self.select_heat_supplier("heat pump")

        if create_vehicle_datasets:

            # add vehicles datasets
            self.add_additional_activities_for_export()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # resize A matrix
            self.A = self.get_A_matrix()

            if "direct air capture" in self.background_configuration:
                if "heat source" in self.background_configuration["direct air capture"]:
                    heat_source = self.background_configuration["direct air capture"][
                        "heat source"
                    ]
                    self.select_heat_supplier(heat_source)
                else:
                    self.background_configuration["direct air capture"][
                        "heat source"
                    ] = "heat pump"
                    self.select_heat_supplier("heat pump")
            else:
                self.background_configuration["direct air capture"] = {
                    "heat source": "heat pump"
                }
                self.select_heat_supplier("heat pump")

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            self.set_inputs_in_A_matrix_for_export(self.array.values)

        else:

            # Create electricity and fuel market datasets
            self.create_electricity_market_for_fuel_prep()

            # Create electricity market dataset for battery production
            self.create_electricity_market_for_battery_production()

            # Create fuel markets
            self.fuel_blends = {}
            self.define_fuel_blends()
            self.set_actual_range()

            self.set_inputs_in_A_matrix(self.array.values)

        # Add rows for fleet vehicles, if any
        if isinstance(self.fleet, xr.core.dataarray.DataArray):
            print("Building fleet average vehicles...")
            self.build_fleet_vehicles()

            # Update dictionary
            self.rev_inputs = self.get_rev_dict_input()

            # Update number of cars
            self.number_of_cars += len(self.scope["year"]) * len(
                self.scope["powertrain"]
            )

        # Remove vehicles not compliant or available
        self.resize_A_matrix_for_export()

        filepath = ExportInventory(
            self.A, self.rev_inputs, db_name=filename or "carculator db"
        ).write_lci_to_excel(
            directory=directory,
            ecoinvent_version=ecoinvent_version,
            software_compatibility=software_compatibility,
            filename=filename,
            export_format=export_format,
            vehicle_specs=self.specs,
        )
        return filepath

    def get_country_of_use(self):

        if "country" not in self.background_configuration:
            self.background_configuration["country"] = "RER"

        return self.background_configuration["country"]

    def define_electricity_mix_for_fuel_prep(self):
        """
        This function defines a fuel mix based either on user-defined mix, or on default mixes for a given country.
        The mix is calculated as the average mix, weighted by the distribution of annually driven kilometers.
        :return:
        """
        try:
            losses_to_low = float(self.bs.losses[self.country]["LV"])
        except KeyError:
            # If losses for the country are not found, assume EU average
            losses_to_low = float(self.bs.losses["RER"]["LV"])

        if "custom electricity mix" in self.background_configuration:
            # If a special electricity mix is specified, we use it
            mix = self.background_configuration["custom electricity mix"]

            if np.shape(mix)[0] != len(self.scope["year"]):
                raise ValueError(
                    "The number of electricity mixes ({}) must match with the "
                    "number of years ({}).".format(
                        np.shape(mix)[0], len(self.scope["year"])
                    )
                )

            if not np.allclose(np.sum(mix, 1), np.ones(len(self.scope["year"]))):
                print(
                    "The sum of the electricity mix share does "
                    "not equal to 1 for each year."
                )

        else:
            use_year = (
                (
                    self.array.values[self.array_inputs["lifetime kilometers"]]
                    / self.array.values[self.array_inputs["kilometers per year"]]
                )
                .reshape(
                    self.iterations,
                    len(self.scope["powertrain"]),
                    len(self.scope["size"]),
                    len(self.scope["year"]),
                )
                .mean(axis=(0, 1, 2))
            )

            if self.country not in self.bs.electricity_mix.country.values:
                print(
                    f"The electricity mix for {self.country} could not be found."
                    "Average European electricity mix is used instead."
                )
                country = "RER"
            else:
                country = self.country

            mix = [
                self.bs.electricity_mix.sel(
                    country=country,
                    variable=self.elec_tech,
                )
                .interp(
                    year=np.arange(year, year + use_year[y]),
                    kwargs={"fill_value": "extrapolate"},
                )
                .mean(axis=0)
                .values
                if y + use_year[y] <= 2050
                else self.bs.electricity_mix.sel(
                    country=country,
                    variable=self.elec_tech,
                )
                .interp(
                    year=np.arange(year, 2051), kwargs={"fill_value": "extrapolate"}
                )
                .mean(axis=0)
                .values
                for y, year in enumerate(self.scope["year"])
            ]

        return np.clip(mix, 0, 1) / np.clip(mix, 0, 1).sum(axis=1)[:, None]

    def define_renewable_rate_in_mix(self):

        try:
            losses_to_low = float(self.bs.losses[self.country]["LV"])
        except KeyError:
            # If losses for the country are not found, assume EU average
            losses_to_low = float(self.bs.losses["RER"]["LV"])

        category_name = (
            "climate change"
            if self.method == "recipe"
            else "climate change - climate change total"
        )

        if self.method_type != "endpoint":
            if self.scenario != "static":
                year = self.scope["year"]
                co2_intensity_tech = (
                    self.B.sel(
                        category=category_name,
                        activity=[self.elec_map[t] for t in self.elec_tech],
                    )
                    .interp(year=year, kwargs={"fill_value": "extrapolate"})
                    .values
                    * losses_to_low
                ) * 1000
            else:
                year = 2020
                co2_intensity_tech = np.resize(
                    (
                        self.B.sel(
                            category=category_name,
                            activity=[self.elec_map[t] for t in self.elec_tech],
                            year=year,
                        ).values
                        * losses_to_low
                        * 1000
                    ),
                    (len(self.scope["year"]), 21),
                )
        else:
            co2_intensity_tech = np.zeros((len(self.scope["year"]), 21))

        sum_renew = [
            np.sum([self.mix[x][i] for i in [0, 3, 4, 5, 8]])
            for x in range(0, len(self.mix))
        ]

        return sum_renew, co2_intensity_tech

    def create_electricity_market_for_fuel_prep(self):
        """This function fills the electricity market that supplies battery charging operations
        and hydrogen production through electrolysis.
        """

        try:
            losses_to_low = float(self.bs.losses[self.country]["LV"])
        except KeyError:
            # If losses for the country are not found, assume EU average
            losses_to_low = float(self.bs.losses["RER"]["LV"])

        # Fill the electricity markets for battery charging and hydrogen production
        for y, year in enumerate(self.scope["year"]):
            m = np.array(self.mix[y], dtype=object).reshape((-1, len(self.mix[y]), 1))
            # Add electricity technology shares

            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    [self.inputs[self.elec_map[t]] for t in self.elec_tech],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "electricity market for fuel preparation" in i[0]
                    ],
                )
            ] = (
                m * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                6.58e-9 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                1.86e-8 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                3.17e-10 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                8.74e-8 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                (5.4e-8 + 2.99e-9) * -1 * losses_to_low
            )

            # Add SF_6 leakage

            self.A[
                :,
                self.inputs[("Sulfur hexafluoride", ("air",), "kilogram")],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "electricity market for fuel preparation" in i[0]
                ],
            ] = (
                (5.4e-8 + 2.99e-9) * -1 * losses_to_low
            )

    def create_electricity_market_for_battery_production(self):
        """
        This function fills in the column in `self.A` concerned with the electricity mix used for manufacturing battery cells
        :return:
        """

        battery_origin = self.background_configuration["energy storage"]["electric"][
            "origin"
        ]

        if battery_origin != "custom electricity mix":

            try:
                losses_to_low = float(self.bs.losses[battery_origin]["LV"])
            except KeyError:
                losses_to_low = float(self.bs.losses["CN"]["LV"])

            if battery_origin not in self.bs.electricity_mix.country.values:
                print(
                    "The electricity mix for {} could not be found. Average Chinese electricity mix is used for "
                    "battery manufacture instead.".format(self.country)
                )
                battery_origin = "CN"

            mix_battery_manufacturing = (
                self.bs.electricity_mix.sel(
                    country=battery_origin,
                    variable=self.elec_tech,
                )
                .interp(year=self.scope["year"], kwargs={"fill_value": "extrapolate"})
                .values
            )
        else:
            # electricity mix for battery manufacturing same as `custom electricity mix`
            mix_battery_manufacturing = self.mix
            losses_to_low = 1.1

        # Fill the electricity markets for battery production
        for y, year in enumerate(self.scope["year"]):
            m = np.array(mix_battery_manufacturing[y], dtype=object).reshape(
                (-1, 21, 1)
            )

            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    [self.inputs[self.elec_map[t]] for t in self.elec_tech],
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "electricity market for energy storage production" in i[0]
                    ],
                )
            ] = (
                m * losses_to_low * -1
            )

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
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                6.58e-9 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                1.86e-8 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                3.17e-10 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                8.74e-8 * -1 * losses_to_low
            )

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
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                (5.4e-8 + 2.99e-9) * -1 * losses_to_low
            )

            # Add SF_6 leakage

            self.A[
                :,
                self.inputs[("Sulfur hexafluoride", ("air",), "kilogram")],
                [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "electricity market for energy storage production" in i[0]
                ],
            ] = (
                (5.4e-8 + 2.99e-9) * -1 * losses_to_low
            )

    def get_share_biogasoline(self):
        """Returns average share of biogasoline according to historical IEA stats"""
        share_biofuel = np.squeeze(
            np.clip(
                self.bs.biogasoline.sel(country=self.country)
                .interp(year=self.scope["year"], kwargs={"fill_value": "extrapolate"})
                .values,
                0,
                0.2,
            )
        )
        return share_biofuel

    def get_share_biodiesel(self):
        """Returns average share of biodiesel according to historical IEA stats"""
        share_biofuel = np.squeeze(
            np.clip(
                self.bs.biodiesel.sel(country=self.country)
                .interp(year=self.scope["year"], kwargs={"fill_value": "extrapolate"})
                .values,
                0,
                0.2,
            )
        )
        return share_biofuel

    def get_share_biomethane(self):
        """Returns average share of biomethane according to historical IEA stats"""
        share_biofuel = np.squeeze(
            np.clip(
                self.bs.biomethane.sel(country=self.country)
                .interp(year=self.scope["year"], kwargs={"fill_value": "extrapolate"})
                .values,
                0,
                1,
            )
        )
        return share_biofuel

    def get_share_biofuel(self):
        """Returns average share of biofuel in blend according to IAM REMIND"""
        try:
            region = self.bs.region_map[self.country]["RegionCode"]
        except KeyError:
            print(
                "Could not find biofuel share information for the region/country specified. "
                "Will use European biofuel share instead."
            )
            region = "EUR"
        scenario = self.scenario if self.scenario != "static" else "SSP2-Base"

        share_biofuel = (
            self.bs.biofuel.sel(
                region=region,
                value=0,
                fuel_type="Biomass fuel",
                scenario=scenario,
            )
            .interp(year=self.scope["year"], kwargs={"fill_value": "extrapolate"})
            .values
        )
        return share_biofuel

    def find_fuel_shares(self, fuel_type):

        default_fuels = {
            "petrol": {
                "primary": "petrol",
                "secondary": "bioethanol - wheat straw",
                "third": "bioethanol - maize starch",
            },
            "diesel": {
                "primary": "diesel",
                "secondary": "biodiesel - cooking oil",
                "third": "biodiesel - algae",
            },
            "cng": {
                "primary": "cng",
                "secondary": "biogas - sewage sludge",
                "third": "biogas - biowaste",
            },
            "hydrogen": {
                "primary": "electrolysis",
                "secondary": "smr - natural gas",
                "third": "atr - natural gas",
            },
        }

        tertiary = None
        tertiary_share = None

        if "fuel blend" in self.background_configuration:
            if fuel_type in self.background_configuration["fuel blend"]:
                primary = self.background_configuration["fuel blend"][fuel_type][
                    "primary fuel"
                ]["type"]

                if (
                    "secondary fuel"
                    in self.background_configuration["fuel blend"][fuel_type]
                ):
                    secondary = self.background_configuration["fuel blend"][fuel_type][
                        "secondary fuel"
                    ]["type"]
                else:

                    if default_fuels[fuel_type]["secondary"] != primary:
                        secondary = default_fuels[fuel_type]["secondary"]
                        self.background_configuration["fuel blend"][fuel_type][
                            "secondary fuel"
                        ] = {"type": secondary}
                    else:
                        secondary = default_fuels[fuel_type]["third"]
                        self.background_configuration["fuel blend"][fuel_type][
                            "secondary fuel"
                        ] = {"type": secondary}

                if (
                    "tertiary fuel"
                    in self.background_configuration["fuel blend"][fuel_type]
                ):
                    tertiary = self.background_configuration["fuel blend"][fuel_type][
                        "tertiary fuel"
                    ]["type"]

                primary_share = self.background_configuration["fuel blend"][fuel_type][
                    "primary fuel"
                ]["share"]

                if (
                    "share"
                    not in self.background_configuration["fuel blend"][fuel_type][
                        "secondary fuel"
                    ]
                ):
                    secondary_share = 1 - np.asarray(primary_share)
                else:
                    secondary_share = self.background_configuration["fuel blend"][
                        fuel_type
                    ]["secondary fuel"]["share"]

                if tertiary:
                    tertiary_share = self.background_configuration["fuel blend"][
                        fuel_type
                    ]["tertiary fuel"]["share"]

            else:
                primary = default_fuels[fuel_type]["primary"]
                secondary = default_fuels[fuel_type]["secondary"]

                if primary == "electrolysis":
                    secondary_share = np.zeros_like(self.scope["year"])
                else:
                    if fuel_type == "diesel":
                        if self.country in self.bs.biodiesel.country.values:
                            secondary_share = self.get_share_biodiesel()
                        else:
                            secondary_share = self.get_share_biofuel()

                    elif fuel_type == "petrol":
                        if self.country in self.bs.biogasoline.country.values:
                            secondary_share = self.get_share_biogasoline()
                        else:
                            secondary_share = self.get_share_biofuel()

                    elif fuel_type == "cng":
                        if self.country in self.bs.biomethane.country.values:
                            secondary_share = self.get_share_biomethane()
                        else:
                            secondary_share = self.get_share_biofuel()
                    else:
                        secondary_share = self.get_share_biofuel()

                if secondary_share.shape == ():
                    secondary_share = secondary_share.reshape(1)

                primary_share = 1 - np.array(secondary_share, dtype=object)

        else:
            primary = default_fuels[fuel_type]["primary"]
            secondary = default_fuels[fuel_type]["secondary"]

            if primary == "electrolysis":
                secondary_share = np.zeros_like(self.scope["year"])
            else:
                if fuel_type == "diesel":
                    if self.country in self.bs.biodiesel.country.values:
                        secondary_share = self.get_share_biodiesel()
                    else:
                        secondary_share = self.get_share_biofuel()

                elif fuel_type == "petrol":
                    if self.country in self.bs.biogasoline.country.values:
                        secondary_share = self.get_share_biogasoline()
                    else:
                        secondary_share = self.get_share_biofuel()

                elif fuel_type == "cng":
                    if self.country in self.bs.biomethane.country.values:
                        secondary_share = self.get_share_biomethane()
                    else:
                        secondary_share = self.get_share_biofuel()
                else:
                    secondary_share = self.get_share_biofuel()

            if secondary_share.shape == ():
                secondary_share = secondary_share.reshape(1)
            primary_share = 1 - np.array(secondary_share, dtype=object)

        return (
            primary,
            secondary,
            primary_share,
            secondary_share,
            tertiary,
            tertiary_share,
        )

    def set_actual_range(self):
        """
        Set the actual range considering the blend.
        Liquid bio-fuels and synthetic fuels typically have a lower calorific value. Hence, the need to recalculate
        the vehicle range.
        Modifies parameter `range` of `array` in place
        """

        if {"ICEV-p", "HEV-p", "PHEV-p"}.intersection(set(self.scope["powertrain"])):
            for y, year in enumerate(self.scope["year"]):

                share_primary = self.fuel_blends["petrol"]["primary"]["share"][y]
                lhv_primary = self.fuel_blends["petrol"]["primary"]["lhv"]
                share_secondary = self.fuel_blends["petrol"]["secondary"]["share"][y]
                lhv_secondary = self.fuel_blends["petrol"]["secondary"]["lhv"]

                if "tertiary" in self.fuel_blends["petrol"]:
                    share_tertiary = self.fuel_blends["petrol"]["tertiary"]["share"][y]
                    lhv_tertiary = self.fuel_blends["petrol"]["tertiary"]["lhv"]
                else:
                    share_tertiary = 0
                    lhv_tertiary = 0

                index = self.get_index_vehicle_from_array(
                    ["ICEV-p", "HEV-p", "PHEV-p"], year, method="and"
                )

                self.array.values[self.array_inputs["range"], :, index] = (
                    (
                        (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_primary
                            * lhv_primary
                        )
                        + (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_secondary
                            * lhv_secondary
                        )
                        + (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_tertiary
                            * lhv_tertiary
                        )
                        + (
                            self.array.values[
                                self.array_inputs["electric energy stored"], :, index
                            ]
                            * 3.6
                        )
                    )
                    * 1000
                    / self.array.values[self.array_inputs["TtW energy"], :, index]
                )

                self.array.values[self.array_inputs["LHV fuel MJ per kg"], :, index] = (
                    (share_primary * lhv_primary)
                    + (share_secondary * lhv_secondary)
                    + (share_tertiary * lhv_tertiary)
                )

        if {"ICEV-d", "HEV-d", "PHEV-d"}.intersection(set(self.scope["powertrain"])):
            for y, year in enumerate(self.scope["year"]):
                share_primary = self.fuel_blends["diesel"]["primary"]["share"][y]
                lhv_primary = self.fuel_blends["diesel"]["primary"]["lhv"]
                share_secondary = self.fuel_blends["diesel"]["secondary"]["share"][y]
                lhv_secondary = self.fuel_blends["diesel"]["secondary"]["lhv"]

                if "tertiary" in self.fuel_blends["diesel"]:
                    share_tertiary = self.fuel_blends["diesel"]["tertiary"]["share"][y]
                    lhv_tertiary = self.fuel_blends["diesel"]["tertiary"]["lhv"]
                else:
                    share_tertiary = 0
                    lhv_tertiary = 0

                index = self.get_index_vehicle_from_array(
                    ["ICEV-d", "PHEV-d", "HEV-d"], year, method="and"
                )

                self.array.values[self.array_inputs["range"], :, index] = (
                    (
                        (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_primary
                            * lhv_primary
                        )
                        + (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_secondary
                            * lhv_secondary
                        )
                        + (
                            self.array.values[self.array_inputs["fuel mass"], :, index]
                            * share_tertiary
                            * lhv_tertiary
                        )
                        + (
                            self.array.values[
                                self.array_inputs["electric energy stored"], :, index
                            ]
                            * 3.6
                        )
                    )
                    * 1000
                    / self.array.values[self.array_inputs["TtW energy"], :, index]
                )

                self.array.values[self.array_inputs["LHV fuel MJ per kg"], :, index] = (
                    (share_primary * lhv_primary)
                    + (share_secondary * lhv_secondary)
                    + (share_tertiary * lhv_tertiary)
                )

    def define_fuel_blends(self):
        """
        This function defines fuel blends from what is passed in `background_configuration`.
        It populates a dictionary `self.fuel_blends` that contains the respective shares, lower heating values
        and CO2 emission factors of the fuels used.

        Source for LHV: https://www.bafu.admin.ch/bafu/en/home/topics/climate/state/data/climate-reporting/references.html

        :return:
        """

        if {"ICEV-p", "HEV-p", "PHEV-p"}.intersection(set(self.scope["powertrain"])):
            fuel_type = "petrol"
            (
                primary,
                secondary,
                primary_share,
                secondary_share,
                tertiary,
                tertiary_share,
            ) = self.find_fuel_shares(fuel_type)
            self.create_fuel_markets(
                fuel_type,
                primary,
                secondary,
                tertiary,
                primary_share,
                secondary_share,
                tertiary_share,
            )
            self.fuel_blends[fuel_type] = {
                "primary": {
                    "type": primary,
                    "share": primary_share,
                    "lhv": self.fuel_specs[primary]["lhv"],
                    "CO2": self.fuel_specs[primary]["co2"],
                },
                "secondary": {
                    "type": secondary,
                    "share": secondary_share,
                    "lhv": self.fuel_specs[secondary]["lhv"],
                    "CO2": self.fuel_specs[secondary]["co2"],
                },
            }

            if tertiary:
                self.fuel_blends[fuel_type]["tertiary"] = {
                    "type": tertiary,
                    "share": tertiary_share,
                    "lhv": self.fuel_specs[tertiary]["lhv"],
                    "CO2": self.fuel_specs[tertiary]["co2"],
                }

        if {"ICEV-d", "HEV-d", "PHEV-d"}.intersection(set(self.scope["powertrain"])):
            fuel_type = "diesel"
            (
                primary,
                secondary,
                primary_share,
                secondary_share,
                tertiary,
                tertiary_share,
            ) = self.find_fuel_shares(fuel_type)
            self.create_fuel_markets(
                fuel_type,
                primary,
                secondary,
                tertiary,
                primary_share,
                secondary_share,
                tertiary_share,
            )
            self.fuel_blends[fuel_type] = {
                "primary": {
                    "type": primary,
                    "share": primary_share,
                    "lhv": self.fuel_specs[primary]["lhv"],
                    "CO2": self.fuel_specs[primary]["co2"],
                },
                "secondary": {
                    "type": secondary,
                    "share": secondary_share,
                    "lhv": self.fuel_specs[secondary]["lhv"],
                    "CO2": self.fuel_specs[secondary]["co2"],
                },
            }

            if tertiary:
                self.fuel_blends[fuel_type]["tertiary"] = {
                    "type": tertiary,
                    "share": tertiary_share,
                    "lhv": self.fuel_specs[tertiary]["lhv"],
                    "CO2": self.fuel_specs[tertiary]["co2"],
                }

        if {"ICEV-g"}.intersection(set(self.scope["powertrain"])):
            fuel_type = "cng"
            (
                primary,
                secondary,
                primary_share,
                secondary_share,
                tertiary,
                tertiary_share,
            ) = self.find_fuel_shares(fuel_type)
            self.create_fuel_markets(
                fuel_type,
                primary,
                secondary,
                tertiary,
                primary_share,
                secondary_share,
                tertiary_share,
            )
            self.fuel_blends[fuel_type] = {
                "primary": {
                    "type": primary,
                    "share": primary_share,
                    "lhv": self.fuel_specs[primary]["lhv"],
                    "CO2": self.fuel_specs[primary]["co2"],
                },
                "secondary": {
                    "type": secondary,
                    "share": secondary_share,
                    "lhv": self.fuel_specs[primary]["lhv"],
                    "CO2": self.fuel_specs[primary]["co2"],
                },
            }

            if tertiary:
                self.fuel_blends[fuel_type]["tertiary"] = {
                    "type": tertiary,
                    "share": tertiary_share,
                    "lhv": self.fuel_specs[tertiary]["lhv"],
                    "CO2": self.fuel_specs[tertiary]["co2"],
                }

        if {"FCEV"}.intersection(set(self.scope["powertrain"])):
            fuel_type = "hydrogen"
            (
                primary,
                secondary,
                primary_share,
                secondary_share,
                tertiary,
                tertiary_share,
            ) = self.find_fuel_shares(fuel_type)
            self.create_fuel_markets(
                fuel_type,
                primary,
                secondary,
                tertiary,
                primary_share,
                secondary_share,
                tertiary_share,
            )
            self.fuel_blends[fuel_type] = {
                "primary": {"type": primary, "share": primary_share},
                "secondary": {"type": secondary, "share": secondary_share},
            }

            if tertiary:
                self.fuel_blends[fuel_type]["tertiary"] = {
                    "type": tertiary,
                    "share": tertiary_share,
                }

        if {"BEV", "PHEV-p", "PHEV-d"}.intersection(set(self.scope["powertrain"])):
            fuel_type = "electricity"
            self.create_fuel_markets(fuel_type)

    def get_sulfur_content(self, location, fuel, year):
        """
        Return the sulfur content in the fuel.
        If a region is passed, the average sulfur content over
        the countries the region contains is returned.
        :param year:
        :param location: str. A country or region ISO code
        :param fuel: str. "diesel" or "gasoline
        :return: float. Sulfur content in ppm.
        """

        try:
            int(year)
        except ValueError:
            raise ValueError(
                "The year for which to fetch sulfur concentration values is not valid."
            )

        if location in self.bs.sulfur.country.values:
            sulfur_concentration = (
                self.bs.sulfur.sel(country=location, year=year, fuel=fuel).sum().values
            )
        else:
            # If the geography is not found,
            # we use the European average

            print(
                f"The sulfur content for {fuel} fuel in {location} could not be found."
                "European average sulfur content is used instead."
            )

            sulfur_concentration = (
                self.bs.sulfur.sel(country="RER", year=year, fuel=fuel).sum().values
            )

        return sulfur_concentration

    def create_fuel_dictionary(self):

        for val in self.fuel_specs.values():
            if any(
                i in val["name"][0].lower()
                for i in ("synthetic", "hydrogen", "ethanol", "biodiesel")
            ):
                val["additional electricity"] = self.find_inputs(
                    "kilowatt hour", val["name"][0], "unit"
                )
            else:
                val["additional electricity"] = 0

        for val in self.fuel_specs.values():
            if any(
                i in val["name"][0].lower() for i in ("synthetic", "hydrogen", "bio")
            ):
                self.find_inputs(
                    "kilowatt hour", val["name"][0], "unit", zero_out_input=True
                )

    def create_fuel_markets(
        self,
        fuel_type,
        primary=None,
        secondary=None,
        tertiary=None,
        primary_share=None,
        secondary_share=None,
        tertiary_share=None,
    ):
        """
        This function creates markets for fuel, considering a given blend, a given fuel type and a given year.
        It also adds separate electricity input in case hydrogen from electrolysis is needed somewhere in the fuel supply chain.
        :return:
        """

        d_dataset_name = {
            "petrol": "fuel supply for gasoline vehicles, ",
            "diesel": "fuel supply for diesel vehicles, ",
            "cng": "fuel supply for gas vehicles, ",
            "hydrogen": "fuel supply for hydrogen vehicles, ",
            "electricity": "electricity supply for electric vehicles, ",
        }

        if fuel_type != "electricity":
            for y, year in enumerate(self.scope["year"]):
                dataset_name = d_dataset_name[fuel_type] + str(year)
                fuel_market_index = [
                    self.inputs[i] for i in self.inputs if i[0] == dataset_name
                ][0]

                try:
                    primary_fuel_activity_index = self.inputs[
                        tuple(self.fuel_specs[primary]["name"])
                    ]
                    secondary_fuel_activity_index = self.inputs[
                        tuple(self.fuel_specs[secondary]["name"])
                    ]
                except KeyError:
                    raise KeyError(
                        "One of the primary or secondary fuels specified in "
                        "the fuel blend for {} is not valid.".format(fuel_type)
                    )

                if tertiary:

                    if ~np.isclose(
                        primary_share[y] + secondary_share[y] + tertiary_share[y],
                        1,
                        rtol=1e-3,
                    ):
                        sum_blend = (
                            primary_share[y] + secondary_share[y] + tertiary_share[y]
                        )
                        print(
                            f"The fuel blend for {fuel_type} in {year} is not equal to 1, but {sum_blend}."
                            f"The primary fuel share is adjusted so that the fuel blend equals 1."
                        )
                        primary_share[y] = 1 - (secondary_share[y] + tertiary_share[y])
                else:

                    if ~np.isclose(primary_share[y] + secondary_share[y], 1, rtol=1e-3):
                        sum_blend = primary_share[y] + secondary_share[y]
                        print(
                            f"The fuel blend for {fuel_type} in {year} is not equal to 1, but {sum_blend}."
                            f"The primary fuel share is adjusted so that the fuel blend equals 1."
                        )
                        primary_share[y] = 1 - secondary_share[y]

                self.A[:, primary_fuel_activity_index, fuel_market_index] = (
                    -1 * primary_share[y]
                )
                self.A[:, secondary_fuel_activity_index, fuel_market_index] = (
                    -1 * secondary_share[y]
                )

                def learning_rate_fuel(fuel, year, share, val):
                    if fuel == "electrolysis":
                        # apply some learning rate for electrolysis
                        electrolysis = -0.3538 * (float(year) - 2010) + 58.589
                        electricity = (val - 58 + electrolysis) * share
                    elif fuel == "synthetic gasoline - energy allocation":
                        # apply some learning rate for electrolysis
                        h2 = 0.338
                        electrolysis = -0.3538 * (float(year) - 2010) + 58.589
                        electricity = val - (h2 * 58)
                        electricity += electrolysis * h2
                        electricity *= share

                    elif fuel == "synthetic gasoline - economic allocation":
                        # apply some learning rate for electrolysis
                        h2 = 0.6385
                        electrolysis = -0.3538 * (float(year) - 2010) + 58.589
                        electricity = val - (h2 * 58)
                        electricity += electrolysis * h2
                        electricity *= share

                    elif fuel == "synthetic diesel - energy allocation":
                        # apply some learning rate for electrolysis
                        h2 = 0.42
                        electrolysis = -0.3538 * (float(year) - 2010) + 58.589
                        electricity = val - (h2 * 58)
                        electricity += electrolysis * h2
                        electricity *= share

                    elif fuel == "synthetic diesel - economic allocation":
                        # apply some learning rate for electrolysis
                        h2 = 0.183
                        electrolysis = -0.3538 * (float(year) - 2010) + 58.589
                        electricity = val - (h2 * 58)
                        electricity += electrolysis * h2
                        electricity *= share

                    else:
                        electricity = val * share
                    return electricity

                additional_electricity_primary = learning_rate_fuel(
                    primary,
                    year,
                    primary_share[y],
                    self.fuel_specs[primary]["additional electricity"],
                )

                additional_electricity_secondary = learning_rate_fuel(
                    secondary,
                    year,
                    secondary_share[y],
                    self.fuel_specs[secondary]["additional electricity"],
                )

                additional_electricity = (
                    additional_electricity_primary + additional_electricity_secondary
                )

                if tertiary:
                    tertiary_fuel_activity_index = self.inputs[
                        self.fuel_specs[tertiary]["name"]
                    ]
                    self.A[:, tertiary_fuel_activity_index, fuel_market_index] = (
                        -1 * tertiary_share[y]
                    )
                    additional_electricity += learning_rate_fuel(
                        tertiary,
                        year,
                        tertiary_share[y],
                        self.fuel_specs[tertiary]["additional electricity"],
                    )

                if additional_electricity > 0:
                    electricity_mix_index = [
                        self.inputs[i]
                        for i in self.inputs
                        if i[0]
                        == "electricity market for fuel preparation, " + str(year)
                    ][0]
                    self.A[:, electricity_mix_index, fuel_market_index] = (
                        -1 * additional_electricity
                    )
        else:
            for year in self.scope["year"]:
                dataset_name = d_dataset_name[fuel_type] + str(year)
                electricity_market_index = [
                    self.inputs[i] for i in self.inputs if i[0] == dataset_name
                ][0]
                electricity_mix_index = [
                    self.inputs[i]
                    for i in self.inputs
                    if i[0] == "electricity market for fuel preparation, " + str(year)
                ][0]
                self.A[:, electricity_mix_index, electricity_market_index] = -1

    def find_inputs(
        self, value_in, value_out, find_input_by="name", zero_out_input=False
    ):
        """
        Finds the exchange inputs to a specified functional unit
        :param zero_out_input:
        :param find_input_by: can be 'name' or 'unit'
        :param value_in: value to look for
        :param value_out: functional unit output
        :return: indices of all inputs to FU, indices of inputs of interest
        :rtype: tuple
        """

        if isinstance(value_out, str):
            value_out = [value_out]

        index_output = [
            self.inputs[i]
            for val in value_out
            for i in self.inputs
            if val.lower() in i[0].lower()
        ]

        f_vector = np.zeros((np.shape(self.A)[1]))

        f_vector[index_output] = 1

        X = np.float32(sparse.linalg.spsolve(sparse.csr_matrix(self.A[0]), f_vector.T))

        ind_inputs = np.nonzero(X)[0]

        if find_input_by == "name":
            ins = [
                i
                for i in ind_inputs
                if value_in.lower() in self.rev_inputs[i][0].lower()
            ]

        if find_input_by == "unit":
            ins = [
                i
                for i in ind_inputs
                if value_in.lower() in self.rev_inputs[i][2].lower()
            ]

        outs = [i for i in ind_inputs if i not in ins]

        sum_supplied = X[ins].sum()

        if zero_out_input:
            # zero out initial inputs
            self.A[np.ix_(np.arange(0, self.A.shape[0]), ins, outs)] *= 0

        else:
            return sum_supplied

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
                ("Glider lightweighting", "GLO", "kilogram", "glider lightweighting")
            ],
            -self.number_of_cars :,
        ] = (
            (
                array[self.array_inputs["lightweighting"], :]
                * array[self.array_inputs["glider base mass"], :]
            )
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
        ] = (
            array[self.array_inputs["curb mass"], :] / 1240 / 150000 * -1
        )

        # Glider EoL + fuel tank
        self.A[
            :,
            self.inputs[
                (
                    "treatment of used glider, passenger car, shredding",
                    "GLO",
                    "kilogram",
                    "used glider, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            (
                array[self.array_inputs["glider base mass"], :]
                * (1 - array[self.array_inputs["lightweighting"], :])
            )
            + array[self.array_inputs["fuel tank mass"], :]
        ) / array[
            self.array_inputs["lifetime kilometers"], :
        ]

        # Battery EoL
        self.A[
            :,
            self.inputs[
                (
                    "market for used Li-ion battery",
                    "GLO",
                    "kilogram",
                    "used Li-ion battery",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[
                [
                    self.array_inputs[l]
                    for l in ["battery cell mass", "battery BoP mass"]
                ],
                :,
            ].sum(axis=0)
            / array[self.array_inputs["lifetime kilometers"], :]
        )

        # Combustion engine + powertrain EoL
        self.A[
            :,
            self.inputs[
                (
                    "treatment of used internal combustion engine, passenger car, shredding",
                    "GLO",
                    "kilogram",
                    "used internal combustion engine, passenger car",
                )
            ],
            -self.number_of_cars :,
        ] = (
            array[
                [
                    self.array_inputs[l]
                    for l in ["combustion engine mass", "powertrain mass"]
                ],
                :,
            ].sum(axis=0)
            / array[self.array_inputs["lifetime kilometers"], :]
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

        # Electric engine + powertrain EoL

        l_elec_pt = [
            "charger mass",
            "converter mass",
            "inverter mass",
            "power distribution unit mass",
            # "powertrain mass",
            "electric engine mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
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
                    "internal combustion engine, passenger car",
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
            * (1 + array[self.array_inputs["fuel cell lifetime replacements"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[("Essential BoP", "GLO", "kilogram", "Essential BoP")],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["fuel cell essential BoP mass"], :]
            * (1 + array[self.array_inputs["fuel cell lifetime replacements"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        self.A[
            :,
            self.inputs[("Stack", "GLO", "kilowatt", "Stack")],
            -self.number_of_cars :,
        ] = (
            array[self.array_inputs["fuel cell stack mass"], :]
            / 0.51
            * (1 + array[self.array_inputs["fuel cell lifetime replacements"], :])
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Start of printout

        print(
            "****************** IMPORTANT BACKGROUND PARAMETERS ******************",
            end="\n * ",
        )

        # Energy storage

        print(
            f"The country of use is {self.country}",
            end="\n * ",
        )

        battery_tech = self.background_configuration["energy storage"]["electric"][
            "type"
        ]
        battery_origin = self.background_configuration["energy storage"]["electric"][
            "origin"
        ]

        print(
            f"Power and energy batteries produced in {battery_origin} using {battery_tech} chemistry",
            end="\n * ",
        )

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
            f"Battery cell, {battery_tech}",
            "GLO",
            "kilogram",
            "Battery cell",
        )

        self.A[:, self.inputs[battery_cell_label], -self.number_of_cars :,] = (
            (
                array[self.array_inputs["battery cell mass"], :]
                * (1 + array[self.array_inputs["battery lifetime replacements"], :])
            )
            / array[self.array_inputs["lifetime kilometers"], :]
            * -1
        )

        # Fetch the overall input of electricity per kg of battery cell
        electricity_batt = self.find_inputs(
            "kilowatt hour", f"Battery cell, {battery_tech}", "unit"
        )

        for y in self.scope["year"]:
            index = self.get_index_vehicle_from_array(y)

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
                        if str(y) in i[0] and "transport, passenger" in i[0]
                    ],
                )
            ] = (
                (
                    electricity_batt
                    * array[self.array_inputs["battery cell mass"], :, index]
                    * (
                        1
                        + array[
                            self.array_inputs["battery lifetime replacements"],
                            :,
                            index,
                        ]
                    )
                )
                / array[self.array_inputs["lifetime kilometers"], :, index]
                * -1
            ).T[
                :, None, :
            ]

        self.find_inputs(
            "kilowatt hour",
            f"Battery cell, {battery_tech}",
            "unit",
            zero_out_input=True,
        )

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if any(
                ele in c[0]
                for ele in ["ICEV-d", "ICEV-p", "HEV-p", "PHEV-p", "PHEV-d", "HEV-d"]
            )
        ]
        index = self.get_index_vehicle_from_array(
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

        index = self.get_index_vehicle_from_array("ICEV-g")
        self.A[
            :,
            self.inputs[
                (
                    "Fuel tank, compressed natural gas, 200 bar",
                    "RER",
                    "kilogram",
                    "Fuel tank, compressed natural gas, 200 bar",
                )
            ],
            [
                self.inputs[i]
                for i in self.inputs
                if i[0].startswith("transport, passenger") and "ICEV-g" in i[0]
            ],
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        if "hydrogen" in self.background_configuration["energy storage"]:
            # If a customization dict is passed
            hydro_tank_technology = self.background_configuration["energy storage"][
                "hydrogen"
            ]["type"]
        else:
            hydro_tank_technology = "carbon fiber"

        dict_tank_map = {
            "carbon fiber": (
                "Fuel tank, compressed hydrogen gas, 700bar",
                "GLO",
                "kilogram",
                "Fuel tank, compressed hydrogen gas, 700bar",
            ),
            "hdpe": (
                "Fuel tank, compressed hydrogen gas, 700bar, with HDPE liner",
                "RER",
                "kilogram",
                "Hydrogen tank",
            ),
            "aluminium": (
                "Fuel tank, compressed hydrogen gas, 700bar, with aluminium liner",
                "RER",
                "kilogram",
                "Hydrogen tank",
            ),
        }

        index = self.get_index_vehicle_from_array("FCEV")
        self.A[
            :,
            self.inputs[dict_tank_map[hydro_tank_technology]],
            [
                self.inputs[i]
                for i in self.inputs
                if i[0].startswith("transport, passenger") and "FCEV" in i[0]
            ],
        ] = (
            array[self.array_inputs["fuel tank mass"], :, index]
            / array[self.array_inputs["lifetime kilometers"], :, index]
            * -1
        ).T

        try:
            sum_renew, co2_intensity_tech = self.define_renewable_rate_in_mix()

        except AttributeError:
            sum_renew = [0] * len(self.scope["year"])
            co2_intensity_tech = [0] * len(self.scope["year"])

        for y, year in enumerate(self.scope["year"]):

            if y + 1 == len(self.scope["year"]):
                end_str = "\n * "
            else:
                end_str = "\n \t * "

            print(
                "in "
                + str(year)
                + ", % of renewable: "
                + str(np.round(sum_renew[y] * 100, 0))
                + "%"
                + ", GHG intensity per kWh: "
                + str(int(np.sum(co2_intensity_tech[y] * self.mix[y])))
                + " g. CO2-eq.",
                end=end_str,
            )

        if any(
            True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in self.scope["powertrain"]
        ):
            for y in self.scope["year"]:
                index = self.get_index_vehicle_from_array(
                    ["BEV", "PHEV-p", "PHEV-d"], y, method="and"
                )

                self.A[
                    np.ix_(
                        np.arange(self.iterations),
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "electricity supply for electric vehicles" in i[0]
                        ],
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "transport, passenger" in i[0]
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

            index = self.get_index_vehicle_from_array("FCEV")

            if "tertiary" in self.fuel_blends["hydrogen"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["hydrogen"]["primary"]["type"],
                        self.fuel_blends["hydrogen"]["secondary"]["type"],
                        self.fuel_blends["hydrogen"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    f"{self.fuel_blends['hydrogen']['primary']['type']} "
                    f"is completed by {self.fuel_blends['hydrogen']['secondary']['type']}.",
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["hydrogen"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["tertiary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                # Fuel supply

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and "FCEV" in i[0]
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "fuel supply for hydrogen vehicles" in i[0]
                    ],
                    ind_A,
                ] = (
                    array[self.array_inputs["fuel mass"], :, ind_array]
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if "ICEV-g" in self.scope["powertrain"]:
            index = self.get_index_vehicle_from_array("ICEV-g")

            if "tertiary" in self.fuel_blends["cng"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["cng"]["primary"]["type"],
                        self.fuel_blends["cng"]["secondary"]["type"],
                        self.fuel_blends["cng"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    "{} is completed by {}.".format(
                        self.fuel_blends["cng"]["primary"]["type"],
                        self.fuel_blends["cng"]["secondary"]["type"],
                    ),
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["cng"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["secondary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["tertiary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["secondary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                # Primary fuel share

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and "ICEV-g" in i[0]
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Supply of gas, including 0.4% of gas input as pump-to-tank leakage
                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0] and "fuel supply for gas vehicles" in i[0]
                    ],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        / array[self.array_inputs["range"], :, ind_array]
                    )
                    * (
                        1
                        + array[
                            self.array_inputs["CNG pump-to-tank leakage"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                # Gas leakage to air

                self.A[
                    :,
                    self.inputs[("Methane, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        / array[self.array_inputs["range"], :, ind_array]
                    )
                    * array[self.array_inputs["CNG pump-to-tank leakage"], :, ind_array]
                    * -1
                ).T

                # Fuel-based emissions from CNG, CO2
                # The share and CO2 emissions factor of CNG is retrieved, if used

                share_fossil = 0
                CO2_fossil = 0

                if self.fuel_blends["cng"]["primary"]["type"] == "cng":
                    share_fossil = self.fuel_blends["cng"]["primary"]["share"][y]
                    CO2_fossil = (
                        self.fuel_blends["cng"]["primary"]["CO2"]
                        * self.fuel_blends["cng"]["primary"]["share"][y]
                    )

                if self.fuel_blends["cng"]["secondary"]["type"] == "cng":
                    share_fossil += self.fuel_blends["cng"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["cng"]["secondary"]["CO2"]
                        * self.fuel_blends["cng"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["cng"]:
                    if self.fuel_blends["cng"]["tertiary"]["type"] == "cng":
                        share_fossil += self.fuel_blends["cng"]["tertiary"]["share"][y]
                        CO2_fossil += (
                            self.fuel_blends["cng"]["tertiary"]["CO2"]
                            * self.fuel_blends["cng"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    array[self.array_inputs["fuel mass"], :, ind_array]
                    * CO2_fossil
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # Fuel-based CO2 emission from alternative CNG
                # The share of non-fossil gas in the blend is retrieved
                # As well as the CO2 emission factor of the fuel

                share_non_fossil = 0
                CO2_non_fossil = 0

                if self.fuel_blends["cng"]["primary"]["type"] != "cng":
                    share_non_fossil = self.fuel_blends["cng"]["primary"]["share"][y]
                    CO2_non_fossil = (
                        self.fuel_blends["cng"]["primary"]["share"][y]
                        * self.fuel_blends["cng"]["primary"]["CO2"]
                    )

                if self.fuel_blends["cng"]["secondary"]["type"] != "cng":
                    share_non_fossil += self.fuel_blends["cng"]["secondary"]["share"][y]
                    CO2_non_fossil += (
                        self.fuel_blends["cng"]["secondary"]["share"][y]
                        * self.fuel_blends["cng"]["secondary"]["CO2"]
                    )

                if "tertiary" in self.fuel_blends["cng"]:
                    if self.fuel_blends["cng"]["tertiary"]["type"] != "cng":
                        share_non_fossil += self.fuel_blends["cng"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["cng"]["tertiary"]["share"][y]
                            * self.fuel_blends["cng"]["tertiary"]["CO2"]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    (
                        (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * CO2_non_fossil
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if [i for i in self.scope["powertrain"] if i in ["ICEV-d", "PHEV-d", "HEV-d"]]:
            index = self.get_index_vehicle_from_array(["ICEV-d", "PHEV-d", "HEV-d"])

            if "tertiary" in self.fuel_blends["diesel"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["diesel"]["primary"]["type"],
                        self.fuel_blends["diesel"]["secondary"]["type"],
                        self.fuel_blends["diesel"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    "{} is completed by {}.".format(
                        self.fuel_blends["diesel"]["primary"]["type"],
                        self.fuel_blends["diesel"]["secondary"]["type"],
                    ),
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["diesel"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["diesel"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["diesel"]["tertiary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["diesel"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and any(x in i[0] for x in ["ICEV-d", "PHEV-d", "HEV-d"])
                ]

                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Fuel supply

                fuel_amount = (
                    array[
                        self.array_inputs["TtW energy, combustion mode"], :, ind_array
                    ]
                    / (
                        self.array.values[
                            self.array_inputs["LHV fuel MJ per kg"], :, ind_array
                        ]
                        * 1000
                    )
                    * (
                        1
                        - self.array.values[
                            self.array_inputs["electric utility factor"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "fuel supply for diesel vehicles" in i[0]
                    ],
                    ind_A,
                ] = fuel_amount

                share_fossil = 0
                CO2_fossil = 0

                # Fuel-based CO2 emission from conventional diesel
                if self.fuel_blends["diesel"]["primary"]["type"] == "diesel":
                    share_fossil += self.fuel_blends["diesel"]["primary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["diesel"]["primary"]["CO2"]
                        * self.fuel_blends["diesel"]["primary"]["share"][y]
                    )

                if self.fuel_blends["diesel"]["secondary"]["type"] == "diesel":
                    share_fossil += self.fuel_blends["diesel"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["diesel"]["secondary"]["CO2"]
                        * self.fuel_blends["diesel"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["diesel"]:
                    if self.fuel_blends["diesel"]["tertiary"]["type"] == "diesel":
                        share_fossil += self.fuel_blends["diesel"]["tertiary"]["share"][
                            y
                        ]
                        CO2_fossil += (
                            self.fuel_blends["diesel"]["tertiary"]["CO2"]
                            * self.fuel_blends["diesel"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount * CO2_fossil
                )

                # Fuel-based SO2 emissions
                # Sulfur concentration value for a given country, a given year, as concentration ratio

                sulfur_concentration = self.get_sulfur_content(
                    self.country, "diesel", year
                )

                self.A[
                    :,
                    self.inputs[("Sulfur dioxide", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount
                    * sulfur_concentration
                    * (64 / 32)  # molar mass of SO2/molar mass of O2
                )

                share_non_fossil = 0
                CO2_non_fossil = 0

                # Fuel-based CO2 emission from alternative petrol
                # The share of non-fossil fuel in the blend is retrieved
                # As well as the CO2 emission factor of the fuel
                if self.fuel_blends["diesel"]["primary"]["type"] != "diesel":
                    share_non_fossil += self.fuel_blends["diesel"]["primary"]["share"][
                        y
                    ]
                    CO2_non_fossil += (
                        self.fuel_blends["diesel"]["primary"]["CO2"]
                        * self.fuel_blends["diesel"]["primary"]["share"][y]
                    )

                if self.fuel_blends["diesel"]["secondary"]["type"] != "diesel":
                    share_non_fossil += self.fuel_blends["diesel"]["secondary"][
                        "share"
                    ][y]
                    CO2_non_fossil += (
                        self.fuel_blends["diesel"]["secondary"]["share"][y]
                        * self.fuel_blends["diesel"]["secondary"]["CO2"]
                    )

                if "tertiary" in self.fuel_blends["diesel"]:
                    if self.fuel_blends["diesel"]["tertiary"]["type"] != "diesel":
                        share_non_fossil += self.fuel_blends["diesel"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["diesel"]["tertiary"]["share"][y]
                            * self.fuel_blends["diesel"]["tertiary"]["CO2"]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    fuel_amount * CO2_non_fossil
                )

        if [i for i in self.scope["powertrain"] if i in ["ICEV-p", "HEV-p", "PHEV-p"]]:
            index = self.get_index_vehicle_from_array(["ICEV-p", "HEV-p", "PHEV-p"])

            if "tertiary" in self.fuel_blends["petrol"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["petrol"]["primary"]["type"],
                        self.fuel_blends["petrol"]["secondary"]["type"],
                        self.fuel_blends["petrol"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    "{} is completed by {}.".format(
                        self.fuel_blends["petrol"]["primary"]["type"],
                        self.fuel_blends["petrol"]["secondary"]["type"],
                    ),
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["petrol"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["petrol"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["petrol"]["tertiary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["petrol"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and any(x in i[0] for x in ["ICEV-p", "HEV-p", "PHEV-p"])
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Fuel supply

                fuel_amount = (
                    array[
                        self.array_inputs["TtW energy, combustion mode"], :, ind_array
                    ]
                    / (
                        self.array.values[
                            self.array_inputs["LHV fuel MJ per kg"], :, ind_array
                        ]
                        * 1000
                    )
                    * (
                        1
                        - self.array.values[
                            self.array_inputs["electric utility factor"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "fuel supply for gasoline vehicles" in i[0]
                    ],
                    ind_A,
                ] = fuel_amount

                share_fossil = 0
                CO2_fossil = 0

                # Fuel-based CO2 emission from conventional petrol
                if self.fuel_blends["petrol"]["primary"]["type"] == "petrol":
                    share_fossil = self.fuel_blends["petrol"]["primary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["petrol"]["primary"]["CO2"]
                        * self.fuel_blends["petrol"]["primary"]["share"][y]
                    )

                if self.fuel_blends["petrol"]["secondary"]["type"] == "petrol":
                    share_fossil += self.fuel_blends["petrol"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["petrol"]["secondary"]["CO2"]
                        * self.fuel_blends["petrol"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["petrol"]:
                    if self.fuel_blends["petrol"]["tertiary"]["type"] == "petrol":
                        share_fossil += self.fuel_blends["petrol"]["tertiary"]["share"][
                            y
                        ]
                        CO2_fossil += (
                            self.fuel_blends["petrol"]["tertiary"]["CO2"]
                            * self.fuel_blends["petrol"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount * CO2_fossil
                )

                # Fuel-based SO2 emissions
                # Sulfur concentration value for a given country, a given year, as a concentration ratio

                sulfur_concentration = self.get_sulfur_content(
                    self.country, "petrol", year
                )

                self.A[
                    :,
                    self.inputs[("Sulfur dioxide", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount
                    * sulfur_concentration
                    * (64 / 32)  # molar mass of SO2/molar mass of O2
                )

                share_non_fossil = 0
                CO2_non_fossil = 0

                # Fuel-based CO2 emission from alternative petrol
                # The share of non-fossil fuel in the blend is retrieved
                # As well as the CO2 emission factor of the fuel
                if self.fuel_blends["petrol"]["primary"]["type"] != "petrol":
                    share_non_fossil += self.fuel_blends["petrol"]["primary"]["share"][
                        y
                    ]
                    CO2_non_fossil += (
                        self.fuel_blends["petrol"]["primary"]["CO2"]
                        * self.fuel_blends["petrol"]["primary"]["share"][y]
                    )

                if self.fuel_blends["petrol"]["secondary"]["type"] != "petrol":
                    share_non_fossil += self.fuel_blends["petrol"]["secondary"][
                        "share"
                    ][y]
                    CO2_non_fossil += (
                        self.fuel_blends["petrol"]["secondary"]["CO2"]
                        * self.fuel_blends["petrol"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["petrol"]:
                    if self.fuel_blends["petrol"]["tertiary"]["type"] != "petrol":
                        share_non_fossil += self.fuel_blends["petrol"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["petrol"]["tertiary"]["share"][y]
                            * self.fuel_blends["petrol"]["tertiary"]["CO2"]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    fuel_amount * CO2_non_fossil
                )

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
        ] = array[self.array_inputs["road wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

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
        ] = array[self.array_inputs["tire wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

        # Brake wear emissions
        ind_A = [self.inputs[i] for i in self.inputs if "transport, passenger" in i[0]]

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
            ind_A,
        ] = array[self.array_inputs["brake wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

        # Infrastructure
        self.A[
            :,
            self.inputs[("market for road", "GLO", "meter-year", "road")],
            -self.number_of_cars :,
        ] = (
            5.37e-7 * array[self.array_inputs["driving mass"], :] * -1
        )

        # Infrastructure maintenance
        self.A[
            :,
            self.inputs[
                ("market for road maintenance", "RER", "meter-year", "road maintenance")
            ],
            -self.number_of_cars :,
        ] = (
            1.29e-3 * -1
        )

        # Exhaust emissions
        # Non-fuel based emissions
        self.A[:, self.index_emissions, -self.number_of_cars :] = (
            array[
                [
                    self.array_inputs[self.map_fuel_emissions[self.rev_inputs[x]]]
                    for x in self.index_emissions
                ]
            ]
            * -1
        ).transpose([1, 0, 2])

        # Noise emissions
        self.A[:, self.index_noise, -self.number_of_cars :] = (
            array[
                [
                    self.array_inputs[self.map_noise_emissions[self.rev_inputs[x]]]
                    for x in self.index_noise
                ]
            ]
            * -1
        ).transpose([1, 0, 2])

        # Emissions of air conditioner refrigerant r134a
        # Leakage assumed to amount to 750g/lifetime according to
        # https://treeze.ch/fileadmin/user_upload/downloads/Publications/Case_Studies/Mobility/544-LCI-Road-NonRoad-Transport-Services-v2.0.pdf
        # but only to cars with an AC system (meaning, with a cooling energy consumption)
        self.A[
            :,
            self.inputs[
                ("Ethane, 1,1,1,2-tetrafluoro-, HFC-134a", ("air",), "kilogram")
            ],
            -self.number_of_cars :,
        ] = (
            0.750 / self.array.values[self.array_inputs["lifetime kilometers"]] * -1
        ) * self.array.values[
            self.array_inputs["cooling energy consumption"]
        ]

        self.A[
            :,
            self.inputs[
                ("market for refrigerant R134a", "GLO", "kilogram", "refrigerant R134a")
            ],
            -self.number_of_cars :,
        ] = (
            (0.75 + 0.55)
            / self.array.values[self.array_inputs["lifetime kilometers"]]
            * -1
        ) * self.array.values[
            self.array_inputs["cooling energy consumption"]
        ]

        print("*********************************************************************")

    def set_inputs_in_A_matrix_for_export(self, array):
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (array[self.array_inputs["glider base mass"], :]) * -1

        self.A[
            :,
            self.inputs[
                ("Glider lightweighting", "GLO", "kilogram", "glider lightweighting")
            ],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["lightweighting"], :]
            * array[self.array_inputs["glider base mass"], :]
        ) * -1

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
            [self.inputs[i] for i in self.inputs if "transport, passenger car" in i[0]],
        ] = (
            array[self.array_inputs["curb mass"], :] / 1240 / 150000 * -1
        )

        # Glider EoL + fuel tank
        self.A[
            :,
            self.inputs[
                (
                    "treatment of used glider, passenger car, shredding",
                    "GLO",
                    "kilogram",
                    "used glider, passenger car",
                )
            ],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["glider base mass"], :]
            * (1 - array[self.array_inputs["lightweighting"], :])
        ) + array[
            self.array_inputs["fuel tank mass"], :
        ]

        # Battery EoL
        self.A[
            :,
            self.inputs[
                (
                    "market for used Li-ion battery",
                    "GLO",
                    "kilogram",
                    "used Li-ion battery",
                )
            ],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = array[
            [self.array_inputs[l] for l in ["battery cell mass", "battery BoP mass"]], :
        ].sum(
            axis=0
        )

        # Combustion engine EoL
        self.A[
            :,
            self.inputs[
                (
                    "treatment of used internal combustion engine, passenger car, shredding",
                    "GLO",
                    "kilogram",
                    "used internal combustion engine, passenger car",
                )
            ],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = array[
            [
                self.array_inputs[l]
                for l in ["combustion engine mass", "powertrain mass"]
            ],
            :,
        ].sum(
            axis=0
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["charger mass"], :] * -1
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["converter mass"], :] * -1
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["electric engine mass"], :] * -1
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["inverter mass"], :] * -1
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["power distribution unit mass"], :] * -1
        )

        l_elec_pt = [
            "charger mass",
            "converter mass",
            "inverter mass",
            "power distribution unit mass",
            # "powertrain mass",
            "electric engine mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
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
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)

        self.A[
            :,
            self.inputs[
                (
                    "market for internal combustion engine, passenger car",
                    "GLO",
                    "kilogram",
                    "internal combustion engine, passenger car",
                )
            ],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[
                [
                    self.array_inputs[l]
                    for l in ["combustion engine mass", "powertrain mass"]
                ],
                :,
            ].sum(axis=0)
        ) * -1

        self.A[
            :,
            self.inputs[("Ancillary BoP", "GLO", "kilogram", "Ancillary BoP")],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["fuel cell ancillary BoP mass"], :] * -1
        )

        self.A[
            :,
            self.inputs[("Essential BoP", "GLO", "kilogram", "Essential BoP")],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["fuel cell essential BoP mass"], :] * -1
        )

        self.A[
            :,
            self.inputs[("Stack", "GLO", "kilowatt", "Stack")],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["fuel cell stack mass"], :] * -1
        )

        # Start of printout

        print(
            "****************** IMPORTANT BACKGROUND PARAMETERS ******************",
            end="\n * ",
        )

        # Energy storage

        print(
            "The country of use is " + self.country,
            end="\n * ",
        )

        battery_tech = self.background_configuration["energy storage"]["electric"][
            "type"
        ]
        battery_origin = self.background_configuration["energy storage"]["electric"][
            "origin"
        ]

        print(
            "Power and energy batteries produced in "
            + battery_origin
            + " using "
            + battery_tech
            + " chemistry.",
            end="\n * ",
        )

        # Use the NMC inventory of Schmidt et al. 2019
        self.A[
            :,
            self.inputs[("Battery BoP", "GLO", "kilogram", "Battery BoP")],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["battery BoP mass"], :]
            * (1 + array[self.array_inputs["battery lifetime replacements"], :])
        ) * -1

        battery_cell_label = (
            "Battery cell, " + battery_tech,
            "GLO",
            "kilogram",
            "Battery cell",
        )

        self.A[
            :,
            self.inputs[battery_cell_label],
            [self.inputs[i] for i in self.inputs if i[0].startswith("Passenger car")],
        ] = (
            array[self.array_inputs["battery cell mass"], :]
            * (1 + array[self.array_inputs["fuel cell lifetime replacements"], :])
        ) * -1

        # Set an input of electricity, given the country of manufacture
        electricity_batt = self.find_inputs(
            "kilowatt hour",
            f"Battery cell, {battery_tech}",
            "unit",
        )

        for y in self.scope["year"]:
            index = self.get_index_vehicle_from_array(y)

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
                        if str(y) in i[0] and i[0].startswith("Passenger car")
                    ],
                )
            ] = (
                electricity_batt
                * (
                    (
                        array[self.array_inputs["battery cell mass"], :, index]
                        * (
                            1
                            + array[
                                self.array_inputs["battery lifetime replacements"],
                                :,
                                index,
                            ]
                        )
                    )
                    * -1
                )
            ).T[
                :, None, :
            ]

        self.find_inputs(
            "kilowatt hour",
            f"Battery cell, {battery_tech}",
            "unit",
            zero_out_input=True,
        )

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if any(
                ele in c[0]
                for ele in ["ICEV-d", "ICEV-p", "HEV-p", "PHEV-p", "PHEV-d", "HEV-d"]
            )
            and c[0].startswith("Passenger car")
        ]
        index = self.get_index_vehicle_from_array(
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
        ] = (array[self.array_inputs["fuel tank mass"], :, index] * -1).T

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if "ICEV-g" in c[0] and c[0].startswith("Passenger car")
        ]

        index = self.get_index_vehicle_from_array("ICEV-g")
        self.A[
            :,
            self.inputs[
                (
                    "Fuel tank, compressed natural gas, 200 bar",
                    "RER",
                    "kilogram",
                    "Fuel tank, compressed natural gas, 200 bar",
                )
            ],
            index_A,
        ] = (array[self.array_inputs["fuel tank mass"], :, index] * -1).T

        if "hydrogen" in self.background_configuration["energy storage"]:
            # If a customization dict is passed
            hydro_tank_technology = self.background_configuration["energy storage"][
                "hydrogen"
            ]["type"]
        else:
            hydro_tank_technology = "carbon fiber"

        dict_tank_map = {
            "carbon fiber": (
                "Fuel tank, compressed hydrogen gas, 700bar",
                "GLO",
                "kilogram",
                "Fuel tank, compressed hydrogen gas, 700bar",
            ),
            "hdpe": (
                "Fuel tank, compressed hydrogen gas, 700bar, with HDPE liner",
                "RER",
                "kilogram",
                "Hydrogen tank",
            ),
            "aluminium": (
                "Fuel tank, compressed hydrogen gas, 700bar, with aluminium liner",
                "RER",
                "kilogram",
                "Hydrogen tank",
            ),
        }

        index_A = [
            self.inputs[c]
            for c in self.inputs
            if "FCEV" in c[0] and c[0].startswith("Passenger car")
        ]

        index = self.get_index_vehicle_from_array("FCEV")
        self.A[:, self.inputs[dict_tank_map[hydro_tank_technology]], index_A] = (
            array[self.array_inputs["fuel tank mass"], :, index] * -1
        ).T

        # END of vehicle building

        self.A[
            :,
            [self.inputs[c] for c in self.inputs if c[0].startswith("Passenger car")],
            [self.inputs[c] for c in self.inputs if "transport, passenger car" in c[0]],
        ] = (
            -1 / array[self.array_inputs["lifetime kilometers"]]
        )

        try:
            sum_renew, co2_intensity_tech = self.define_renewable_rate_in_mix()

        except AttributeError:
            sum_renew = [0] * len(self.scope["year"])
            co2_intensity_tech = [0] * len(self.scope["year"])

        for y, year in enumerate(self.scope["year"]):

            if y + 1 == len(self.scope["year"]):
                end_str = "\n * "
            else:
                end_str = "\n \t * "

            print(
                "in "
                + str(year)
                + ", % of renewable: "
                + str(np.round(sum_renew[y] * 100, 0))
                + "%"
                + ", GHG intensity per kWh: "
                + str(int(np.sum(co2_intensity_tech[y] * self.mix[y])))
                + " g. CO2-eq.",
                end=end_str,
            )

        if any(
            True for x in ["BEV", "PHEV-p", "PHEV-d"] if x in self.scope["powertrain"]
        ):
            for y in self.scope["year"]:
                index = self.get_index_vehicle_from_array(
                    ["BEV", "PHEV-p", "PHEV-d"], y, method="and"
                )

                self.A[
                    np.ix_(
                        np.arange(self.iterations),
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "electricity supply for electric vehicles" in i[0]
                        ],
                        [
                            self.inputs[i]
                            for i in self.inputs
                            if str(y) in i[0]
                            and "transport, passenger" in i[0]
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

            index = self.get_index_vehicle_from_array("FCEV")

            if "tertiary" in self.fuel_blends["hydrogen"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["hydrogen"]["primary"]["type"],
                        self.fuel_blends["hydrogen"]["secondary"]["type"],
                        self.fuel_blends["hydrogen"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    "{} is completed by {}.".format(
                        self.fuel_blends["hydrogen"]["primary"]["type"],
                        self.fuel_blends["hydrogen"]["secondary"]["type"],
                    ),
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["hydrogen"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["tertiary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["hydrogen"]["secondary"]["share"][y]
                                * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                # Fuel supply

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and "FCEV" in i[0]
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "fuel supply for hydrogen vehicles" in i[0]
                    ],
                    ind_A,
                ] = (
                    array[self.array_inputs["fuel mass"], :, ind_array]
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if "ICEV-g" in self.scope["powertrain"]:
            index = self.get_index_vehicle_from_array("ICEV-g")

            if "tertiary" in self.fuel_blends["cng"]:
                print(
                    "{} is completed by {} and {}.".format(
                        self.fuel_blends["cng"]["primary"]["type"],
                        self.fuel_blends["cng"]["secondary"]["type"],
                        self.fuel_blends["cng"]["tertiary"]["type"],
                    ),
                    end="\n \t * ",
                )

            else:

                print(
                    "{} is completed by {}.".format(
                        self.fuel_blends["cng"]["primary"]["type"],
                        self.fuel_blends["cng"]["secondary"]["type"],
                    ),
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["cng"]:
                    print(
                        "in "
                        + str(year)
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["secondary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%"
                        + " _________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["tertiary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )
                else:
                    print(
                        "in "
                        + str(year)
                        + " _________________________________________ "
                        + str(
                            np.round(
                                self.fuel_blends["cng"]["secondary"]["share"][y] * 100,
                                0,
                            )
                        )
                        + "%",
                        end=end_str,
                    )

                # Primary fuel share

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and "ICEV-g" in i[0]
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Supply of gas, including 4% of gas input as pump-to-tank leakage
                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0] and "fuel supply for gas vehicles" in i[0]
                    ],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        / array[self.array_inputs["range"], :, ind_array]
                    )
                    * (
                        1
                        + array[
                            self.array_inputs["CNG pump-to-tank leakage"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                # Gas leakage to air

                self.A[
                    :,
                    self.inputs[("Methane, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    (
                        array[self.array_inputs["fuel mass"], :, ind_array]
                        / array[self.array_inputs["range"], :, ind_array]
                    )
                    * array[self.array_inputs["CNG pump-to-tank leakage"], :, ind_array]
                    * -1
                ).T

                # Fuel-based emissions from CNG, CO2
                # The share and CO2 emissions factor of CNG is retrieved, if used

                share_fossil = 0
                CO2_fossil = 0

                if self.fuel_blends["cng"]["primary"]["type"] == "cng":
                    share_fossil = self.fuel_blends["cng"]["primary"]["share"][y]
                    CO2_fossil = (
                        self.fuel_blends["cng"]["primary"]["CO2"]
                        * self.fuel_blends["cng"]["primary"]["share"][y]
                    )

                if self.fuel_blends["cng"]["secondary"]["type"] == "cng":
                    share_fossil += self.fuel_blends["cng"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["cng"]["secondary"]["CO2"]
                        * self.fuel_blends["cng"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["cng"]:
                    if self.fuel_blends["cng"]["tertiary"]["type"] == "cng":
                        share_fossil += self.fuel_blends["cng"]["tertiary"]["share"][y]
                        CO2_fossil += (
                            self.fuel_blends["cng"]["tertiary"]["CO2"]
                            * self.fuel_blends["cng"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    array[self.array_inputs["fuel mass"], :, ind_array]
                    * CO2_fossil
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

                # Fuel-based CO2 emission from alternative CNG
                # The share of non-fossil gas in the blend is retrieved
                # As well as the CO2 emission factor of the fuel

                share_non_fossil = 0
                CO2_non_fossil = 0

                if self.fuel_blends["cng"]["primary"]["type"] != "cng":
                    share_non_fossil = self.fuel_blends["cng"]["primary"]["share"][y]
                    CO2_non_fossil = (
                        self.fuel_blends["cng"]["primary"]["CO2"]
                        * self.fuel_blends["cng"]["primary"]["share"][y]
                    )

                if self.fuel_blends["cng"]["secondary"]["type"] != "cng":
                    share_non_fossil += self.fuel_blends["cng"]["secondary"]["share"][y]
                    CO2_non_fossil += (
                        self.fuel_blends["cng"]["secondary"]["CO2"]
                        * self.fuel_blends["cng"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["cng"]:
                    if self.fuel_blends["cng"]["tertiary"]["type"] != "cng":
                        share_non_fossil += self.fuel_blends["cng"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["cng"]["tertiary"]["CO2"]
                            * self.fuel_blends["cng"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    (
                        (
                            array[self.array_inputs["fuel mass"], :, ind_array]
                            * CO2_non_fossil
                        )
                    )
                    / array[self.array_inputs["range"], :, ind_array]
                    * -1
                ).T

        if [i for i in self.scope["powertrain"] if i in ["ICEV-d", "PHEV-d", "HEV-d"]]:
            index = self.get_index_vehicle_from_array(["ICEV-d", "PHEV-d", "HEV-d"])

            if "tertiary" in self.fuel_blends["diesel"]:
                print(
                    f"{self.fuel_blends['diesel']['primary']['type']} "
                    f"is completed by {self.fuel_blends['diesel']['secondary']['type']} "
                    f"and {self.fuel_blends['diesel']['tertiary']['type']}.",
                    end="\n \t * ",
                )

            else:
                print(
                    f"{self.fuel_blends['diesel']['primary']['type']} "
                    f"is completed by {self.fuel_blends['diesel']['secondary']['type']}.",
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["diesel"]:

                    pct1 = np.round(
                        self.fuel_blends["diesel"]["secondary"]["share"][y] * 100,
                        0,
                    )
                    pct2 = np.round(
                        self.fuel_blends["diesel"]["tertiary"]["share"][y] * 100,
                        0,
                    )
                    print(
                        f"in {year} _________________ {pct1}% _________________ {pct2}%",
                        end=end_str,
                    )
                else:
                    pct1 = np.round(
                        self.fuel_blends["diesel"]["secondary"]["share"][y] * 100,
                        0,
                    )
                    print(
                        f"in {year} ___________________________________________ {pct1}%",
                        end=end_str,
                    )

                ind_A = [
                    self.inputs[i]
                    for i in self.inputs
                    if str(year) in i[0]
                    and "transport, passenger car" in i[0]
                    and any(x in i[0] for x in ["ICEV-d", "PHEV-d", "HEV-d"])
                ]

                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Fuel supply

                fuel_amount = (
                    array[
                        self.array_inputs["TtW energy, combustion mode"], :, ind_array
                    ]
                    / (
                        self.array.values[
                            self.array_inputs["LHV fuel MJ per kg"], :, ind_array
                        ]
                        * 1000
                    )
                    * (
                        1
                        - self.array.values[
                            self.array_inputs["electric utility factor"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                self.A[
                    :,
                    [
                        self.inputs[i]
                        for i in self.inputs
                        if str(year) in i[0]
                        and "fuel supply for diesel vehicles" in i[0]
                    ],
                    ind_A,
                ] = fuel_amount

                share_fossil = 0
                CO2_fossil = 0

                # Fuel-based CO2 emission from conventional diesel
                if self.fuel_blends["diesel"]["primary"]["type"] == "diesel":
                    share_fossil = self.fuel_blends["diesel"]["primary"]["share"][y]
                    CO2_fossil = (
                        self.fuel_blends["diesel"]["primary"]["CO2"]
                        * self.fuel_blends["diesel"]["primary"]["share"][y]
                    )

                if self.fuel_blends["diesel"]["secondary"]["type"] == "diesel":
                    share_fossil += self.fuel_blends["diesel"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["diesel"]["secondary"]["CO2"]
                        * self.fuel_blends["diesel"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["diesel"]:
                    if self.fuel_blends["diesel"]["tertiary"]["type"] == "diesel":
                        share_fossil += self.fuel_blends["diesel"]["tertiary"]["share"][
                            y
                        ]
                        CO2_fossil += (
                            self.fuel_blends["diesel"]["tertiary"]["CO2"]
                            * self.fuel_blends["diesel"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount * CO2_fossil
                )

                # Fuel-based SO2 emissions
                # Sulfur concentration value for a given country, a given year, as concentration ratio

                sulfur_concentration = self.get_sulfur_content(
                    self.country, "diesel", year
                )

                self.A[
                    :,
                    self.inputs[("Sulfur dioxide", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount
                    * sulfur_concentration
                    * (64 / 32)  # molar mass of SO2/molar mass of O2
                )

                share_non_fossil = 0
                CO2_non_fossil = 0

                # Fuel-based CO2 emission from alternative petrol
                # The share of non-fossil fuel in the blend is retrieved
                # As well as the CO2 emission factor of the fuel
                if self.fuel_blends["diesel"]["primary"]["type"] != "diesel":
                    share_non_fossil = self.fuel_blends["diesel"]["primary"]["share"][y]
                    CO2_non_fossil = (
                        self.fuel_blends["diesel"]["primary"]["CO2"]
                        * self.fuel_blends["diesel"]["primary"]["share"][y]
                    )

                if self.fuel_blends["diesel"]["secondary"]["type"] != "diesel":
                    share_non_fossil += self.fuel_blends["diesel"]["secondary"][
                        "share"
                    ][y]
                    CO2_non_fossil += (
                        self.fuel_blends["diesel"]["secondary"]["share"][y]
                        * self.fuel_blends["diesel"]["secondary"]["CO2"]
                    )

                if "tertiary" in self.fuel_blends["diesel"]:
                    if self.fuel_blends["diesel"]["tertiary"]["type"] != "diesel":
                        share_non_fossil += self.fuel_blends["diesel"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["diesel"]["tertiary"]["share"][y]
                            * self.fuel_blends["diesel"]["tertiary"]["CO2"]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    fuel_amount * CO2_non_fossil
                )

        if [i for i in self.scope["powertrain"] if i in ["ICEV-p", "HEV-p", "PHEV-p"]]:
            index = self.get_index_vehicle_from_array(["ICEV-p", "HEV-p", "PHEV-p"])

            if "tertiary" in self.fuel_blends["petrol"]:
                print(
                    f"{self.fuel_blends['petrol']['primary']['type']} "
                    f"is completed by {self.fuel_blends['petrol']['secondary']['type']} "
                    f"and {self.fuel_blends['petrol']['tertiary']['type']}.",
                    end="\n \t * ",
                )

            else:
                print(
                    f"{self.fuel_blends['petrol']['primary']['type']} "
                    f"is completed by {self.fuel_blends['petrol']['secondary']['type']}.",
                    end="\n \t * ",
                )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                if "tertiary" in self.fuel_blends["petrol"]:
                    pct1 = np.round(
                        self.fuel_blends["petrol"]["secondary"]["share"][y] * 100,
                        0,
                    )
                    pct2 = np.round(
                        self.fuel_blends["petrol"]["tertiary"]["share"][y] * 100,
                        0,
                    )
                    print(
                        f"in {year} _________________ {pct1}% _________________ {pct2}%",
                        end=end_str,
                    )
                else:
                    pct1 = np.round(
                        self.fuel_blends["petrol"]["secondary"]["share"][y] * 100,
                        0,
                    )
                    print(
                        f"in {year} ___________________________________________ {pct1}%",
                        end=end_str,
                    )

                ind_A = [
                    j
                    for i, j in self.inputs.items()
                    if str(year) in i[0]
                    and "transport, passenger" in i[0]
                    and any(x in i[0] for x in ["ICEV-p", "HEV-p", "PHEV-p"])
                ]
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Fuel supply

                fuel_amount = (
                    array[
                        self.array_inputs["TtW energy, combustion mode"], :, ind_array
                    ]
                    / (
                        self.array.values[
                            self.array_inputs["LHV fuel MJ per kg"], :, ind_array
                        ]
                        * 1000
                    )
                    * (
                        1
                        - self.array.values[
                            self.array_inputs["electric utility factor"], :, ind_array
                        ]
                    )
                    * -1
                ).T

                self.A[
                    :,
                    [
                        j
                        for i, j in self.inputs.items()
                        if str(year) in i[0]
                        and "fuel supply for gasoline vehicles" in i[0]
                    ],
                    ind_A,
                ] = fuel_amount

                share_fossil = 0
                CO2_fossil = 0

                # Fuel-based CO2 emission from conventional petrol
                if self.fuel_blends["petrol"]["primary"]["type"] == "petrol":
                    share_fossil = self.fuel_blends["petrol"]["primary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["petrol"]["primary"]["CO2"]
                        * self.fuel_blends["petrol"]["primary"]["share"][y]
                    )

                if self.fuel_blends["petrol"]["secondary"]["type"] == "petrol":
                    share_fossil += self.fuel_blends["petrol"]["secondary"]["share"][y]
                    CO2_fossil += (
                        self.fuel_blends["petrol"]["secondary"]["CO2"]
                        * self.fuel_blends["petrol"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["petrol"]:
                    if self.fuel_blends["petrol"]["tertiary"]["type"] == "petrol":
                        share_fossil += self.fuel_blends["petrol"]["tertiary"]["share"][
                            y
                        ]
                        CO2_fossil += (
                            self.fuel_blends["petrol"]["tertiary"]["CO2"]
                            * self.fuel_blends["petrol"]["tertiary"]["share"][y]
                        )

                self.A[
                    :,
                    self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount * CO2_fossil
                )

                # Fuel-based SO2 emissions
                # Sulfur concentration value for a given country, a given year, as a concentration ratio

                sulfur_concentration = self.get_sulfur_content(
                    self.country, "petrol", year
                )

                self.A[
                    :,
                    self.inputs[("Sulfur dioxide", ("air",), "kilogram")],
                    ind_A,
                ] = (
                    fuel_amount
                    * sulfur_concentration
                    * (64 / 32)  # molar mass of SO2/molar mass of O2
                )

                share_non_fossil = 0
                CO2_non_fossil = 0

                # Fuel-based CO2 emission from alternative petrol
                # The share of non-fossil fuel in the blend is retrieved
                # As well as the CO2 emission factor of the fuel
                if self.fuel_blends["petrol"]["primary"]["type"] != "petrol":
                    share_non_fossil += self.fuel_blends["petrol"]["primary"]["share"][
                        y
                    ]
                    CO2_non_fossil += (
                        self.fuel_blends["petrol"]["primary"]["CO2"]
                        * self.fuel_blends["petrol"]["primary"]["share"][y]
                    )

                if self.fuel_blends["petrol"]["secondary"]["type"] != "petrol":
                    share_non_fossil += self.fuel_blends["petrol"]["secondary"][
                        "share"
                    ][y]
                    CO2_non_fossil += (
                        self.fuel_blends["petrol"]["secondary"]["CO2"]
                        * self.fuel_blends["petrol"]["secondary"]["share"][y]
                    )

                if "tertiary" in self.fuel_blends["petrol"]:
                    if self.fuel_blends["petrol"]["tertiary"]["type"] != "petrol":
                        share_non_fossil += self.fuel_blends["petrol"]["tertiary"][
                            "share"
                        ][y]
                        CO2_non_fossil += (
                            self.fuel_blends["petrol"]["tertiary"]["share"][y]
                            * self.fuel_blends["petrol"]["tertiary"]["CO2"]
                        )

                self.A[
                    :,
                    self.inputs[
                        (
                            "Carbon dioxide, non-fossil",
                            ("air",),
                            "kilogram",
                        )
                    ],
                    ind_A,
                ] = (
                    fuel_amount * CO2_non_fossil
                )

        # Non-exhaust emissions

        # Road wear emissions + 33.3% of re-suspended road dust
        ind_A = [self.inputs[i] for i in self.inputs if "transport, passenger" in i[0]]
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
            ind_A,
        ] = array[self.array_inputs["road wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

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
            ind_A,
        ] = array[self.array_inputs["tire wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

        # Brake wear emissions
        # BEVs only emit 20% of what a combustion engine vehicle emit according to
        # https://link.springer.com/article/10.1007/s11367-014-0792-4
        ind_A = [self.inputs[i] for i in self.inputs if "transport, passenger" in i[0]]

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
            ind_A,
        ] = array[self.array_inputs["brake wear emissions"], :] + (
            0.333 * array[self.array_inputs["road dust emissions"], :]
        )

        # Infrastructure
        self.A[
            :,
            self.inputs[("market for road", "GLO", "meter-year", "road")],
            [self.inputs[i] for i in self.inputs if "transport, passenger car" in i[0]],
        ] = (
            5.37e-7 * array[self.array_inputs["driving mass"], :] * -1
        )

        # Infrastructure maintenance
        self.A[
            :,
            self.inputs[
                ("market for road maintenance", "RER", "meter-year", "road maintenance")
            ],
            [self.inputs[i] for i in self.inputs if "transport, passenger car" in i[0]],
        ] = (
            1.29e-3 * -1
        )

        # Exhaust emissions
        # Non-fuel based emissions

        self.A[
            np.ix_(
                np.arange(self.iterations),
                self.index_emissions,
                [
                    j
                    for i, j in self.inputs.items()
                    if "transport, passenger car" in i[0]
                ],
            )
        ] = (
            array[
                [
                    self.array_inputs[self.map_fuel_emissions[self.rev_inputs[x]]]
                    for x in self.index_emissions
                ]
            ]
            * -1
        ).transpose(
            [1, 0, 2]
        )

        # Noise emissions
        self.A[
            np.ix_(
                np.arange(self.iterations),
                self.index_noise,
                [
                    j
                    for i, j in self.inputs.items()
                    if "transport, passenger car" in i[0]
                ],
            )
        ] = (
            array[
                [
                    self.array_inputs[self.map_noise_emissions[self.rev_inputs[x]]]
                    for x in self.index_noise
                ]
            ]
            * -1
        ).transpose(
            [1, 0, 2]
        )

        # Emissions of air conditioner refrigerant r134a
        # Leakage assumed to amount to 53g according to
        # https://treeze.ch/fileadmin/user_upload/downloads/Publications/Case_Studies/Mobility/544-LCI-Road-NonRoad-Transport-Services-v2.0.pdf
        # but only to cars with an AC system (meaning, with a cooling energy consumption)

        self.A[
            :,
            self.inputs[
                ("Ethane, 1,1,1,2-tetrafluoro-, HFC-134a", ("air",), "kilogram")
            ],
            [j for i, j in self.inputs.items() if "transport, passenger car" in i[0]],
        ] = (
            0.75 / self.array.values[self.array_inputs["lifetime kilometers"]] * -1
        ) * self.array.values[
            self.array_inputs["cooling energy consumption"]
        ]

        self.A[
            :,
            self.inputs[
                ("market for refrigerant R134a", "GLO", "kilogram", "refrigerant R134a")
            ],
            [j for i, j in self.inputs.items() if "transport, passenger car" in i[0]],
        ] = (
            (0.75 + 0.55)
            / self.array.values[self.array_inputs["lifetime kilometers"]]
            * -1
        ) * self.array.values[
            self.array_inputs["cooling energy consumption"]
        ]

        print("*********************************************************************")

    def select_heat_supplier(self, heat_supplier):
        """
        The heat supply is an important aspect of direct air capture.
        Here, we can change the supplier of heat.
        :param heat_supplier: by default "waste heat". Must be one of "waste heat", "biomass heat",
        "natural gas heat", "market heat".
        :type heat_supplier: str
        :return:
        """

        d_heat_suppliers = {
            "waste heat": (
                "heat, from municipal waste incineration to generic market for heat district or industrial, other than natural gas",
                "CH",
                "megajoule",
                "heat, district or industrial, other than natural gas",
            ),
            "biomass heat": (
                "heat production, hardwood chips from forest, at furnace 1000kW, state-of-the-art 2014",
                "CH",
                "megajoule",
                "heat, district or industrial, other than natural gas",
            ),
            "natural gas heat": (
                "market group for heat, central or small-scale, natural gas",
                "RER",
                "megajoule",
                "heat, central or small-scale, natural gas",
            ),
            "market heat": (
                "steam production, as energy carrier, in chemical industry",
                "RoW",
                "megajoule",
                "heat, from steam, in chemical industry",
            ),
        }

        air_capture = self.inputs[
            (
                "carbon dioxide, captured from atmosphere",
                "RER",
                "kilogram",
                "carbon dioxide, captured from the atmosphere",
            )
        ]

        methanol_distillation = self.inputs[
            (
                "methanol distillation, hydrogen from electrolysis, CO2 from DAC",
                "RER",
                "kilogram",
                "methanol, purified",
            )
        ]

        all_inds = [self.inputs[i] for i in list(d_heat_suppliers.values())]

        # amount of heat for DAC
        heat_amount = (
            self.A[np.ix_(range(self.A.shape[0]), all_inds, [air_capture])] * -1
        ).max() * -1

        # zero out the heat input
        self.A[np.ix_(range(self.A.shape[0]), all_inds, [air_capture])] = 0

        if heat_supplier == "heat pump":

            # we convert the need for heat into electricity
            # later on, this need for electricity will be provided
            # by the country-specific mix
            # we assume a CoP of 2.9

            id_elec = self.inputs[
                (
                    "market group for electricity, low voltage",
                    "ENTSO-E",
                    "kilowatt hour",
                    "electricity, low voltage",
                )
            ]

            self.A[np.ix_(range(self.A.shape[0]), [id_elec], [air_capture])] = (
                heat_amount / 2.9 / 3.6
            )

        else:
            # find index of the new supplier and set the amount
            ind = self.inputs[d_heat_suppliers[heat_supplier]]
            self.A[np.ix_(range(self.A.shape[0]), [ind], [air_capture])] = heat_amount

        # Methanol distillation
        heat_amount = (
            self.A[np.ix_(range(self.A.shape[0]), all_inds, [methanol_distillation])]
            * -1
        ).max() * -1

        # zero out the heat input
        self.A[np.ix_(range(self.A.shape[0]), all_inds, [methanol_distillation])] = 0

        if heat_supplier == "heat pump":
            # we convert the need for heat into electricity
            # later on, this need for electricity will be provided
            # by the country-specific mix
            # we assume a CoP of 2.9

            id_elec = self.inputs[
                (
                    "market group for electricity, low voltage",
                    "ENTSO-E",
                    "kilowatt hour",
                    "electricity, low voltage",
                )
            ]

            self.A[
                np.ix_(range(self.A.shape[0]), [id_elec], [methanol_distillation])
            ] = (heat_amount / 2.9 / 3.6)

        else:
            # find index of the new supplier and set the amount
            ind = self.inputs[d_heat_suppliers[heat_supplier]]
            self.A[
                np.ix_(range(self.A.shape[0]), [ind], [methanol_distillation])
            ] = heat_amount
