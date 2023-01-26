"""
inventory.py contains InventoryCalculation which provides all methods to solve inventories.
"""

import csv
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import xarray as xr
import yaml
from scipy import sparse

from . import DATA_DIR
from .background_systems import BackgroundSystemModel
from .export import ExportInventory

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

IAM_FILES_DIR = DATA_DIR / "IAM"


RANGE_PARAM = {
    "two-wheeler": "range",
    "car": "range",
    "bus": "daily distance",
    "truck": "target range",
}


def check_func_unit(func_unit):
    """Check if func_unit is a valid functional unit."""
    if func_unit not in ["vkm", "pkm", "tkm"]:
        raise ValueError(
            f"Functional unit must be one of "
            f"'vkm', 'pkm', 'tkm', "
            f"not {func_unit}"
        )
    return func_unit


def check_scenario(scenario):
    """Check if scenario is a valid scenario."""
    valid_scenarios = ["SSP2-NPi", "SSP2-PkBudg1150", "SSP2-PkBudg500", "static"]
    if scenario not in valid_scenarios:
        raise ValueError(
            f"Scenario must be one of " f"{valid_scenarios}, " f"not {scenario}"
        )
    return scenario


def get_noise_emission_flows() -> dict:
    """Get noise emission flows from the noise emission file."""

    return {
        (
            f"noise, octave {i}, day time, {comp}",
            (f"octave {i}", "day time", comp),
            "joule",
        ): f"noise, octave {i}, day time, {comp}"
        for i in range(1, 9)
        for comp in ["urban", "suburban", "rural"]
    }


def get_exhaust_emission_flows() -> dict:
    with open(
        DATA_DIR / "emission_factors" / "exhaust_and_noise_flows.yaml",
        "r",
        encoding="utf-8",
    ) as stream:
        flows = yaml.safe_load(stream)["exhaust"]

    d_comp = {
        "urban": "urban air close to ground",
        "suburban": "non-urban air or from high stacks",
        "rural": "low population density, long-term",
    }

    return {
        (v, ("air", d_comp[comp]), "kilogram"): f"{k} direct emissions, {comp}"
        for k, v in flows.items()
        for comp in ["urban", "suburban", "rural"]
    }


def get_dict_impact_categories(method, indicator) -> dict:
    """
    Load a dictionary with available impact assessment
    methods as keys, and assessment level and categories as values.

    :return: dictionary
    :rtype: dict
    """
    filename = "dict_impact_categories.csv"
    filepath = DATA_DIR / "lcia" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of impact categories could not be found."
        )

    csv_dict = {}

    with open(filepath, encoding="utf-8") as f:
        input_dict = csv.reader(f, delimiter=";")
        for row in input_dict:
            if row[0] == method and row[3] == indicator:
                csv_dict[row[2]] = {
                    "method": row[1],
                    "category": row[2],
                    "type": row[3],
                    "abbreviation": row[4],
                    "unit": row[5],
                    "source": row[6],
                }

    return csv_dict


def get_dict_input() -> dict:
    """
    Load a dictionary with tuple ("name of activity", "location", "unit",
    "reference product") as key, row/column
    indices as values.

    :return: dictionary with `label:index` pairs.
    :rtype: dict

    """
    filename = f"dict_inputs_A_matrix.csv"
    filepath = DATA_DIR / "IAM" / filename
    if not filepath.is_file():
        raise FileNotFoundError("The dictionary of activity labels could not be found.")

    with open(filepath, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=";")
        raw = list(reader)
        for _r, r in enumerate(raw):

            if len(r) == 3:
                r[1] = eval(r[1])
            raw[_r] = tuple(r)

        return {j: i for i, j in enumerate(list(raw))}


class Inventory:
    """
    Build and solve the inventory for results characterization and inventory export

    :ivar vm: object from the VehicleModel class
    :ivar background_configuration: dictionary that contains choices for background system
    :ivar scenario: IAM energy scenario to use (
        "SSP2-NPi": Nationally implemented policies, limits temperature increase by 2100 to 3.3 degrees Celsius,
        "SSP2-PkBudg1150": limits temperature increase by 2100 to 2 degrees Celsius,
        "SSP2-PkBudg500": limits temperature increase by 2100 to 1.5 degrees Celsius,
        "static": no forward-looking modification of the background inventories).
        "SSP2-NPi" selected by default.)

    """

    def __init__(
        self,
        vm,
        background_configuration: dict = None,
        scenario: str = "SSP2-NPi",
        method: str = "recipe",
        indicator: str = "midpoint",
        functional_unit: str = "vkm",
    ) -> None:

        self.vm = vm

        self.scope = {
            "size": vm.array.coords["size"].values.tolist(),
            "powertrain": vm.array.coords["powertrain"].values.tolist(),
            "year": vm.array.coords["year"].values.tolist(),
        }
        self.scenario = check_scenario(scenario)
        self.func_unit = check_func_unit(functional_unit)

        self.method = method
        self.indicator = indicator if method == "recipe" else "midpoint"

        self.array = vm.array.stack(desired=["size", "powertrain", "year"])
        self.iterations = len(vm.array.value.values)

        self.number_of_vehicles = (self.vm["TtW energy"] > 0).sum().values

        self.array_inputs = {
            x: i for i, x in enumerate(list(self.array.parameter.values))
        }
        self.array_powertrains = {
            x: i for i, x in enumerate(list(self.array.powertrain.values))
        }

        self.background_configuration = {}
        self.background_configuration.update(background_configuration or {})

        self.inputs = get_dict_input()

        self.bs = BackgroundSystemModel()
        self.add_additional_activities()
        self.rev_inputs = {v: k for k, v in self.inputs.items()}

        with open(
            DATA_DIR / "electricity" / "elec_tech_map.yaml", "r", encoding="utf-8"
        ) as stream:
            self.elec_map = yaml.safe_load(stream)
            self.elec_map = {k: tuple(v) for k, v in self.elec_map.items()}

        self.electricity_technologies = list(self.elec_map.keys())

        self.A = self.get_A_matrix()
        # Create electricity and fuel market datasets
        self.mix = self.define_electricity_mix_for_fuel_prep()
        self.create_electricity_mix_for_fuel_prep()
        # Create electricity market dataset for battery production
        self.create_electricity_market_for_battery_production()
        self.rev_inputs = {v: k for k, v in self.inputs.items()}
        self.create_fuel_markets()

        self.exhaust_emissions = get_exhaust_emission_flows()
        self.noise_emissions = get_noise_emission_flows()

        self.list_cat, self.split_indices = self.get_split_indices()

        self.impact_categories = get_dict_impact_categories(
            method=self.method, indicator=self.indicator
        )

        # Create the B matrix
        self.B = self.get_B_matrix()
        self.rev_inputs = {v: k for k, v in self.inputs.items()}

        self.fill_in_A_matrix()
        self.remove_uncompliant_vehicles()

    def get_results_table(self, sensitivity: bool = False) -> xr.DataArray:
        """
        Format a xarray.DataArray array to receive the results.

        :param sensitivity: if True, the results table will
        be formatted to receive sensitivity analysis results
        :return: xarrray.DataArray
        """

        params = [a for a in self.array.value.values]
        response = xr.DataArray(
            np.zeros(
                (
                    len(self.impact_categories),
                    len(self.scope["size"]),
                    len(self.scope["powertrain"]),
                    len(self.scope["year"]),
                    len(self.list_cat),
                    self.iterations,
                )
            ),
            coords=[
                list(self.impact_categories.keys()),
                self.scope["size"],
                self.scope["powertrain"],
                self.scope["year"],
                self.list_cat,
                np.arange(0, self.iterations) if not sensitivity else params,
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

        if sensitivity:
            # remove the `impact` dimension
            response = response.squeeze("impact")

        return response

    def get_split_indices(self):
        """
        Return list of indices to split the results into categories.

        :return: list of indices
        :rtype: list
        """
        # read `impact_source_categories.yml` file
        with open(
            DATA_DIR / "lcia" / "impact_source_categories.yml", "r", encoding="utf-8"
        ) as stream:
            source_cats = yaml.safe_load(stream)

        idx_cats = defaultdict(list)

        for cat, name in source_cats.items():
            for n in name:
                if self.find_input_indices((n,)):
                    idx_cats[cat].extend(self.find_input_indices((n,)))

        # add flows corresponding to `exhaust - direct`
        idx_cats["direct - exhaust"] = [
            self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
            self.inputs[("Carbon dioxide, non-fossil", ("air",), "kilogram")],
        ]
        idx_cats["direct - exhaust"].extend(
            [self.inputs[i] for i in self.exhaust_emissions]
        )
        idx_cats["direct - exhaust"].extend(
            [self.inputs[i] for i in self.noise_emissions]
        )

        idx_cats["direct - non-exhaust"].append(
            self.inputs[("Methane, fossil", ("air",), "kilogram")],
        )

        # idx for an input that has no burden
        # oxygen in this case
        extra_idx = [j for i, j in self.inputs.items() if i[0].lower() == "oxygen"][0]

        list_ind = [val for val in idx_cats.values()]
        maxLen = max(map(len, list_ind))
        for row in list_ind:
            while len(row) < maxLen:
                row.append(extra_idx)

        return list(idx_cats.keys()), list_ind

    def get_load_factor(self):

        # If the FU is in passenger-km, we normalize the results by
        # the number of passengers
        if self.func_unit == "vkm":
            load_factor = 1
        elif self.func_unit == "pkm":
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
        else:
            # ton kilometers
            load_factor = np.resize(
                self.array[self.array_inputs["cargo mass"]].values / 1000,
                (
                    1,
                    len(self.scope["size"]),
                    len(self.scope["powertrain"]),
                    len(self.scope["year"]),
                    1,
                    1,
                ),
            )

        return load_factor

    def calculate_impacts(self, sensitivity=False):

        if self.scenario != "static":
            B = self.B.interp(
                year=self.scope["year"], kwargs={"fill_value": "extrapolate"}
            ).values
        else:
            B = self.B[0].values

        # Prepare an array to store the results
        results = self.get_results_table(sensitivity=sensitivity)

        new_arr = np.zeros((self.A.shape[1], self.B.shape[1], len(self.scope["year"])))

        f_vector = np.zeros((np.shape(self.A)[1]))

        # Collect indices of activities contributing to the first level
        idx_car_trspt = self.find_input_indices(
            (f"transport, {self.vm.vehicle_type}, ",)
        )
        idx_cars = self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",))

        idx_others = [
            i
            for i in range(len(self.inputs))
            if i not in idx_car_trspt and i not in idx_cars
        ]

        arr = (
            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    idx_others,
                    np.concatenate([idx_cars, idx_car_trspt]),
                )
            ]
            .sum(axis=0)
            .sum(axis=1)
        )

        nonzero_idx = np.argwhere(arr)

        for a in nonzero_idx:
            f_vector[:] = 0
            f_vector[a] = 1
            X = sparse.linalg.spsolve(sparse.csr_matrix(self.A[0]), f_vector.T)

            new_arr[a] = (X * B).sum(axis=-1).T

        new_arr = new_arr.transpose(1, 0, 2)

        arr = (
            self.A[:, :, idx_car_trspt].reshape(
                self.iterations,
                -1,
                len(self.scope["size"]),
                len(self.scope["powertrain"]),
                len(self.scope["year"]),
            )
            * new_arr[:, None, :, None, None, :]
            * -1
        )

        arr += (
            self.A[:, :, idx_cars].reshape(
                self.iterations,
                -1,
                len(self.scope["size"]),
                len(self.scope["powertrain"]),
                len(self.scope["year"]),
            )
            * new_arr[:, None, :, None, None, :]
            * self.A[:, idx_cars, idx_car_trspt].reshape(
                self.iterations,
                -1,
                len(self.scope["size"]),
                len(self.scope["powertrain"]),
                len(self.scope["year"]),
            )
        )

        arr = arr[:, :, self.split_indices].sum(axis=3)

        # reshape the array to match the dimensions of the results table
        arr = arr.transpose(0, 3, 4, 5, 2, 1)

        if sensitivity:
            results[...] = arr.sum(axis=-2)
            results /= results.sel(parameter="reference")
        else:
            results[...] = arr

        return results / self.get_load_factor()

    def add_additional_activities(self):
        # Add as many rows and columns as cars to consider
        # Also add additional columns and rows for electricity markets
        # for fuel preparation and energy battery production

        maximum = max(self.inputs.values())

        for year in self.scope["year"]:
            for fuel in ["petrol", "diesel", "hydrogen", "cng"]:
                maximum += 1
                self.inputs[
                    (
                        f"fuel supply for {fuel} vehicles, {year}",
                        self.vm.country,
                        "kilogram",
                        "fuel",
                    )
                ] = maximum

            for electricity_source in [
                (f"electricity supply for electric vehicles, {year}", self.vm.country),
                (f"electricity supply for fuel preparation, {year}", self.vm.country),
                (
                    f"electricity supply for battery production, {year}",
                    self.vm.energy_storage["origin"],
                ),
            ]:
                maximum += 1
                self.inputs[
                    (
                        electricity_source[0],
                        electricity_source[1],
                        "kilowatt hour",
                        "electricity, low voltage",
                    )
                ] = maximum

        with open(DATA_DIR / "emission_factors" / "euro_classes.yaml", "r") as stream:
            euro_classes = yaml.safe_load(stream)[self.vm.vehicle_type]

        list_years = np.clip(
            self.scope["year"],
            min(euro_classes.keys()),
            max(euro_classes.keys()),
        )

        list_euro_classes = [euro_classes[y] for y in list(list_years)]

        for size in self.scope["size"]:
            for powertrain in self.scope["powertrain"]:
                for euro_class, year in zip(list_euro_classes, self.scope["year"]):

                    maximum += 1

                    if self.func_unit == "vkm":
                        unit = "kilometer"
                    elif self.func_unit == "pkm":
                        unit = "passenger kilometer"
                    else:
                        unit = "ton kilometer"

                    if powertrain in ["BEV", "BEV-opp", "BEV-depot", "BEV-motion"]:
                        chemistry = self.vm.energy_storage["electric"][
                            (powertrain, size, year)
                        ]
                        name = f"transport, {self.vm.vehicle_type}, {powertrain}, {chemistry} battery, {size}, {year}"
                        ref = f"transport, {self.vm.vehicle_type}"

                    elif powertrain in ["FCEV", "Human"]:
                        name = f"transport, {self.vm.vehicle_type}, {powertrain}, {size}, {year}"
                        ref = f"transport, {self.vm.vehicle_type}"

                    else:
                        name = f"transport, {self.vm.vehicle_type}, {powertrain}, {size}, {year}, Euro-{euro_class}"
                        ref = f"transport, {self.vm.vehicle_type}, Euro-{euro_class}"

                    # add transport activity
                    self.inputs[(name, self.vm.country, unit, ref)] = maximum
                    maximum += 1
                    # add vehicle
                    self.inputs[
                        (
                            name.replace(
                                f"transport, {self.vm.vehicle_type}",
                                self.vm.vehicle_type.capitalize(),
                            ),
                            self.vm.country,
                            "unit",
                            ref.replace(
                                f"transport, {self.vm.vehicle_type}",
                                self.vm.vehicle_type.capitalize(),
                            ),
                        )
                    ] = maximum

    def get_A_matrix(self):
        """
        Load the A matrix. The matrix contains exchanges of products (rows)
        between activities (columns).

        :return: A matrix with three dimensions of shape (number of values,
        number of products, number of activities).
        :rtype: numpy.ndarray

        """

        filename = "A_matrix.csv"
        filepath = DATA_DIR / "IAM" / filename
        if not filepath.is_file():
            raise FileNotFoundError("The IAM files could not be found.")

        # build matrix A from coordinates
        A_coords = np.genfromtxt(filepath, delimiter=";")
        I = A_coords[:, 0].astype(int)
        J = A_coords[:, 1].astype(int)
        initial_A = sparse.csr_matrix((A_coords[:, 2], (I, J))).toarray()
        new_A = np.identity(len(self.inputs))
        new_A[0 : np.shape(initial_A)[0], 0 : np.shape(initial_A)[0]] = initial_A

        # Resize the matrix to fit the number of iterations in `array`
        new_A = np.resize(
            new_A,
            (self.array.shape[1], new_A.shape[0], new_A.shape[1]),
        )
        return new_A

    def get_B_matrix(self) -> xr.DataArray:
        """
        Load the B matrix. The B matrix contains impact assessment
        figures for a give impact assessment method,
        per unit of activity. Its length column-wise equals
        the length of the A matrix row-wise.
        Its length row-wise equals the number of
        impact assessment methods.

        :return: an array with impact values per unit
        of activity for each method.
        :rtype: numpy.ndarray

        """

        filepaths = [
            str(fp)
            for fp in list(Path(IAM_FILES_DIR).glob("*.csv"))
            if all(x in str(fp) for x in [self.method, self.indicator, self.scenario])
        ]

        filepaths = sorted(filepaths, key=lambda x: int(re.findall("\d+", x)[1]))

        B = np.zeros((len(filepaths), len(self.impact_categories), len(self.inputs)))

        for f, filepath in enumerate(filepaths):
            initial_B = np.genfromtxt(filepath, delimiter=";")
            new_B = np.zeros(
                (
                    initial_B.shape[0],
                    len(self.inputs),
                )
            )

            new_B[0 : initial_B.shape[0], 0 : initial_B.shape[1]] = initial_B
            B[f, :, :] = new_B

        return xr.DataArray(
            B,
            coords=[
                [2005, 2010, 2020, 2030, 2040, 2050]
                if self.scenario != "static"
                else [2020],
                list(self.impact_categories.keys()),
                list(self.inputs.keys()),
            ],
            dims=["year", "category", "activity"],
        )

    def get_index_vehicle_from_array(
        self, items_to_look_for, items_to_look_for_also=[], method="or"
    ):
        """
        Return list of row/column indices of self.vm.array of labels
        that contain the string defined in `items_to_look_for`.

        :param items_to_look_for_also:
        :param method:
        :param items_to_look_for: string to search for
        :return: list
        """
        if not isinstance(items_to_look_for, list):
            items_to_look_for = [items_to_look_for]

        if not isinstance(items_to_look_for_also, list):
            items_to_look_for_also = [items_to_look_for_also]

        list_vehicles = self.array.desired.values

        return (
            [
                c
                for c, v in enumerate(list_vehicles)
                if set(items_to_look_for).intersection(v)
            ]
            if method == "or"
            else [
                c
                for c, v in enumerate(list_vehicles)
                if set(items_to_look_for).intersection(v)
                and set(items_to_look_for_also).intersection(v)
            ]
        )

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

    def define_electricity_mix_for_fuel_prep(self) -> np.ndarray:
        """
        This function defines a fuel mix based either on user-defined mix,
        or on default mixes for a given country.
        The mix is calculated as the average mix, weighted by the
        distribution of annually driven kilometers.
        :return:
        """
        try:
            losses_to_low = float(self.bs.losses[self.vm.country]["LV"])
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

            if self.vm.country not in self.bs.electricity_mix.country.values:
                print(
                    f"The electricity mix for {self.vm.country} could not be found."
                    "Average European electricity mix is used instead."
                )
                country = "RER"
            else:
                country = self.vm.country

            mix = [
                self.bs.electricity_mix.sel(
                    country=country,
                    variable=self.electricity_technologies,
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
                    variable=self.electricity_technologies,
                )
                .interp(
                    year=np.arange(year, 2051), kwargs={"fill_value": "extrapolate"}
                )
                .mean(axis=0)
                .values
                for y, year in enumerate(self.scope["year"])
            ]

        return np.clip(mix, 0, 1) / np.clip(mix, 0, 1).sum(axis=1)[:, None]

    def define_renewable_rate_in_mix(self) -> list:
        """
        This function returns the renewable rate in the electricity mix
        for each year.
        """

        sum_renew = [
            np.sum([mix[i] for i in [0, 3, 4, 5, 8, 9, 10, 11, 14, 15, 18, 19]])
            for mix in self.mix
        ]

        return sum_renew

    # @lru_cache
    def find_input_indices(self, contains: [tuple, str], excludes: tuple = ()) -> list:
        """
        This function finds the indices of the inputs in the A matrix
        that contain the strings in the contains list, and do not
        contain the strings in the excludes list.
        :param contains: list of strings
        :param excludes: list of strings
        :return: list of indices
        """
        indices = []

        if not isinstance(contains, tuple):
            contains = tuple(contains)

        if not isinstance(excludes, tuple):
            excludes = tuple(excludes)

        for i, input in enumerate(self.inputs):
            if all([c in input[0] for c in contains]) and not any(
                [e in input[0] for e in excludes]
            ):
                indices.append(i)

        return indices

    def add_electricity_infrastructure(self, dataset, losses):

        for y, year in enumerate(self.scope["year"]):
            # Add transmission network for high and medium voltage
            for input in [
                (
                    "transmission network construction, electricity, high voltage",
                    dataset,
                    6.58e-9 * -1 * losses,
                ),
                (
                    "transmission network construction, electricity, medium voltage",
                    dataset,
                    1.86e-8 * -1 * losses,
                ),
                (
                    "transmission network construction, long-distance",
                    dataset,
                    3.17e-10 * -1 * losses,
                ),
                (
                    "distribution network construction, electricity, low voltage",
                    dataset,
                    8.74e-8 * -1 * losses,
                ),
                (
                    "market for sulfur hexafluoride, liquid",
                    dataset,
                    (5.4e-8 + 2.99e-9) * -1 * losses,
                ),
                ("Sulfur hexafluoride", dataset, (5.4e-8 + 2.99e-9) * -1 * losses),
            ]:
                self.A[
                    np.ix_(
                        np.arange(self.iterations),
                        self.find_input_indices(
                            input[0],
                        ),
                        self.find_input_indices((input[1], str(year))),
                    )
                ] = input[2]

    def create_electricity_mix_for_fuel_prep(self):
        """
        This function fills the electricity market that
        supplies battery charging operations
        and hydrogen production through electrolysis.
        """

        try:
            losses_to_low = float(self.bs.losses[self.vm.country]["LV"])
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
                    [
                        self.inputs[self.elec_map[t]]
                        for t in self.electricity_technologies
                    ],
                    self.find_input_indices(
                        ("electricity supply for fuel preparation", str(year))
                    ),
                )
            ] = (
                m * -1 * losses_to_low
            )

        self.add_electricity_infrastructure(
            "electricity supply for fuel preparation", losses_to_low
        )

    def create_electricity_market_for_battery_production(self):
        """
        This function fills in the column in `self.A` concerned
        with the electricity mix used for manufacturing battery cells
        :return:
        """

        battery_origin = self.vm.energy_storage.get("origin", "CN")

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
                    variable=self.electricity_technologies,
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
                    [
                        self.inputs[self.elec_map[t]]
                        for t in self.electricity_technologies
                    ],
                    self.find_input_indices(
                        ("electricity supply for battery production", str(year))
                    ),
                )
            ] = (
                m * losses_to_low * -1
            )

        self.add_electricity_infrastructure(
            "electricity supply for battery production", losses_to_low
        )

    def get_sulfur_content(self, location, fuel, year):
        """
        Return the sulfur content in the fuel.
        If a region is passed, the average sulfur content over
        the countries the region contains is returned.

        :param year:
        :param location: str. A country or region ISO code
        :param fuel: str. "diesel" or "petrol"
        :return: float. Sulfur content in ppm.
        """

        if fuel not in self.bs.sulfur.fuel.values:
            return 0

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

    def learning_rate_fuel(self, fuel, year, share, val):
        """
        This function calculates the learning rate for hydrogen
        production by electrolysis.
        """

        amount_h2 = {
            "electrolysis": 1,
            "synthetic gasoline - energy allocation": 0.338,
            "synthetic gasoline - economic allocation": 0.6385,
            "synthetic diesel - energy allocation": 0.42,
            "synthetic diesel - economic allocation": 0.183,
        }
        year = float(year)
        electrolysis_electricity = -0.3538 * (year - 2010) + 58.589
        electricity = val - (amount_h2.get(fuel, 0) * 58)
        electricity += electrolysis_electricity * amount_h2.get(fuel, 0)
        electricity *= share

        return electricity

    def create_fuel_markets(self):
        """
        This function creates markets for fuel, considering a given blend,
        a given fuel type and a given year.
        It also adds separate electricity input in case hydrogen
        from electrolysis is needed somewhere in the fuel supply chain.
        :return:
        """

        d_dataset_name = {
            "petrol": "fuel supply for petrol vehicles, ",
            "diesel": "fuel supply for diesel vehicles, ",
            "cng": "fuel supply for cng vehicles, ",
            "hydrogen": "fuel supply for hydrogen vehicles, ",
        }

        # electricity dataset
        for y, year in enumerate(self.scope["year"]):
            self.A[
                :,
                self.find_input_indices(
                    ("electricity supply for fuel preparation, ", str(year))
                ),
                self.find_input_indices(
                    (f"electricity supply for electric vehicles, {year}",)
                ),
            ] = -1

        # fuel datasets
        for fuel_type in self.vm.fuel_blend:
            primary = self.vm.fuel_blend[fuel_type]["primary"]["type"]
            secondary = self.vm.fuel_blend[fuel_type]["secondary"]["type"]

            electricity_requirement_primary = self.find_input_requirement(
                value_in="kilowatt hour",
                value_out=self.vm.fuel_blend[fuel_type]["primary"]["name"][0],
                find_input_by="unit",
                zero_out_input=True,
            )
            electricity_requirement_secondary = self.find_input_requirement(
                value_in="kilowatt hour",
                value_out=self.vm.fuel_blend[fuel_type]["secondary"]["name"][0],
                find_input_by="unit",
                zero_out_input=True,
            )

            for y, year in enumerate(self.scope["year"]):
                primary_share = self.vm.fuel_blend[fuel_type]["primary"]["share"][y]
                secondary_share = self.vm.fuel_blend[fuel_type]["secondary"]["share"][y]
                dataset_name = f"{d_dataset_name[fuel_type]}{year}"
                fuel_market_index = self.find_input_indices((dataset_name,))

                try:
                    primary_fuel_activity_index = self.inputs[
                        self.vm.fuel_blend[fuel_type]["primary"]["name"]
                    ]
                    secondary_fuel_activity_index = self.inputs[
                        self.vm.fuel_blend[fuel_type]["secondary"]["name"]
                    ]
                except KeyError:
                    raise KeyError(
                        "One of the primary or secondary fuels specified in "
                        "the fuel blend for {} is not valid.".format(fuel_type)
                    )

                self.A[:, primary_fuel_activity_index, fuel_market_index] = (
                    -1 * primary_share
                )
                self.A[:, secondary_fuel_activity_index, fuel_market_index] = (
                    -1 * secondary_share
                )

                additional_electricity_primary = self.learning_rate_fuel(
                    primary, year, primary_share, electricity_requirement_primary
                )

                additional_electricity_secondary = self.learning_rate_fuel(
                    secondary,
                    year,
                    secondary_share,
                    electricity_requirement_secondary,
                )

                additional_electricity = (
                    additional_electricity_primary + additional_electricity_secondary
                )

                if additional_electricity > 0:
                    elec_idx = self.find_input_indices(
                        (f"electricity supply for fuel preparation, {year}",)
                    )

                    self.A[:, elec_idx, fuel_market_index] = -1 * additional_electricity

    def find_input_requirement(
        self,
        value_in,
        value_out,
        find_input_by="name",
        zero_out_input=False,
        filter_activities=None,
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
            value_out = (value_out,)

        index_output = self.find_input_indices(value_out)

        f_vector = np.zeros((np.shape(self.A)[1]))
        f_vector[index_output] = 1

        X = sparse.linalg.spsolve(sparse.csr_matrix(self.A[0]), f_vector.T)

        ind_inputs = np.nonzero(X)[0]

        if find_input_by == "name":
            ins = [
                i
                for i in ind_inputs
                if value_in.lower() in self.rev_inputs[i][0].lower()
            ]

        elif find_input_by == "unit":

            ins = [
                i
                for i in ind_inputs
                if value_in.lower() in self.rev_inputs[i][2].lower()
            ]
        else:
            raise ValueError("find_input_by must be 'name' or 'unit'")

        outs = [i for i in ind_inputs if i not in ins]

        if filter_activities:
            outs = [
                i
                for e in filter_activities
                for i in outs
                if e.lower() in self.rev_inputs[i][0].lower()
            ]

        ins = [
            i
            for i in ins
            if self.A[np.ix_(np.arange(0, self.A.shape[0]), [i], outs)].sum() != 0
        ]

        sum_supplied = X[ins].sum()

        if zero_out_input:
            # zero out initial inputs
            self.A[np.ix_(np.arange(0, self.A.shape[0]), ins, outs)] = 0

        return sum_supplied

    def get_fuel_blend_carbon_intensity(
        self, fuel_type: str
    ) -> [np.ndarray, np.ndarray]:
        """
        Returns the carbon intensity of a fuel blend.
        :param fuel_type: fuel type
        :return: carbon intensity of fuel blend fossil, and biogenic
        """
        primary_share = self.vm.fuel_blend[fuel_type]["primary"]["share"]
        secondary_share = self.vm.fuel_blend[fuel_type]["secondary"]["share"]

        primary_CO2 = self.vm.fuel_blend[fuel_type]["primary"]["CO2"]
        secondary_CO2 = self.vm.fuel_blend[fuel_type]["secondary"]["CO2"]

        primary_biogenic_share = self.vm.fuel_blend[fuel_type]["primary"][
            "biogenic share"
        ]
        secondary_biogenic_share = self.vm.fuel_blend[fuel_type]["secondary"][
            "biogenic share"
        ]

        return (
            primary_share * primary_CO2 * (1 - primary_biogenic_share)
            + secondary_share * secondary_CO2 * (1 - secondary_biogenic_share),
            primary_share * primary_CO2 * primary_biogenic_share
            + secondary_share * secondary_CO2 * secondary_biogenic_share,
        )

    def fill_in_A_matrix(self):
        """
        Fill-in the A matrix. Does not return anything. Modifies in place.
        Shape of the A matrix (values, products, activities).

        :param array: :attr:`array` from :class:`CarModel` class
        """

        pass

    def add_fuel_cell_stack(self):

        self.A[
            :,
            self.find_input_indices(("Ancillary BoP",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = (
            self.array[self.array_inputs["fuel cell ancillary BoP mass"], :] * -1
        )

        self.A[
            :,
            self.find_input_indices(("Essential BoP",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = (
            self.array[self.array_inputs["fuel cell essential BoP mass"], :] * -1
        )

        # note: `Stack`refers to the power of the stack, not mass
        self.A[
            :,
            self.find_input_indices(contains=("Stack",), excludes=("PEM",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = (
            self.array[self.array_inputs["fuel cell power"], :]
            * (1 + self.array[self.array_inputs["fuel cell lifetime replacements"], :])
            * -1
        )

    def add_hydrogen_tank(self):

        hydro_tank_type = self.vm.energy_storage.get(
            "hydrogen", {"tank type": "carbon fiber"}
        )["tank type"]

        dict_tank_map = {
            "carbon fiber": "Fuel tank, compressed hydrogen gas, 700bar, with carbon fiber",
            "hdpe": "Fuel tank, compressed hydrogen gas, 700bar, with HDPE liner",
            "aluminium": "Fuel tank, compressed hydrogen gas, 700bar, with aluminium liner",
        }

        index = self.get_index_vehicle_from_array(["FCEV"])

        self.A[
            :,
            self.find_input_indices((dict_tank_map[hydro_tank_type],)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ", "FCEV")),
        ] = (
            self.array[self.array_inputs["fuel tank mass"], :, index] * -1
        )

    def add_battery(self):

        # Start of printout
        print(
            "****************** IMPORTANT BACKGROUND PARAMETERS ******************",
            end="\n * ",
        )

        # Energy storage
        print(f"The country of use is {self.vm.country}.", end="\n * ")

        battery_tech = list(set(list(self.vm.energy_storage["electric"].values())))
        battery_origin = self.vm.energy_storage.get("origin", "CN")

        print(
            "Power and energy batteries produced "
            f"in {battery_origin} using {battery_tech} chemistry/ies",
            end="\n * ",
        )

        # Use the NMC inventory of Schmidt et al. 2019
        self.A[
            :,
            self.find_input_indices(("Battery BoP",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = (
            self.array[self.array_inputs["battery BoP mass"], :]
            * (1 + self.array[self.array_inputs["battery lifetime replacements"], :])
        ) * -1

        self.A[
            :,
            self.find_input_indices((f"Battery cell, {battery_tech[0]}",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = (
            self.array[self.array_inputs["battery cell mass"], :]
            * (1 + self.array[self.array_inputs["battery lifetime replacements"], :])
        ) * -1

        # Set an input of electricity, given the country of manufacture
        print(f"Battery cell, {battery_tech[0]}")

        electricity_batt = self.find_input_requirement(
            value_in="kilowatt hour",
            value_out=f"Battery cell, {battery_tech[0]}",
            find_input_by="unit",
            filter_activities=["NMC", "LFP", "LTO", "NCA"],
            zero_out_input=True,
        )

        for year in self.scope["year"]:
            index = self.get_index_vehicle_from_array(year)

            self.A[
                np.ix_(
                    np.arange(self.iterations),
                    self.find_input_indices(
                        (f"electricity supply for battery production, {year}",)
                    ),
                    self.find_input_indices(
                        (f"{self.vm.vehicle_type.capitalize()}, ", str(year))
                    ),
                )
            ] = (
                electricity_batt
                * (
                    (
                        self.array[self.array_inputs["battery cell mass"], :, index]
                        * (
                            1
                            + self.array[
                                self.array_inputs["battery lifetime replacements"],
                                :,
                                index,
                            ]
                        )
                    )
                    * -1
                ).values[:, np.newaxis, :]
            )

        # Battery EoL
        self.A[
            :,
            self.find_input_indices(("market for used Li-ion battery",)),
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
        ] = self.array[
            [self.array_inputs[l] for l in ["battery cell mass", "battery BoP mass"]], :
        ].sum(
            axis=0
        ) * (
            1 + self.array[self.array_inputs["battery lifetime replacements"], :]
        )

    def add_cng_tank(self):

        index = self.get_index_vehicle_from_array("ICEV-g")
        self.A[
            :,
            self.find_input_indices(
                contains=("Fuel tank, compressed natural gas, 200 bar",)
            ),
            self.find_input_indices(
                contains=(f"{self.vm.vehicle_type.capitalize()}, ", "ICEV-g")
            ),
        ] = (
            self.array[self.array_inputs["fuel tank mass"], :, index] * -1
        )

    def add_vehicle_to_transport_dataset(self):

        self.A[
            :,
            self.find_input_indices((f"{self.vm.vehicle_type.capitalize()}, ",)),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = (
            -1 / self.array[self.array_inputs["lifetime kilometers"]]
        )

    def display_renewable_rate_in_mix(self):

        sum_renew = self.define_renewable_rate_in_mix()

        for y, year in enumerate(self.scope["year"]):
            if y + 1 == len(self.scope["year"]):
                end_str = "\n * "
            else:
                end_str = "\n \t * "

            print(
                f"in {year}, % of renewable: {np.round(sum_renew[y] * 100)}.",
                end=end_str,
            )

    def add_electricity_to_electric_vehicles(self) -> None:

        electric_powertrains = [
            "BEV",
            "BEV-opp",
            "BEV-motion",
            "BEV-depot",
            "PHEV-p",
            "PHEV-d",
        ]

        if any(True for x in electric_powertrains if x in self.scope["powertrain"]):
            for year in self.scope["year"]:
                index = self.get_index_vehicle_from_array(
                    electric_powertrains, year, method="and"
                )

                self.A[
                    np.ix_(
                        np.arange(self.iterations),
                        self.find_input_indices(
                            (f"electricity supply for electric vehicles, {year}",)
                        ),
                        self.find_input_indices(
                            contains=(
                                f"transport, {self.vm.vehicle_type}, ",
                                str(year),
                            ),
                            excludes=("ICEV", "FCEV", ", HEV", "Human"),
                        ),
                    )
                ] = (
                    self.array[self.array_inputs["electricity consumption"], :, index]
                    * -1
                ).values[
                    :, None, :
                ]

    def add_hydrogen_to_fuel_cell_vehicles(self) -> None:

        if "FCEV" in self.scope["powertrain"]:

            index = self.get_index_vehicle_from_array("FCEV")

            print(
                "{} is completed by {}.".format(
                    self.vm.fuel_blend["hydrogen"]["primary"]["type"],
                    self.vm.fuel_blend["hydrogen"]["secondary"]["type"],
                ),
                end="\n \t * ",
            )

            for y, year in enumerate(self.scope["year"]):
                if y + 1 == len(self.scope["year"]):
                    end_str = "\n * "
                else:
                    end_str = "\n \t * "

                print(
                    f"in {year} _________________________________________ "
                    f"{np.round(self.vm.fuel_blend['hydrogen']['secondary']['share'][y]* 100)}%",
                    end=end_str,
                )

                # Fuel supply
                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                self.A[
                    :,
                    self.find_input_indices(
                        ("fuel supply for hydrogen vehicles", str(year))
                    ),
                    self.find_input_indices(
                        (f"transport, {self.vm.vehicle_type}, ", "FCEV", str(year))
                    ),
                ] = (
                    self.array[self.array_inputs["fuel mass"], :, ind_array]
                    / self.array[
                        self.array_inputs[RANGE_PARAM[self.vm.vehicle_type]],
                        :,
                        ind_array,
                    ]
                    * -1
                )

    def display_fuel_blend(self, fuel) -> None:

        print(
            "{} is completed by {}.".format(
                self.vm.fuel_blend[fuel]["primary"]["type"],
                self.vm.fuel_blend[fuel]["secondary"]["type"],
            ),
            end="\n \t * ",
        )

        for y, year in enumerate(self.scope["year"]):
            if y + 1 == len(self.scope["year"]):
                end_str = "\n * "
            else:
                end_str = "\n \t * "

            print(
                f"in {year} _________________________________________ {np.round(self.vm.fuel_blend[fuel]['secondary']['share'][y] * 100)}%",
                end=end_str,
            )

    def add_carbon_dioxide_emissions(
        self, year, powertrain, powertrain_short, fossil_co2, biogenic_co2
    ) -> None:

        ind_array = [
            x
            for x in self.get_index_vehicle_from_array(year)
            if x in self.get_index_vehicle_from_array(powertrain)
        ]
        idx = [f"transport, {self.vm.vehicle_type}, ", str(year)]
        idx.extend(powertrain_short)

        self.A[
            :,
            self.inputs[("Carbon dioxide, fossil", ("air",), "kilogram")],
            self.find_input_indices(contains=tuple(idx), excludes=("battery",)),
        ] = (
            self.array[self.array_inputs["fuel mass"], :, ind_array]
            * fossil_co2
            / self.array[
                self.array_inputs[RANGE_PARAM[self.vm.vehicle_type]], :, ind_array
            ]
            * -1
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
            self.find_input_indices(contains=tuple(idx), excludes=("battery",)),
        ] = (
            ((self.array[self.array_inputs["fuel mass"], :, ind_array] * biogenic_co2))
            / self.array[
                self.array_inputs[RANGE_PARAM[self.vm.vehicle_type]], :, ind_array
            ]
            * -1
        )

    def add_sulphur_emissions(self, year, fuel, powertrain_short, powertrains) -> None:

        # Fuel-based SO2 emissions
        # Sulfur concentration value for a given country, a given year, as concentration ratio

        sulfur_concentration = self.get_sulfur_content(self.vm.country, fuel, year)
        index = self.get_index_vehicle_from_array(powertrains, [year], method="and")
        idx = [f"transport, {self.vm.vehicle_type}, ", str(year)]
        idx.extend(powertrain_short)

        if sulfur_concentration:

            self.A[
                :,
                self.inputs[("Sulfur dioxide", ("air",), "kilogram")],
                self.find_input_indices(contains=tuple(idx), excludes=("battery",)),
            ] = (
                self.array[self.array_inputs["fuel mass"], :, index]
                / self.array[
                    self.array_inputs[RANGE_PARAM[self.vm.vehicle_type]], :, index
                ]
                * -1
                * sulfur_concentration
                * (64 / 32)  # molar mass of SO2/molar mass of O2
            )

    def add_fuel_to_vehicles(self, fuel, powertrains, powertrains_short) -> None:

        if [i for i in self.scope["powertrain"] if i in powertrains]:
            index = self.get_index_vehicle_from_array(powertrains)
            (
                fuel_blend_fossil_CO2,
                fuel_blend_biogenic_CO2,
            ) = self.get_fuel_blend_carbon_intensity(fuel)

            self.display_fuel_blend(fuel)

            for y, year in enumerate(self.scope["year"]):

                ind_array = [
                    x for x in self.get_index_vehicle_from_array(year) if x in index
                ]

                # Fuel supply
                self.A[
                    :,
                    self.find_input_indices(
                        (f"fuel supply for {fuel} vehicles", str(year))
                    ),
                    self.find_input_indices(
                        contains=(
                            f"transport, {self.vm.vehicle_type}, ",
                            powertrains_short,
                            str(year),
                        ),
                        excludes=("battery",),
                    ),
                ] = (
                    self.array[self.array_inputs["fuel mass"], :, ind_array]
                    / self.array[
                        self.array_inputs[RANGE_PARAM[self.vm.vehicle_type]],
                        :,
                        ind_array,
                    ]
                ) * -1

                self.add_carbon_dioxide_emissions(
                    year,
                    powertrains,
                    [powertrains_short],
                    fuel_blend_fossil_CO2[y],
                    fuel_blend_biogenic_CO2[y],
                )

                self.add_sulphur_emissions(year, fuel, [powertrains_short], powertrains)

    def add_road_maintenance(self) -> None:

        # Infrastructure maintenance
        self.A[
            :,
            self.find_input_indices(("market for road maintenance",)),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = (
            1.29e-3 * -1
        )

    def add_road_construction(self) -> None:

        # Infrastructure
        self.A[
            :,
            self.find_input_indices(
                contains=("market for road",), excludes=("maintenance", "wear")
            ),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = (
            5.37e-7 * self.array[self.array_inputs["driving mass"], :] * -1
        )

    def add_exhaust_emissions(self) -> None:

        # Exhaust emissions
        # Non-fuel based emissions
        self.A[
            np.ix_(
                np.arange(self.iterations),
                [self.inputs[i] for i in self.exhaust_emissions],
                self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
            )
        ] = (
            self.array[[self.array_inputs[x] for x in self.exhaust_emissions.values()]]
            * -1
        ).values.transpose(
            1, 0, 2
        )

    def add_noise_emissions(self) -> None:

        # Noise emissions
        self.A[
            np.ix_(
                np.arange(self.iterations),
                [self.inputs[i] for i in self.noise_emissions],
                self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
            )
        ] = (
            self.array[[self.array_inputs[x] for x in self.noise_emissions.values()]]
            * -1
        ).values.transpose(
            1, 0, 2
        )

    def add_refrigerant_emissions(self) -> None:

        # Emissions of air conditioner refrigerant r134a
        # Leakage assumed to amount to 53g according to
        # https://treeze.ch/fileadmin/user_upload/downloads/Publications/Case_Studies/Mobility/544-LCI-Road-NonRoad-Transport-Services-v2.0.pdf
        # but only to cars with an AC system (meaning, with a cooling energy consumption)
        # and only for vehicles before 2022

        loss_rate = {
            "car": 0.75,
            "bus": 16,
            "truck": 0.94,
        }

        refill_rate = {
            "car": 0.55,
            "bus": 7.5,
            "truck": 1.1,
        }

        if any(y < 2022 for y in self.scope["year"]):
            idx_cars_before_2022 = [
                self.inputs[i]
                for i in self.inputs
                if f"transport, {self.vm.vehicle_type}, " in i[0]
                and int(re.findall(r"20(\w+)", i[0])[0]) < 22
            ]
            index = self.get_index_vehicle_from_array(
                [i for i in self.scope["year"] if i < 2022]
            )

            self.A[
                :,
                self.inputs[
                    ("Ethane, 1,1,1,2-tetrafluoro-, HFC-134a", ("air",), "kilogram")
                ],
                idx_cars_before_2022,
            ] = (
                (
                    loss_rate[self.vm.vehicle_type]
                    / self.array.values[
                        self.array_inputs["lifetime kilometers"], :, index
                    ]
                    * -1
                )
                * self.array.values[
                    self.array_inputs["cooling energy consumption"], :, index
                ]
            ).T

            self.A[
                :,
                self.find_input_indices(("market for refrigerant R134a",)),
                idx_cars_before_2022,
            ] = (
                (
                    (
                        loss_rate[self.vm.vehicle_type]
                        + refill_rate[self.vm.vehicle_type]
                    )
                    / self.array.values[
                        self.array_inputs["lifetime kilometers"], :, index
                    ]
                    * -1
                )
                * self.array.values[
                    self.array_inputs["cooling energy consumption"], :, index
                ]
            ).T

    def add_abrasion_emissions(self) -> None:

        # Non-exhaust emissions

        abrasion_datasets = {
            (
                "road wear",
                "two-wheeler",
            ): "market for road wear emissions, passenger car",
            (
                "brake wear",
                "two-wheeler",
            ): "market for brake wear emissions, passenger car",
            (
                "tire wear",
                "two-wheeler",
            ): "market for tyre wear emissions, passenger car",
            ("road wear", "car"): "market for road wear emissions, passenger car",
            ("brake wear", "car"): "market for brake wear emissions, passenger car",
            ("tire wear", "car"): "market for tyre wear emissions, passenger car",
            ("road wear", "truck"): "treatment of road wear emissions, lorry",
            ("brake wear", "truck"): "treatment of brake wear emissions, lorry",
            ("tire wear", "truck"): "treatment of tyre wear emissions, lorry",
            ("road wear", "bus"): "treatment of road wear emissions, lorry",
            ("brake wear", "bus"): "treatment of brake wear emissions, lorry",
            ("tire wear", "bus"): "treatment of tyre wear emissions, lorry",
        }
        # Road wear emissions + 33.3% of re-suspended road dust
        self.A[
            :,
            self.find_input_indices(
                (abrasion_datasets[("road wear", self.vm.vehicle_type)],)
            ),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = self.array[self.array_inputs["road wear emissions"], :] + (
            0.333 * self.array[self.array_inputs["road dust emissions"], :]
        )

        self.A[
            :,
            self.find_input_indices(
                (abrasion_datasets[("tire wear", self.vm.vehicle_type)],)
            ),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = self.array[self.array_inputs["tire wear emissions"], :] + (
            0.333 * self.array[self.array_inputs["road dust emissions"], :]
        )

        # Brake wear emissions
        # BEVs only emit 20% of what a combustion engine vehicle emit according to
        # https://link.springer.com/article/10.1007/s11367-014-0792-4

        self.A[
            :,
            self.find_input_indices(
                (abrasion_datasets[("brake wear", self.vm.vehicle_type)],)
            ),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] = self.array[self.array_inputs["brake wear emissions"], :] + (
            0.333 * self.array[self.array_inputs["road dust emissions"], :]
        )

    def is_the_vehicle_compliant(self, index):
        """
        Checks if the vehicle has a positive `TtW energy` value
        """

        name = self.rev_inputs[index][0]

        if (
            any([w for w in self.scope["powertrain"] if w in name])
            and any([w for w in self.scope["year"] if str(w) in name])
            and any([w for w in self.scope["size"] if w in name])
        ):

            powertrain = [w for w in self.scope["powertrain"] if w in name][0]
            size = [w for w in self.scope["size"] if w in name][0]
            year = [w for w in self.scope["year"] if str(w) in name][0]

            return np.all(
                self.vm.array.loc[
                    dict(
                        size=size,
                        powertrain=powertrain,
                        year=year,
                        parameter="TtW energy",
                    )
                ]
                > 0
            )

        return True

    def remove_uncompliant_vehicles(self):
        """
        Remove vehicles from self.A that do not have a TtW energy superior to 0.
        """
        # Get the indices of the vehicles that are not compliant
        idx = self.find_input_indices((self.vm.vehicle_type,))
        idx.extend(
            self.find_input_indices(
                (self.vm.vehicle_type.capitalize()),
            )
        )

        idx = [i for i in idx if not self.is_the_vehicle_compliant(i)]

        # Zero-out the rows of the non-compliant vehicles
        self.A[:, :, idx] = 0
        self.A[:, idx, idx] = 1

    def change_functional_unit(self) -> None:

        load_factor = self.get_load_factor()
        idx_cars = self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",))
        idx_others = [i for i in range(self.A.shape[1]) if i not in idx_cars]

        self.A[np.ix_(np.arange(self.iterations), idx_others, idx_cars,)] *= (
            1 / load_factor
        ).reshape(-1, 1, len(idx_cars))

        # iterate through self.inputs and change the unit
        keys_to_modify = {
            key: value
            for key, value in self.inputs.items()
            if key[0].startswith(f"transport, {self.vm.vehicle_type}")
        }

        for key, value in keys_to_modify.items():
            new_key = list(key)
            new_key[2] = self.func_unit
            del self.inputs[key]
            self.inputs[tuple(new_key)] = value

        # update self.rev_inputs
        self.rev_inputs = {v: k for k, v in self.inputs.items()}

    def export_lci(
        self,
        ecoinvent_version="3.8",
        filename=f"carculator_lci",
        directory=None,
        software="brightway2",
        format="bw2io",
    ):
        """
        Export the inventory. Can export to Simapro (as csv), or brightway2 (as bw2io object, file or string).
        :param db_name:
        :param ecoinvent_version: str. "3.5", "3.6", "3.7" or "3.8"
        ::return: inventory, or the filepath where the file is saved..
        :rtype: list
        """

        if self.func_unit != "vkm":
            self.change_functional_unit()

        lci = ExportInventory(
            array=self.A,
            vehicle_model=self.vm,
            indices=self.rev_inputs,
            db_name=f"{filename}_{self.vm.vehicle_type}",
        )

        if software == "brightway2":
            return lci.write_bw2_lci(
                ecoinvent_version=ecoinvent_version,
                directory=directory,
                filename=filename,
                export_format=format,
            )

        else:
            return lci.write_simapro_lci(
                ecoinvent_version=ecoinvent_version,
                directory=directory,
                filename=filename,
            )
