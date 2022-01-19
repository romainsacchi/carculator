"""
export.py contains the class Export, which offers methods to export the inventories
in different formats.
"""

import csv
import datetime
import io
import json
import os
import uuid
from typing import Dict, List, Tuple, Union

import bw2io
import numpy as np
import pyprind
import xarray as xr
from bw2io.export.excel import create_valid_worksheet_name, safe_filename, xlsxwriter

from . import DATA_DIR, __version__


def load_mapping_37_to_36() -> Dict[
    Tuple[str, str, str, str], Tuple[str, str, str, str]
]:
    """Load mapping dictionary between ecoinvent 3.7 and 3.6"""

    # Load the matching dictionary
    filename = "ei37_to_ei36.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of activities flows match between ecoinvent 3.7 and 3.6 could not be found."
        )
    with open(filepath, encoding="utf-8") as f:
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
    (_, _, *header), *data = csv_list

    dict_ei36 = {}
    for row in data:
        (
            name,
            location,
            unit,
            ref_prod,
            name_36,
            location_36,
            unit_36,
            ref_prod_36,
        ) = row
        dict_ei36[(name, location, unit, ref_prod)] = (
            name_36,
            location_36,
            unit_36,
            ref_prod_36,
        )

    return dict_ei36


def load_references() -> Dict[str, Dict[str, str]]:
    """Load a dictionary with references/sources of datasets"""

    # Load the matching dictionary
    filename = "references.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError("The dictionary of references could not be found.")
    with open(filepath, encoding="latin1") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]
    _, *data = csv_list

    dict_reference = {}
    for row in data:
        name, source, description, special_remark, category_1, category_2 = row
        dict_reference[name] = {
            "source": source,
            "description": description,
            "special remark": special_remark,
            "category 1": category_1,
            "category 2": category_2,
        }

    return dict_reference


def load_uvek_transport_distances() -> Dict[str, Dict[str, float]]:
    """Load a dictionary with transport distances for inventory export to UVEK database"""

    # Load the matching dictionary
    filename = "transport_distance_uvek.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary with transport distances could not be found."
        )
    with open(filepath, encoding="latin1") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]
    _, *data = csv_list

    dict_distance = {}
    for row in data:
        name, _, _, train_RER, truck_RER, barge_RER, train_CH, truck_CH, barge_CH = row

        dict_distance[name] = {
            "train RER": float(train_RER),
            "truck RER": float(truck_RER),
            "barge RER": float(barge_RER),
            "train CH": float(train_CH),
            "truck CH": float(truck_CH),
            "barge CH": float(barge_CH),
        }

    return dict_distance


def load_mapping_37_to_35() -> Dict[
    Tuple[str, str, str, str], Tuple[str, str, str, str]
]:
    """Load mapping dictionary between ecoinvent 3.7 and 3.5"""

    # Load the matching dictionary
    filename = "ei37_to_ei35.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of activities flows match between ecoinvent 3.7 and 3.5 could not be found."
        )
    with open(filepath, encoding="utf-8") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]
    _, _, *data = csv_list

    dict_ei35 = {}
    for row in data:
        (
            name,
            location,
            unit,
            ref_prod,
            name_35,
            location_35,
            unit_35,
            ref_prod_35,
        ) = row
        dict_ei35[(name, location, unit, ref_prod)] = (
            name_35,
            location_35,
            unit_35,
            ref_prod_35,
        )

    return dict_ei35


def get_simapro_biosphere() -> Dict[str, str]:

    # Load the matching dictionary between ecoinvent and Simapro biosphere flows
    # for each ecoinvent biosphere flow name, it gives the corresponding Simapro name

    filename = "simapro-biosphere.json"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of biosphere flow match between ecoinvent and Simapro could not be found."
        )
    with open(filepath, encoding="utf-8") as json_file:
        data = json.load(json_file)
    dict_bio = {}
    for d in data:
        dict_bio[d[2]] = d[1]

    return dict_bio


def get_simapro_technosphere() -> Dict[Tuple[str, str], str]:

    # Load the matching dictionary between ecoinvent and Simapro product flows

    filename = "simapro-technosphere-3.5.csv"
    filepath = DATA_DIR / filename
    with open(filepath, encoding="utf-8") as f:
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
    (_, _, *header), *data = csv_list

    dict_tech = {}
    for row in data:
        name, location, simapro_name = row
        simapro_name = simapro_name.split("|")[:2]
        dict_tech[(name, location)] = "|".join(simapro_name)

    return dict_tech


class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices, db_name="carculator export"):
        self.array: xr.DataArray = array
        self.indices: Dict[int, Tuple[str, str, str, str]] = indices

        self.rename_vehicles()

        self.db_name: str = db_name
        self.references = load_references()

        # See https://docs.brightwaylca.org/intro.html#uncertainty-type
        self.uncertainty_id = {
            # scipy.stats distr. params --> stats.array distr. params
            "triang": 5,  # c --> loc (mode), loc --> min, loc + scale --> max
            "weibull_min": 8,  # c --> shape, scale --> scale
            "gamma": 9,  # a --> shape, scale --> scale, loc --> loc
            "beta": 10,  # a --> loc, b --> shape, scale --> scale
            "lognorm": 2,  # s --> scale (std), scale --> loc (exp(mean))
            "norm": 3,  # loc --> loc (mean), scale --> scale (std)
            "uniform": 4,  # loc --> min, loc + scale --> max
            "t": 12,  # df --> shape, loc --> loc, scale --> scale
        }

        self.map_37_to_36 = load_mapping_37_to_36()
        self.map_37_to_35 = load_mapping_37_to_35()
        self.map_36_to_uvek = self.load_mapping_36_to_uvek()
        self.map_36_to_uvek_for_simapro = self.load_mapping_36_to_uvek_for_simapro()
        self.tags = self.load_tags()
        self.uvek_dist = load_uvek_transport_distances()

    def rename_vehicles(self) -> None:

        """
        Rename powertrain acronyms to full length descriptive terms

        """

        d_names = {
            "ICEV-d": "diesel",
            "ICEV-p": "gasoline",
            "ICEV-g": "compressed gas",
            "HEV-p": "gasoline hybrid",
            "HEV-d": "diesel hybrid",
            "PHEV-p": "plugin gasoline hybrid",
            "PHEV-d": "plugin diesel hybrid",
            "BEV": "battery electric",
            "FCEV": "fuel cell electric",
        }

        for k, value in self.indices.items():
            for key, val in d_names.items():
                if key in value[0]:
                    new_val = list(value)
                    new_val[0] = new_val[0].replace(key, val)
                    self.indices[k] = tuple(new_val)

    @staticmethod
    def load_tags() -> Dict[str, str]:
        """Loads dictionary of tags for further use in BW2"""

        filename = "tags.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError("The dictionary of tags could not be found.")
        with open(filepath, encoding="utf-8") as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        data = csv_list

        dict_tags = {}
        for row in data:
            name, tag = row
            dict_tags[name] = tag

        return dict_tags

    @staticmethod
    def load_mapping_36_to_uvek() -> Dict[
        Tuple[str, str, str, str], Tuple[str, str, str, str]
    ]:
        """Load mapping dictionary between ecoinvent 3.6 and UVEK"""

        # Load the matching dictionary between ecoinvent and Simapro biosphere flows
        filename = "uvek_mapping.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of activities flows match between ecoinvent 3.6 and UVEK could not be found."
            )
        with open(filepath, encoding="utf-8") as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        (_, _, *header), *data = csv_list

        dict_uvek = {}
        for row in data:
            (
                name,
                location,
                unit,
                ref_prod,
                uvek_name,
                uvek_loc,
                uvek_unit,
                uvek_ref_prod,
                simapro_name,
            ) = row
            dict_uvek[(name, ref_prod, unit, location)] = (
                uvek_name,
                uvek_ref_prod,
                uvek_unit,
                uvek_loc,
            )
        return dict_uvek

    @staticmethod
    def load_mapping_36_to_uvek_for_simapro() -> Dict[Tuple[str, str, str, str], str]:
        """Load mapping dictionary between ecoinvent 3.6 and UVEK for Simapro name format"""

        # Load the matching dictionary between ecoinvent and Simapro biosphere flows
        filename = "uvek_mapping.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of activities flows match between ecoinvent 3.6 and UVEK could not be found."
            )
        with open(filepath, encoding="utf-8") as f:
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        (_, _, *header), *data = csv_list

        dict_uvek = {}
        for row in data:
            (
                name,
                location,
                unit,
                ref_prod,
                _,
                _,
                _,
                _,
                simapro_name,
            ) = row
            dict_uvek[(name, location, unit, ref_prod)] = simapro_name

        return dict_uvek

    def write_lci(
        self,
        presamples: bool,
        ecoinvent_version: str,
        vehicle_specs: xr.DataArray,
    ) -> List[Dict]:
        """
        Return the inventory as a dictionary
        If if there several values for one exchange, uncertainty information is generated.
        If `presamples` is True, returns the inventory as well as a `presamples` matrix.
        If `presamples` is False, returns the inventory with characterized uncertainty information.

        :returns: a dictionary that contains all the exchanges
        :rtype: dict
        """

        uvek_activities_to_remove = [
            "market for activated carbon, granular",
            "market for iodine",
            "market for molybdenum trioxide",
            "market for fly ash and scrubber sludge",
        ]

        ei35_activities_to_remove = ["latex production"]

        list_act = []

        if presamples:
            presamples_matrix = []

        # List of coordinates for non-zero values
        non_zeroes = np.nonzero(self.array[0, :, :])
        # List of coordinates where activities present more than once (to filter out "empty" activities, that is,
        # activities with only one reference product exchange)
        u, c = np.unique(non_zeroes[1], return_counts=True)
        dup = u[c > 1]

        # Filter out coordinates of "empty" activities
        coords = np.column_stack(
            (
                non_zeroes[0][np.isin(non_zeroes[1], dup)],
                non_zeroes[1][np.isin(non_zeroes[1], dup)],
            )
        )

        # Iterate through activities
        bar = pyprind.ProgBar(len(dup))
        for d in dup:
            bar.update(item_id=d)
            list_exc = []
            for row, col in coords[coords[:, 1] == d]:
                tuple_output = self.indices[col]
                tuple_input = self.indices[row]
                mult_factor = 1

                if ecoinvent_version == "3.6":
                    tuple_output = self.map_37_to_36.get(tuple_output, tuple_output)
                    tuple_input = self.map_37_to_36.get(tuple_input, tuple_input)

                if ecoinvent_version == "3.5":
                    tuple_output = self.map_37_to_35.get(tuple_output, tuple_output)
                    tuple_input = self.map_37_to_35.get(tuple_input, tuple_input)

                    if tuple_output[0] in ei35_activities_to_remove:
                        continue

                    if tuple_input[0] in ei35_activities_to_remove:
                        continue

                if ecoinvent_version == "uvek":

                    tuple_output = self.map_36_to_uvek.get(tuple_output, tuple_output)

                    if tuple_input[0] in uvek_activities_to_remove:
                        continue

                    tuple_input = self.map_36_to_uvek.get(tuple_input, tuple_input)

                if len(self.array[:, row, col]) == 1:
                    # No uncertainty, only one value
                    amount = self.array[0, row, col] * mult_factor
                    uncertainty = [("uncertainty type", 1)]

                elif np.all(
                    np.isclose(self.array[:, row, col], self.array[0, row, col])
                ):
                    # Several values, but all the same, so no uncertainty
                    amount = self.array[0, row, col] * mult_factor
                    uncertainty = [("uncertainty type", 1)]

                else:
                    # Uncertainty

                    # Generate pre-sampled values
                    amount = np.median(self.array[:, row, col]) * mult_factor
                    uncertainty = [("uncertainty type", 1)]
                    if len(tuple_input) > 3:
                        type_exc = "technosphere"
                    else:
                        type_exc = "biosphere"

                    presamples_matrix.append(
                        (
                            self.array[:, row, col] * -1,
                            [(tuple_input, tuple_output, type_exc)],
                            type_exc,
                        )
                    )

                # Look for a tag, if any
                tag = [
                    self.tags[t] for t in list(self.tags.keys()) if t in tuple_input[0]
                ]
                if len(tag) > 0:
                    tag = tag[0]
                else:
                    tag = "other"

                # If reference product
                if tuple_output == tuple_input:

                    list_exc.append(
                        {
                            "name": tuple_output[0],
                            "database": self.db_name,
                            "amount": amount,
                            "unit": tuple_output[2],
                            "type": "production",
                            "location": tuple_output[1],
                            "reference product": tuple_output[3],
                            "tag": tag,
                        }
                    )
                    list_exc[-1].update(uncertainty)

                # If not, if input is technosphere exchange
                elif len(tuple_input) > 3:
                    list_exc.append(
                        {
                            "name": tuple_input[0],
                            "database": self.db_name,
                            "amount": amount * -1,
                            "unit": tuple_input[2],
                            "type": "technosphere",
                            "location": tuple_input[1],
                            "reference product": tuple_input[3],
                            "tag": tag,
                        }
                    )
                    list_exc[-1].update(uncertainty)

                # If not, then input is biosphere exchange
                else:
                    list_exc.append(
                        {
                            "name": tuple_input[0],
                            "database": "biosphere3",
                            "amount": amount * -1,
                            "unit": tuple_input[2],
                            "type": "biosphere",
                            "categories": tuple_input[1],
                            "tag": tag,
                        }
                    )
                    list_exc[-1].update(uncertainty)

            # Look for a tag, if any
            tag = [self.tags[t] for t in list(self.tags.keys()) if t in tuple_output[0]]
            if len(tag) > 0:
                tag = tag[0]
            else:
                tag = "other"

            if tuple_output[0] in self.references:
                source = self.references[tuple_output[0]]["source"]
                description = self.references[tuple_output[0]]["description"]
                special_remark = self.references[tuple_output[0]]["special remark"]
            else:

                try:
                    key = [
                        k
                        for k in self.references
                        if k.lower() in tuple_output[0].lower()
                    ][0]
                except IndexError:
                    print(tuple_output[0])
                source = self.references[key]["source"]
                description = self.references[key]["description"]
                special_remark = self.references[key]["special remark"]

            string = ""
            if "passenger car" in tuple_output[0].lower():

                d_pwt = {
                    "gasoline": "ICEV-p",
                    "diesel": "ICEV-d",
                    "compressed gas": "ICEV-g",
                    "diesel hybrid": "HEV-d",
                    "gasoline hybrid": "HEV-p",
                    "plugin diesel hybrid": "PHEV-d",
                    "plugin gasoline hybrid": "PHEV-p",
                    "battery electric": "BEV",
                    "fuel cell electric": "FCEV",
                }

                d_units = {
                    "lifetime kilometers": "[km]",
                    "kilometers per year": "[km/year]",
                    "range": "[km]",
                    "TtW efficiency": "[%]",
                    "TtW energy": "[kj/km]",
                    "electric energy stored": "[kWh]",
                    "oxidation energy stored": "[kWh]",
                    "combustion power share": "[%]",
                    "combustion power": "[kW]",
                    "electric power": "[kW]",
                    "curb mass": "[kg]",
                    "driving mass": "[kg]",
                    "energy battery mass": "[kg]",
                    "fuel cell system efficiency": "[%]",
                }

                d_names = {
                    "lifetime kilometers": "Km over lifetime",
                    "kilometers per year": "Yearly mileage",
                    "range": "Autonomy on a full tank/battery",
                    "TtW efficiency": "Tank-to-wheel efficiency",
                    "TtW energy": "Tank-to-wheel energy consumption",
                    "electric energy stored": "Battery capacity",
                    "oxidation energy stored": "Fuel tank capacity",
                    "combustion power share": "Power share from combustion engine",
                    "combustion power": "Combustion engine power",
                    "electric power": "Electric motor power",
                    "curb mass": "Curb mass (excl. driver and cargo)",
                    "driving mass": "Driving mass (incl. driver and cargo)",
                    "energy battery mass": "Mass of battery",
                    "fuel cell system efficiency": "Fuel cell system efficiency",
                }

                split_name = [t.strip() for t in tuple_output[0].split(",")]

                pwt, size, year = (None, None, None)

                if "fleet average" not in tuple_output[0]:
                    if "transport, " in tuple_output[0]:
                        if "battery electric" in split_name[2]:
                            _, _, pwt, _, size, year = split_name
                        elif "fuel cell" in split_name[2]:
                            _, _, pwt, size, year = split_name
                        else:
                            _, _, pwt, size, year, _ = split_name
                    else:
                        if "battery electric" in split_name[1]:
                            _, pwt, _, size, year = split_name
                        elif "fuel cell" in split_name[1]:
                            _, pwt, size, year = split_name
                        else:
                            _, pwt, size, year, _ = split_name
                else:
                    if split_name[3] == "fleet average":

                        if len(split_name) == 6:

                            if split_name[4] in [
                                "Micro",
                                "Mini",
                                "Small",
                                "Lower medium",
                                "Medium",
                                "Medium SUV",
                                "Large",
                                "Large SUV",
                                "Van",
                            ]:
                                _, _, _, _, size, year = [
                                    t.strip() for t in tuple_output[0].split(",")
                                ]
                                pwt = None
                            else:
                                _, _, _, _, pwt, year = [
                                    t.strip() for t in tuple_output[0].split(",")
                                ]
                                size = None

                        if len(split_name) == 5:
                            _, _, _, _, year = [
                                t.strip() for t in tuple_output[0].split(",")
                            ]
                            size = None
                            pwt = None

                    else:
                        if split_name[1] == "battery electric":

                            _, pwt, _, size, year, _ = [
                                t.strip() for t in tuple_output[0].split(",")
                            ]
                        elif split_name[1] == "fuel cell electric":
                            _, pwt, size, year, _ = [
                                t.strip() for t in tuple_output[0].split(",")
                            ]
                        else:
                            _, pwt, size, year, _, _ = [
                                t.strip() for t in tuple_output[0].split(",")
                            ]
                if not size is None and pwt is not None:
                    pwt = d_pwt[pwt]

                    if not vehicle_specs is None:

                        for p in vehicle_specs.parameter.values:

                            try:
                                val = vehicle_specs.sel(
                                    powertrain=pwt,
                                    size=size,
                                    year=int(year),
                                    value=0,
                                    parameter=p,
                                ).values

                                if val != 0:

                                    if p in (
                                        "TtW efficiency",
                                        "combustion power share",
                                        "capacity utilization",
                                        "fuel cell system efficiency",
                                    ):
                                        val = int(val * 100)
                                    else:
                                        val = int(val)

                                    string += f"{d_names[p]}: {val} {d_units[p]}. "
                            except KeyError:
                                print(
                                    f"Could not find vehicle specs for {pwt} {size} {year}"
                                )
                else:
                    if not pwt is None:
                        pwt = d_pwt[pwt]
                        string = f"Fleet average {pwt} vehicle in {year}, all sizes considered."

                    else:
                        string = "Fleet average vehicle, all sizes and powertrains considered."

            # Added transport distances if the inventory
            # is meant for the UVEK database
            if ecoinvent_version == "uvek":
                dist_train, dist_truck, dist_barge = (0, 0, 0)
                if tuple_output[1] in (
                    "RER",
                    "Europe without Switzerland",
                    "SE",
                    "GLO",
                    "DE",
                    "JP",
                    "CN",
                ):
                    for exc in list_exc:
                        if exc["name"] in self.uvek_dist:
                            dist_train += (
                                self.uvek_dist[exc["name"]]["train RER"]
                                * float(exc["amount"])
                                / 1000
                            )
                            dist_truck += (
                                self.uvek_dist[exc["name"]]["truck RER"]
                                * float(exc["amount"])
                                / 1000
                            )
                            dist_barge += (
                                self.uvek_dist[exc["name"]]["barge RER"]
                                * float(exc["amount"])
                                / 1000
                            )

                    if dist_train > 0:

                        list_exc.append(
                            {
                                "name": "market for transport, freight train",
                                "database": "ecoinvent",
                                "amount": dist_train,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "Europe without Switzerland",
                                "reference product": "transport, freight train",
                            }
                        )
                    if dist_truck > 0:

                        list_exc.append(
                            {
                                "name": "market for transport, freight, lorry >32 metric ton, EURO4",
                                "database": "ecoinvent",
                                "amount": dist_truck,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "RER",
                                "reference product": "transport, freight, lorry >32 metric ton, EURO4",
                            }
                        )
                    if dist_barge > 0:

                        list_exc.append(
                            {
                                "name": "market for transport, freight, inland waterways, barge",
                                "database": "ecoinvent",
                                "amount": dist_barge,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "RER",
                                "reference product": "transport, freight, inland waterways, barge",
                            }
                        )

                elif tuple_output[1] == "CH":

                    for exc in list_exc:
                        if exc["name"] in self.uvek_dist:
                            dist_train += (
                                self.uvek_dist[exc["name"]]["train CH"]
                                * float(exc["amount"])
                                / 1000
                            )
                            dist_truck += (
                                self.uvek_dist[exc["name"]]["truck CH"]
                                * float(exc["amount"])
                                / 1000
                            )
                            dist_barge += (
                                self.uvek_dist[exc["name"]]["barge CH"]
                                * float(exc["amount"])
                                / 1000
                            )

                    if dist_train > 0:
                        list_exc.append(
                            {
                                "name": "market for transport, freight train",
                                "database": "ecoinvent",
                                "amount": dist_train,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "CH",
                                "reference product": "transport, freight train",
                            }
                        )
                    if dist_truck > 0:
                        list_exc.append(
                            {
                                "name": "market for transport, freight, lorry >32 metric ton, EURO4",
                                "database": "ecoinvent",
                                "amount": dist_truck,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "RER",
                                "reference product": "transport, freight, lorry >32 metric ton, EURO4",
                            }
                        )
                    if dist_barge > 0:
                        list_exc.append(
                            {
                                "name": "market for transport, freight, inland waterways, barge",
                                "database": "ecoinvent",
                                "amount": dist_barge,
                                "unit": "ton kilometer",
                                "type": "technosphere",
                                "location": "RER",
                                "reference product": "transport, freight, inland waterways, barge",
                            }
                        )

            list_act.append(
                {
                    "production amount": 1,
                    "database": self.db_name,
                    "name": tuple_output[0],
                    "unit": tuple_output[2],
                    "location": tuple_output[1],
                    "exchanges": list_exc,
                    "reference product": tuple_output[3],
                    "type": "process",
                    "code": str(uuid.uuid1()),
                    "tag": tag,
                    "source": source,
                    "description": description,
                    "special remark": special_remark,
                    "comment": string,
                }
            )
        if presamples:
            return list_act, presamples_matrix

        return list_act

    def write_lci_to_excel(
        self,
        ecoinvent_version: str,
        software_compatibility: str,
        vehicle_specs: xr.DataArray = None,
        directory: str = None,
        filename: str = None,
        export_format: str = "file",
    ) -> Union[bytes, str]:
        """
        Export a file that can be consumed by the software defined in `software_compatibility`.
        Alternatively, exports a string representation of the file (in case the invenotry should be downloaded
        from a browser, for example)

        :param vehicle_specs:
        :param filename:
        :param export_format:
        :param directory: str. path to export the file to.
        :type directory: str or pathlib.Path
        :param ecoinvent_version: str. "3.5", "3.6" or "uvek"
        :type ecoinvent_version: str
        :param software_compatibility: str. "brightway2" or "simapro"
        :type software_compatibility: str
        If "string", returns a string.
        :returns: returns the file path of the exported inventory.
        :rtype: str.
        """

        if software_compatibility == "brightway2":
            if filename is None:
                safe_name = (
                    safe_filename(
                        f"carculator_inventory_export_{datetime.date.today()}_brightway2",
                        False,
                    )
                    + ".xlsx"
                )
            else:
                safe_name = (
                    safe_filename(
                        filename,
                        False,
                    )
                    + ".xlsx"
                )
        else:
            safe_name = (
                safe_filename(
                    f"carculator_inventory_export_{datetime.date.today()}_simapro",
                    False,
                )
                + ".csv"
            )

        if directory is None:
            filepath_export = safe_name
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)
            filepath_export = os.path.join(directory, safe_name)

        list_act = self.write_lci(
            presamples=False,
            ecoinvent_version=ecoinvent_version,
            vehicle_specs=vehicle_specs,
        )

        if software_compatibility == "brightway2":
            data = self.format_data_for_lci_for_bw2(list_act)

            if export_format == "file":

                workbook = xlsxwriter.Workbook(filepath_export)
                bold = workbook.add_format({"bold": True})
                bold.set_font_size(12)
                highlighted = {
                    "Activity",
                    "Database",
                    "Exchanges",
                    "Parameters",
                    "Database parameters",
                    "Project parameters",
                }
                frmt = lambda x: bold if row[0] in highlighted else None
                sheet = workbook.add_worksheet(create_valid_worksheet_name("test"))

                for row_index, row in enumerate(data):
                    for col_index, value in enumerate(row):
                        if value is None:
                            continue
                        elif isinstance(value, float):
                            sheet.write_number(row_index, col_index, value, frmt(value))
                        else:
                            sheet.write_string(row_index, col_index, value, frmt(value))
                print(f"Inventories exported to {filepath_export}.")
                workbook.close()

            if export_format == "string":
                output = io.BytesIO()
                workbook = xlsxwriter.Workbook(output, {"in_memory": True})
                bold = workbook.add_format({"bold": True})
                bold.set_font_size(12)
                highlighted = {
                    "Activity",
                    "Database",
                    "Exchanges",
                    "Parameters",
                    "Database parameters",
                    "Project parameters",
                }

                frmt = lambda x: bold if row[0] in highlighted else None
                sheet = workbook.add_worksheet("carculator export")

                for row_index, row in enumerate(data):
                    for col_index, value in enumerate(row):
                        if value is None:
                            continue
                        elif isinstance(value, float):
                            sheet.write_number(row_index, col_index, value, frmt(value))
                        else:
                            sheet.write_string(row_index, col_index, value, frmt(value))

                workbook.close()
                output.seek(0)
                return output.read()

        else:
            if export_format == "file":
                with open(
                    filepath_export, "w", newline="", encoding="latin1"
                ) as csvFile:
                    writer = csv.writer(csvFile, delimiter=";")
                    rows = self.format_data_for_lci_for_simapro(
                        data=list_act, ei_version=ecoinvent_version
                    )
                    for row in rows:
                        writer.writerow(row)
                csvFile.close()
                print(f"Inventories exported to {filepath_export}.")

            if export_format == "string":
                csvFile = io.StringIO()
                writer = csv.writer(
                    csvFile,
                    delimiter=";",
                    quoting=csv.QUOTE_NONE,
                    quotechar="",
                    escapechar="\\",
                )
                rows = self.format_data_for_lci_for_simapro(list_act, ecoinvent_version)
                for row in rows:
                    writer.writerow(row)
                csvFile.seek(0)
                return csvFile.read()

    def format_data_for_lci_for_bw2(self, data):

        rows = []
        rows.extend((["Database", self.db_name], ("format", "Excel spreadsheet")))
        rows.append([])

        for k in data:
            if k.get("exchanges"):
                rows.extend(
                    (
                        ["Activity", k["name"]],
                        ("location", k["location"]),
                        ("production amount", float(k["production amount"])),
                        ("reference product", k.get("reference product")),
                        ("type", "process"),
                        ("unit", k["unit"]),
                        ("worksheet name", "None"),
                        ("source", k["source"]),
                        ("description", k["description"]),
                        ("special remark", k["special remark"]),
                        ("comment", k["comment"]),
                        ["Exchanges"],
                        [
                            "name",
                            "amount",
                            "database",
                            "location",
                            "unit",
                            "categories",
                            "type",
                            "reference product",
                            "tag",
                        ],
                    )
                )

                for e in k["exchanges"]:
                    rows.append(
                        [
                            e["name"],
                            float(e["amount"]),
                            e["database"],
                            e.get("location", "None"),
                            e["unit"],
                            "::".join(e.get("categories", ())),
                            e["type"],
                            e.get("reference product"),
                            e.get("tag", "other"),
                        ]
                    )
            else:
                rows.extend(
                    (
                        ["Activity", k["name"]],
                        ("type", "biosphere"),
                        ("unit", k["unit"]),
                        ("worksheet name", "None"),
                    )
                )
            rows.append([])

        return rows

    def format_data_for_lci_for_simapro(
        self, data: List[Dict], ei_version: str
    ) -> List[List]:

        # not all biosphere flows exist in simapro
        simapro_biosphere_flows_to_remove = [
            "Gangue, in ground",
            "Water, turbine use, unspecified natural origin",
            "Oxygen",
            "Volume occupied, reservoir",
            "Xenon-135",
            "Noble gases, radioactive, unspecified",
            "Radon-222",
            "Xenon-133",
            "Hydrogen-3, Tritium",
            "Radon-222",
            "Radon-220",
            "Oxygen",
            "Occupation, traffic area, road network",
            "Energy, gross calorific value, in biomass, primary forest",
            "Carbon-14",
        ]

        headers = [
            "{SimaPro 9.1.1.1}",
            "{processes}",
            "{Project: carculator import"
            + f"{datetime.datetime.today():%d.%m.%Y}"
            + "}",
            "{CSV Format version: 9.0.0}",
            "{CSV separator: Semicolon}",
            "{Decimal separator: .}",
            "{Date separator: .}",
            "{Short date format: dd.MM.yyyy}",
            "{Export platform IDs: No}",
            "{Skip empty fields: No}",
            "{Convert expressions to constants: No}",
            "{Selection: Selection(1)}",
            "{Related objects(system descriptions, substances, units, etc.): Yes}",
            "{Include sub product stages and processes: Yes}",
        ]
        fields = [
            "Process",
            "Category type",
            "Type",
            "Process name",
            "Time Period",
            "Geography",
            "Technology",
            "Comment",
            "Representativeness",
            "Cut off rules",
            "Capital goods",
            "Date",
            "Boundary with nature",
            "Infrastructure",
            "Record",
            "Generator",
            "Literature references",
            "External documents",
            "Collection method",
            "Data treatment",
            "Verification",
            "System description",
            "Allocation rules",
            "Products",
            "Waste treatment",
            "Materials/fuels",
            "Resources",
            "Emissions to air",
            "Emissions to water",
            "Emissions to soil",
            "Final waste flows",
            "Non material emission",
            "Social issues",
            "Economic issues",
            "Waste to treatment",
            "End",
        ]
        simapro_units = {
            "kilogram": "kg",
            "cubic meter": "m3",
            "kilowatt hour": "kWh",
            "kilometer": "km",
            "ton kilometer": "tkm",
            "megajoule": "MJ",
            "unit": "p",
            "square meter": "m2",
            "kilowatt": "p",
            "hour": "hr",
            "square meter-year": "m2a",
            "meter": "m",
            "vehicle-kilometer": "vkm",
            "person-kilometer": "personkm",
            "meter-year": "my",
        }

        dict_tech = get_simapro_technosphere()
        dict_bio = get_simapro_biosphere()

        rows = []

        for item in headers:
            rows.append([item])
        rows.append([])

        list_own_datasets = []

        for a in data:
            list_own_datasets.append(
                a["name"].capitalize() + " {" + a.get("location", "GLO") + "}"
            )

        # We loop through the activities
        for a in data:

            # We fetch teh main and sub categories (sub category is in fact a path)
            if a["name"] in self.references:
                main_category = self.references[a["name"]]["category 1"]
                category = self.references[a["name"]]["category 2"]
                source = self.references[a["name"]]["source"]
                description = self.references[a["name"]]["description"]
                special_remark = self.references[a["name"]]["special remark"]
            else:
                # if we cannot find it, it's because some keys are more general
                key = [
                    k for k in self.references.keys() if k.lower() in a["name"].lower()
                ][0]
                main_category = self.references[key]["category 1"]
                category = self.references[key]["category 2"]
                source = self.references[key]["source"]
                description = self.references[key]["description"]
                special_remark = self.references[key]["special remark"]

            # We loop through the fields SimaPro expects to see
            for item in fields:

                # If it is a waste treatment activity, we skip the field `Products`
                if main_category == "waste treatment" and item == "Products":
                    continue

                # It is not a waste treatment activity, we skip the field `Waste treatment`
                if main_category != "waste treatment" and item == "Waste treatment":
                    continue

                rows.append([item])

                if item == "Process name":

                    if ei_version in ("3.5", "3.6"):
                        name = (
                            a["name"].capitalize()
                            + " {"
                            + a.get("location", "GLO")
                            + "}"
                            + "| Cut-off, U"
                        )

                    if ei_version == "uvek":
                        name = a["name"] + "/" + a["location"] + " U"

                    rows.append([name])

                if item == "Type":
                    rows.append(["Unit process"])

                if item == "Comment":

                    if a["comment"] != "":
                        string = a["comment"]
                    else:
                        string = ""

                    string += "Originally published in: "
                    string += source

                    if description != "":
                        string += " Description: "
                        string += description

                    if special_remark != "":
                        string += " Special remark(s): "
                        string += special_remark

                    rows.append([string])

                if item == "Category type":
                    rows.append([main_category])

                if item == "Generator":
                    rows.append(["carculator " + str(__version__)])

                if item == "Geography":
                    rows.append([a["location"]])

                if item == "Time Period":
                    rows.append(
                        ["Between 2010 and 2020. Extrapolated to the selected years."]
                    )

                if item == "Date":
                    rows.append([f"{datetime.datetime.today():%d.%m.%Y}"])

                if item in (
                    "Cut off rules",
                    "Capital goods",
                    "Technology",
                    "Representativeness",
                    "Boundary with nature",
                ):
                    rows.append(["Unspecified"])

                if item == "Infrastructure":
                    rows.append(["Yes"])

                if item == "External documents":
                    rows.append(["https://carculator.psi.ch"])

                if item in "System description":
                    rows.append(["carculator"])

                if item in "Allocation rules":
                    rows.append(
                        [
                            "In the instance of joint-production, allocation of process burden based on"
                            "economic relative revenue of each co-product."
                        ]
                    )

                if item == "Literature references":
                    rows.append(["Sacchi et al. 2020"])

                if item == "Collection method":
                    rows.append(
                        [
                            "Modeling and assumptions: https://carculator.readthedocs.io/en/latest/modeling.html"
                        ]
                    )

                if item == "Verification":
                    rows.append(["In review. Susceptible to change."])

                if item == "Waste treatment":
                    if ei_version in ("3.5", "3.6"):
                        rows.append(
                            [
                                dict_tech.get((a["name"], a["location"]), name),
                                simapro_units[a["unit"]],
                                1.0,
                                "not defined",
                                category,
                            ]
                        )

                    if ei_version == "uvek":
                        rows.append(
                            [
                                a["name"] + "/" + a["location"] + " U",
                                simapro_units[a["unit"]],
                                1.0,
                                "not defined",
                                category,
                            ]
                        )

                if item == "Products":
                    for e in a["exchanges"]:
                        if e["type"] == "production":
                            name = (
                                e["name"].capitalize()
                                + " {"
                                + e.get("location", "GLO")
                                + "}"
                                + "| Cut-off, U"
                            )

                            if ei_version in ("3.5", "3.6"):
                                rows.append(
                                    [
                                        dict_tech.get((a["name"], a["location"]), name),
                                        simapro_units[a["unit"]],
                                        1.0,
                                        "100%",
                                        "not defined",
                                        category,
                                    ]
                                )

                            if ei_version == "uvek":
                                rows.append(
                                    [
                                        a["name"] + "/" + a["location"] + " U",
                                        simapro_units[a["unit"]],
                                        1.0,
                                        "100%",
                                        "not defined",
                                        category,
                                    ]
                                )

                if item == "Materials/fuels":
                    for e in a["exchanges"]:
                        if e["type"] == "technosphere":
                            if ei_version in ("3.5", "3.6"):
                                if (
                                    not any(
                                        i.lower() in e["name"].lower()
                                        for i in (
                                            "waste",
                                            "emissions",
                                            "treatment",
                                            "scrap",
                                            "used powertrain",
                                            "disposal",
                                            "sludge",
                                            "used li-ion",
                                            "mineral oil storage",
                                        )
                                    )
                                    or any(
                                        i in e["name"]
                                        for i in [
                                            "from municipal waste incineration",
                                            "municipal solid waste, incineration",
                                            "Biomethane",
                                            "biogas upgrading",
                                            "anaerobic digestion, with biogenic carbon uptake",
                                        ]
                                    )
                                    or any(
                                        i.lower() in e["reference product"].lower()
                                        for i in [
                                            "electricity",
                                        ]
                                    )
                                ):

                                    if ei_version == "3.6":
                                        (
                                            e["name"],
                                            e["location"],
                                            e["unit"],
                                            e["reference product"],
                                        ) = self.map_37_to_36.get(
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                        )
                                    if ei_version == "3.5":
                                        (
                                            e["name"],
                                            e["location"],
                                            e["unit"],
                                            e["reference product"],
                                        ) = self.map_37_to_35.get(
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                        )

                                    name = (
                                        e["name"].capitalize()
                                        + " {"
                                        + e.get("location", "GLO")
                                        + "}"
                                    )

                                    if name not in list_own_datasets:
                                        name = (
                                            e["reference product"].capitalize()
                                            + " {"
                                            + e.get("location", "GLO")
                                            + "}"
                                        )

                                        if "market" in e["name"]:
                                            name += (
                                                "| market for "
                                                + e["reference product"].lower()
                                                + " "
                                            )
                                        if "market group" in e["name"]:
                                            name += (
                                                "| market group for "
                                                + e["reference product"].lower()
                                                + " "
                                            )

                                        if "production" in e["name"]:
                                            if (
                                                len(e["reference product"].split(", "))
                                                > 1
                                            ):
                                                name += (
                                                    "| "
                                                    + e["reference product"].split(
                                                        ", "
                                                    )[0]
                                                    + " production, "
                                                    + e["reference product"].split(
                                                        ", "
                                                    )[1]
                                                    + " "
                                                )

                                    rows.append(
                                        [
                                            dict_tech.get(
                                                (e["name"], e["location"]), name
                                            )
                                            + "| Cut-off, U",
                                            simapro_units[e["unit"]],
                                            "{:.3E}".format(e["amount"]),
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                            if ei_version == "uvek":
                                if (
                                    not any(
                                        i.lower() in e["name"].lower()
                                        for i in (
                                            "waste",
                                            "emissions",
                                            "treatment",
                                            "scrap",
                                            "used powertrain",
                                            "disposal",
                                            "used passenger car",
                                            "used electric passenger car",
                                            "anaerobic digestion, with biogenic carbon uptake",
                                            "mineral oil storage",
                                        )
                                    )
                                    or any(
                                        i in e["name"]
                                        for i in [
                                            "from municipal waste incineration",
                                            "aluminium scrap, new",
                                            "brake wear emissions",
                                            "tyre wear emissions",
                                            "road wear emissions",
                                            "used powertrain from electric passenger car",
                                        ]
                                    )
                                    or (
                                        "municipal solid waste, incineration"
                                        in e["name"]
                                        and e["unit"] == "kilowatt hour"
                                    )
                                ):

                                    if e["name"] not in [i["name"] for i in data]:

                                        try:
                                            name = self.map_36_to_uvek_for_simapro[
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ]
                                        except:
                                            print(
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            )
                                            name = ""

                                    else:
                                        name = e["name"] + "/" + e["location"] + " U"

                                    uvek_multiplication_factors = {
                                        "market for heat, from steam, in chemical industry": 1
                                        / 2.257,
                                        "steam production, as energy carrier, in chemical industry": 1
                                        / 2.257,
                                        "market group for natural gas, high pressure": 0.842,
                                        "market for natural gas, high pressure": 0.842,
                                        "market for natural gas, high pressure, vehicle grade": 0.842,
                                        "market for chemical factory": 1 / 12.6e6,
                                        "market for used powertrain from electric passenger car, manual dismantling": -1,
                                    }

                                    uvek_units = {
                                        "market for chemical factory": "unit",
                                        "market for heat, from steam, in chemical industry": "kilogram",
                                        "steam production, as energy carrier, in chemical industry": "kilogram",
                                        "market for manual dismantling of used electric passenger car": "kilogram",
                                        "market group for natural gas, high pressure": "kilogram",
                                        "market for natural gas, high pressure": "kilogram",
                                        "market for transport, pipeline, onshore, petroleum": "kilometer",
                                        "market for natural gas, high pressure, vehicle grade": "megajoule",
                                    }

                                    if e["name"] in uvek_multiplication_factors:
                                        factor = uvek_multiplication_factors[e["name"]]
                                    else:
                                        factor = 1

                                    if e["name"] in uvek_units:
                                        e["unit"] = uvek_units[e["name"]]

                                    rows.append(
                                        [
                                            name,
                                            simapro_units[e["unit"]],
                                            "{:.3E}".format(e["amount"] * factor),
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                if item == "Resources":
                    for e in a["exchanges"]:
                        if (
                            e["type"] == "biosphere"
                            and e["categories"][0] == "natural resource"
                        ):
                            if e["name"] not in simapro_biosphere_flows_to_remove:
                                rows.append(
                                    [
                                        dict_bio[e["name"]],
                                        "",
                                        simapro_units[e["unit"]],
                                        "{:.3E}".format(e["amount"]),
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )

                if item == "Emissions to air":
                    for e in a["exchanges"]:
                        if (
                            e["type"] == "biosphere" and e["categories"][0] == "air"
                        ) or e["name"] in [
                            "Carbon dioxide, from soil or biomass stock",
                            "Carbon dioxide, to soil or biomass stock",
                        ]:
                            if e["name"] not in simapro_biosphere_flows_to_remove:

                                if e["name"].lower() == "water":
                                    e["unit"] = "kilogram"
                                    e["amount"] /= 1000

                                if e["name"] in [
                                    "Carbon dioxide, to soil or biomass stock"
                                ]:
                                    rows.append(
                                        [
                                            dict_bio.get(e["name"], e["name"]),
                                            "",
                                            simapro_units[e["unit"]],
                                            "{:.3E}".format(e["amount"] * -1),
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                                else:
                                    rows.append(
                                        [
                                            dict_bio.get(e["name"], e["name"]),
                                            "",
                                            simapro_units[e["unit"]],
                                            "{:.3E}".format(e["amount"]),
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                if item == "Emissions to water":
                    for e in a["exchanges"]:
                        if e["type"] == "biosphere" and e["categories"][0] == "water":
                            if e["name"] not in simapro_biosphere_flows_to_remove:
                                if e["name"].lower() == "water":
                                    e["unit"] = "kilogram"
                                    e["amount"] /= 1000

                                rows.append(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        "",
                                        simapro_units[e["unit"]],
                                        "{:.3E}".format(e["amount"]),
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )

                if item == "Emissions to soil":
                    for e in a["exchanges"]:
                        if (
                            e["type"] == "biosphere" and e["categories"][0] == "soil"
                        ) and e["name"] not in [
                            "Carbon dioxide, from soil or biomass stock",
                            "Carbon dioxide, to soil or biomass stock",
                        ]:
                            if e["name"] not in simapro_biosphere_flows_to_remove:
                                rows.append(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        "",
                                        simapro_units[e["unit"]],
                                        "{:.3E}".format(e["amount"]),
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )

                if item == "Waste to treatment":
                    for e in a["exchanges"]:
                        is_waste = False
                        if e["type"] == "technosphere":

                            # We check if this is indeed a waste treatment activity
                            if e["name"] in self.references:
                                if self.references[e["name"]] == "waste treatment":
                                    is_waste = True
                            else:
                                if (
                                    any(
                                        i.lower() in e["name"].lower()
                                        for i in (
                                            " waste ",
                                            "emissions",
                                            "treatment",
                                            "scrap",
                                            "used powertrain",
                                            "used passenger car",
                                            "used electric passenger car",
                                            "municipal solid waste",
                                            "disposal",
                                            "rainwater mineral oil",
                                            "sludge",
                                            "used li-ion",
                                        )
                                    )
                                    and not any(
                                        i.lower() in e["name"].lower()
                                        for i in (
                                            "anaerobic",
                                            "cooking",
                                            "heat",
                                            "manual dismantling",
                                        )
                                    )
                                    and e["unit"] not in ["kilowatt hour", "megajoule"]
                                ):
                                    is_waste = True

                            # Yes, it is a waste treatment activity
                            if is_waste:

                                name = ""

                                # In SimaPro, waste inputs are positive numbers
                                if e["amount"] < 0:
                                    e["amount"] *= -1

                                if ei_version in ("3.5", "3.6"):

                                    if ei_version == "3.6":
                                        (
                                            e["name"],
                                            e["location"],
                                            e["unit"],
                                            e["reference product"],
                                        ) = self.map_37_to_36.get(
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                        )
                                    if ei_version == "3.5":
                                        (
                                            e["name"],
                                            e["location"],
                                            e["unit"],
                                            e["reference product"],
                                        ) = self.map_37_to_35.get(
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                            (
                                                e["name"],
                                                e["location"],
                                                e["unit"],
                                                e["reference product"],
                                            ),
                                        )

                                    name = dict_tech.get(
                                        (e["name"], e["location"]),
                                        e["name"] + " {" + e["location"] + "}",
                                    )

                                    rows.append(
                                        [
                                            name + "| Cut-off, U",
                                            simapro_units[e["unit"]],
                                            "{:.3E}".format(e["amount"]),
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                                if ei_version == "uvek":

                                    if not any(
                                        i in e["name"].lower()
                                        for i in [
                                            "brake wear",
                                            "tyre wear",
                                            "road wear",
                                            "aluminium scrap, new",
                                            "used powertrain from electric passenger car",
                                        ]
                                    ):

                                        uvek_multiplication_factors = {
                                            "market for manual dismantling of used electric passenger car": 1
                                            / 1200,
                                            "manual dismantling of used passenger car with internal combustion engine": 1
                                            / 1200,
                                            "market for manual dismantling of used passenger car with internal combustion engine": 1
                                            / 1200,
                                        }

                                        if e["name"] in uvek_multiplication_factors:
                                            factor = uvek_multiplication_factors[
                                                e["name"]
                                            ]
                                        else:
                                            factor = 1

                                        if e["name"] not in [i["name"] for i in data]:
                                            try:
                                                name = self.map_36_to_uvek_for_simapro[
                                                    e["name"],
                                                    e["location"],
                                                    e["unit"],
                                                    e["reference product"],
                                                ]
                                            except:
                                                print(
                                                    e["name"],
                                                    e["location"],
                                                    e["unit"],
                                                    e["reference product"],
                                                )
                                                name = ""

                                        else:
                                            name = (
                                                e["name"] + "/" + e["location"] + " U"
                                            )

                                        rows.append(
                                            [
                                                name,
                                                simapro_units[e["unit"]],
                                                "{:.3E}".format(e["amount"] * factor),
                                                "undefined",
                                                0,
                                                0,
                                                0,
                                            ]
                                        )

                rows.append([])

        # System description
        rows.append(["System description"])
        rows.append([])
        rows.append(["Name"])
        rows.append(["carculator"])
        rows.append([])
        rows.append(["Category"])
        rows.append(["transport"])
        rows.append([])
        rows.append(["Description"])
        rows.append(
            [
                "Prospective life cycle assessment model for passenger cars developed by PSI"
            ]
        )
        rows.append([])
        rows.append(["Cut-off rules"])
        rows.append(
            [
                "All environmentally-relevant flows are included, as far as the authors knowledge permits."
                "Also, residual material (e.g., biomass residue) and energy (e.g., waste heat) "
                "come free of burden, except for the necessary steps to make it reusable"
                " (transport, conditioning, etc.)."
            ]
        )
        rows.append([])
        rows.append(["Energy model"])
        rows.append(
            [
                "The energy consumption of vehicles calculated based on a physics model, including "
                "inertia, rolling resistance, aerodynamic drag, road gradient, etc."
            ]
        )
        rows.append([])
        rows.append(["Transport model"])
        rows.append(["Based on Sacchi et al. 2020 (in review)"])
        rows.append([])
        rows.append(["Allocation rules"])
        rows.append(
            [
                "The system modeling is attributional. In the instance of joint-production, the allocation of "
                "burden between co-products is generally based on the relative economic revenue of "
                "each product, to align with the underlying database ecoinvent cut-off."
            ]
        )
        rows.append(["End"])
        rows.append([])

        # Literature reference
        rows.append(["Literature reference"])
        rows.append([])
        rows.append(["Name"])
        rows.append(["Sacchi et al. 2020"])
        rows.append([])
        rows.append(["Documentation link"])
        rows.append(["https://www.psi.ch/en/ta/preprint"])
        rows.append([])
        rows.append(["Comment"])
        rows.append(
            ["Pre-print available at: https://www.psi.ch/en/media/57994/download"]
        )
        rows.append([])
        rows.append(["Category"])
        rows.append(["carculator"])
        rows.append([])
        rows.append(["Description"])
        description = (
            "carculator: an open-source tool for prospective environmental and "
            "economic life cycle assessment of vehicles. When, Where and How can battery-electric "
            "vehicles help reduce greenhouse gas emissions?\n"
        )
        description += "Romain Sacchi, Christian Bauer and Brian L. Cox\n"
        description += (
            "Submitted to Environmental Science and Technology on November 17th, 2020"
        )

        rows.append([description])

        return rows

    def write_lci_to_bw(
        self,
        presamples: bool,
        ecoinvent_version: str,
        vehicle_specs: xr.DataArray = None,
    ) -> bw2io.importers.base_lci.LCIImporter:
        """
        Return a LCIImporter object with the inventory as `data` attribute.

        :return: LCIImporter object to be imported in a Brightway2 project
        :rtype: bw2io.base_lci.LCIImporter
        """
        if presamples:
            data, array = self.write_lci(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=vehicle_specs,
            )
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return i, array
        else:
            data = self.write_lci(
                presamples=presamples,
                ecoinvent_version=ecoinvent_version,
                vehicle_specs=vehicle_specs,
            )
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return i

    def best_fit_distribution(self, data, bins=200, ax=None):
        import warnings

        import pandas as pd
        import scipy.stats as st

        """
        Model data by finding best fit distribution to data
        Return the most likely value as well as a list of tuples that contains distribution parameters
        """
        # Get histogram of original data
        y, x = np.histogram(data, bins=bins, density=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0

        # Distributions to check
        DISTRIBUTIONS = [
            # st.beta,
            # st.gamma,
            # st.lognorm,
            st.norm,
            # st.t,
            # st.triang,
            # st.uniform,
            # st.weibull_min,
        ]

        # Best holders
        best_distribution = st.norm
        best_params = (0.0, 1.0)
        best_sse = np.inf

        # Estimate distribution parameters from data
        for distribution in DISTRIBUTIONS:

            # Try to fit the distribution
            try:
                # Ignore warnings from data that can't be fit
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")

                    # fit dist to data
                    params = distribution.fit(data)

                    # Separate parts of parameters
                    arg = params[:-2]
                    loc = params[-2]
                    scale = params[-1]

                    # Calculate fitted PDF and error with fit in distribution
                    pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                    sse = np.sum(np.power(y - pdf, 2.0))

                    # if axis pass in add to plot
                    try:
                        if ax:
                            pd.Series(pdf, x).plot(ax=ax)
                    except Exception:
                        pass

                    # identify if this distribution is better
                    if best_sse > sse > 0:
                        best_distribution = distribution
                        best_params = params
                        best_sse = sse

            except Exception:
                pass

        # Lognormal distribution
        if self.uncertainty_id[best_distribution.name] == 2:
            mu, std = st.norm.fit(data)
            return [("uncertainty type", 2), ("scale", std), ("loc", mu)]

        # Normal distribution
        if self.uncertainty_id[best_distribution.name] == 3:
            return [
                ("uncertainty type", 3),
                ("loc", best_params[0]),
                ("scale", best_params[1]),
            ]

        # Uniform distribution
        if self.uncertainty_id[best_distribution.name] == 4:
            return [
                ("uncertainty type", 4),
                ("minimum", best_params[0]),
                ("maximum", (best_params[0] + best_params[1])),
            ]

        # Triangular distribution
        if self.uncertainty_id[best_distribution.name] == 5:
            return [
                ("uncertainty type", 5),
                ("loc", best_params[1]),
                ("minimum", np.min(data)),
                ("maximum", np.max(data)),
            ]

        # Gamma distribution
        if self.uncertainty_id[best_distribution.name] == 9:
            return [
                ("uncertainty type", 9),
                ("shape", best_params[0]),
                ("scale", best_params[2]),
                ("loc", best_params[1]),
            ]

        # Beta distribution
        if self.uncertainty_id[best_distribution.name] == 10:
            return [
                ("uncertainty type", 10),
                ("loc", best_params[0]),
                ("shape", best_params[1]),
                ("scale", best_params[3]),
            ]

        # Weibull distribution
        if self.uncertainty_id[best_distribution.name] == 8:
            return [
                ("uncertainty type", 8),
                ("shape", best_params[0]),
                ("loc", best_params[1]),
                ("scale", best_params[2]),
            ]

        # Student's T distribution
        if self.uncertainty_id[best_distribution.name] == 12:
            return [
                ("uncertainty type", 12),
                ("shape", best_params[0]),
                ("loc", best_params[1]),
                ("scale", best_params[2]),
            ]

        def make_pdf(self, dist, params, size=10000):
            """Generate distributions's Probability Distribution Function"""
            import pandas as pd

            # Separate parts of parameters
            arg = params[:-2]
            loc = params[-2]
            scale = params[-1]

            # Get same start and end points of distribution
            start = (
                dist.ppf(0.01, *arg, loc=loc, scale=scale)
                if arg
                else dist.ppf(0.01, loc=loc, scale=scale)
            )
            end = (
                dist.ppf(0.99, *arg, loc=loc, scale=scale)
                if arg
                else dist.ppf(0.99, loc=loc, scale=scale)
            )

            # Build PDF and turn into pandas Series
            x = np.linspace(start, end, size)
            y = dist.pdf(x, loc=loc, scale=scale, *arg)
            pdf = pd.Series(y, x)

            return pdf
