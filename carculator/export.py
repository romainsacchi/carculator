from bw2io.export.excel import (
    safe_filename,
    xlsxwriter,
    create_valid_worksheet_name,
)
import bw2io
import numpy as np
import os
import pyprind
import uuid
import datetime
import json
import csv
from . import DATA_DIR, __version__


def load_mapping_37_to_36():
    """Load mapping dictionary between ecoinvent 3.7 and 3.6"""

    # Load the matching dictionary
    filename = "ei37_to_ei36.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of activities flows match between ecoinvent 3.7 and 3.6 could not be found."
        )
    with open(filepath) as f:
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

def load_references():
    """Load a dictionary with references of datasets"""

    # Load the matching dictionary
    filename = "references.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of references could not be found."
        )
    with open(filepath) as f:
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
    header, *data = csv_list

    dict_reference = {}
    for row in data:
        name, source, description, special_remark, category_1, category_2 = row
        dict_reference[name] = {
            "source": source,
            "description": description,
            "special remark": special_remark,
            "category 1": category_1,
            "category 2": category_2
        }

    return dict_reference

def load_mapping_37_to_35():
    """Load mapping dictionary between ecoinvent 3.7 and 3.5"""

    # Load the matching dictionary
    filename = "ei37_to_ei35.csv"
    filepath = DATA_DIR / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of activities flows match between ecoinvent 3.7 and 3.5 could not be found."
        )
    with open(filepath) as f:
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
    (_, _, *header), *data = csv_list

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

class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices, db_name="carculator export"):
        self.array = array
        self.indices = indices
        self.rename_vehicles()
        self.db_name = db_name
        self.references = load_references()
        # See https://docs.brightwaylca.org/intro.html#uncertainty-type
        self.uncertainty_ID = {
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
        self.map_remind_ecoinvent = {
            (
                "market group for electricity, high voltage",
                "EUR",
                "kilowatt hour",
                "electricity, high voltage",
            ): (
                "market group for electricity, high voltage",
                "ENTSO-E",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            (
                "market group for electricity, medium voltage",
                "EUR",
                "kilowatt hour",
                "electricity, medium voltage",
            ): (
                "market group for electricity, medium voltage",
                "ENTSO-E",
                "kilowatt hour",
                "electricity, medium voltage",
            ),
            (
                "market group for electricity, low voltage",
                "EUR",
                "kilowatt hour",
                "electricity, low voltage",
            ): (
                "market group for electricity, low voltage",
                "ENTSO-E",
                "kilowatt hour",
                "electricity, low voltage",
            ),
            (
                "market group for electricity, medium voltage",
                "JPN",
                "kilowatt hour",
                "electricity, medium voltage",
            ): (
                "market for electricity, medium voltage",
                "JP",
                "kilowatt hour",
                "electricity, medium voltage",
            ),
            (
                "market group for electricity, high voltage",
                "World",
                "kilowatt hour",
                "electricity, high voltage",
            ): (
                "market group for electricity, high voltage",
                "GLO",
                "kilowatt hour",
                "electricity, high voltage",
            ),
            (
                "market group for electricity, medium voltage",
                "World",
                "kilowatt hour",
                "electricity, medium voltage",
            ): (
                "market group for electricity, medium voltage",
                "GLO",
                "kilowatt hour",
                "electricity, medium voltage",
            ),
            (
                "market group for electricity, low voltage",
                "World",
                "kilowatt hour",
                "electricity, low voltage",
            ): (
                "market group for electricity, low voltage",
                "GLO",
                "kilowatt hour",
                "electricity, low voltage",
            ),
            (
                "cement production, Portland",
                "EUR",
                "kilogram",
                "cement, Portland",
            ): (
                "cement production, Portland",
                "CH",
                "kilogram",
                "cement, Portland",
            ),

        }
        self.map_ecoinvent_remind= {
            (
                "biogas upgrading - sewage sludge - amine scrubbing - best",
                "CH",
                "kilogram",
                "biogas upgrading - sewage sludge - amine scrubbing - best",
            ): (
                "biogas upgrading - sewage sludge - amine scrubbing - best",
                "RER",
                "kilogram",
                "biogas upgrading - sewage sludge - amine scrubbing - best",
            ),

        }
        self.map_37_to_36 = load_mapping_37_to_36()
        self.map_37_to_35 = load_mapping_37_to_35()
        self.map_36_to_uvek = self.load_mapping_36_to_uvek()
        self.map_36_to_uvek_for_simapro = self.load_mapping_36_to_uvek_for_simapro()
        self.tags = self.load_tags()

    def rename_vehicles(self):

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
            for key in d_names:
                if key in value[0]:
                    new_val = list(value)
                    new_val[0] = new_val[0].replace(key, d_names[key])
                    self.indices[k] = tuple(new_val)

    def load_tags(self):
        """Loads dictionary of tags for further use in BW2"""

        filename = "tags.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of tags could not be found."
            )
        with open(filepath) as f:
            csv_list = [
                [val.strip() for val in r.split(";")] for r in f.readlines()
            ]
        data = csv_list

        dict_tags = {}
        for row in data:
            name, tag = row
            dict_tags[name] = tag

        return dict_tags

    def load_mapping_36_to_uvek(self):
        """Load mapping dictionary between ecoinvent 3.6 and UVEK"""

        # Load the matching dictionary between ecoinvent and Simapro biosphere flows
        filename = "uvek_mapping.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of activities flows match between ecoinvent 3.6 and UVEK could not be found."
            )
        with open(filepath) as f:
            csv_list = [
                [val.strip() for val in r.split(";")] for r in f.readlines()
            ]
        (_, _, *header), *data = csv_list

        dict_uvek = {}
        for row in data:
            name, location, unit, ref_prod, uvek_name, uvek_loc, uvek_unit, uvek_ref_prod, simapro_name = row
            dict_uvek[(name, ref_prod, unit, location)] = (uvek_name, uvek_ref_prod, uvek_unit, uvek_loc)
        return dict_uvek

    def load_mapping_36_to_uvek_for_simapro(self):
        """Load mapping dictionary between ecoinvent 3.6 and UVEK for Simapro name format"""

        # Load the matching dictionary between ecoinvent and Simapro biosphere flows
        filename = "uvek_mapping.csv"
        filepath = DATA_DIR / filename
        if not filepath.is_file():
            raise FileNotFoundError(
                "The dictionary of activities flows match between ecoinvent 3.6 and UVEK could not be found."
            )
        with open(filepath) as f:
            csv_list = [
                [val.strip() for val in r.split(";")] for r in f.readlines()
            ]
        (_, _, *header), *data = csv_list

        dict_uvek = {}
        for row in data:
            name, location, unit, ref_prod, uvek_name, uvek_loc, uvek_unit, uvek_ref_prod, simapro_name = row
            dict_uvek[(name, location, unit, ref_prod)] = simapro_name

        return dict_uvek

    def write_lci(self, presamples, ecoinvent_compatibility, ecoinvent_version, forbidden_activities=None):
        """
        Return the inventory as a dictionary
        If if there several values for one exchange, uncertainty information is generated.
        If `presamples` is True, returns the inventory as well as a `presamples` matrix.
        If `presamples` is False, returns the inventory with characterized uncertainty information.
        If `ecoinvent_compatibility` is True, the inventory is made compatible with ecoinvent. If False,
        the inventory is compatible with the REMIND-ecoinvent hybrid database output of the `rmnd_lca` library.

        :returns: a dictionary that contains all the exchanges
        :rtype: dict
        """

        # List of activities that are already part of the REMIND-ecoinvent database.
        # They should not appear in the exported inventories, otherwise they will be duplicates
        activities_to_be_removed = [
            "algae cultivation | algae broth production",
            "algae harvesting| dry algae production",
            "transport, pipeline, supercritical CO2, 200km w/o recompression",
            "Ethanol from maize starch",
            "Natural Gas provision (at medium pressure grid) {RER}, EU mix",
            "woodchips from forestry residues",
            "Ethanol from wheat straw pellets",
            "straw pellets",
            "Biodiesel from cooking oil",
            "Straw bales | baling of straw",
            "drilling, deep borehole/m",
            "Sugar beet cultivation {RER} | sugar beet production Europe | Alloc Rec, U",
            "Refined Waste Cooking Oil {RER} | Refining of waste cooking oil Europe | Alloc Rec, U",
            "Ethanol from forest residues",
            "Ethanol from sugarbeet",
            "pipeline, supercritical CO2/km",
            "Biodiesel from algae",
            "Maize cultivation, drying and storage {RER} | Maize production Europe | Alloc Rec, U",
            "Fischer Tropsch reactor and upgrading plant, construction",
            "Walls and foundations, for hydrogen refuelling station",
            "container, with pipes and fittings, for diaphragm compressor",
            "RWGS tank construction",
            "storage module, high pressure, at fuelling station",
            "pumps, carbon dioxide capture process",
            "PEM electrolyzer, Operation and Maintenance",
            "heat exchanger, carbon dioxide capture process",
            "biogas upgrading - sewage sludge - amine scrubbing - best",
            "Hydrogen refuelling station, SMR",
            "transformer and rectifier unit, for electrolyzer",
            "PEM electrolyzer, ACDC Converter",
            "carbon dioxide, captured from atmosphere",
            "PEM electrolyzer, Balance of Plant",
            "Sabatier reaction methanation unit",
            "PEM electrolyzer, Stack",
            "hot water tank, carbon dioxide capture process",
            "cooling unit, carbon dioxide capture process",
            "diaphragm compressor module, high pressure",
            "carbon dioxide capture system",
            "Hydrogen dispenser, for gaseous hydrogen",
            "diaphragms, for diaphragm compressor",
            "MTG production facility, construction",
            "Disposal, hydrogen fuelling station",
            "production of 2 wt-% potassium iodide solution",
            "production of nickle-based catalyst for methanation",
            "wiring and tubing, carbon dioxide capture process",
            "control panel, carbon dioxide capture process",
            "adsorption and desorption unit, carbon dioxide capture process",
            "Buffer tank",
            "frequency converter, for diaphragm compressor",
            'Hydrogen, gaseous, 30 bar, from hard coal gasification and reforming, at coal gasification plant',
            #'Methanol distillation',
            'CO2 storage/at H2 production plant, pre, pipeline 200km, storage 1000m',
            #'Syngas, RWGS, Production',
            'softwood forestry, mixed species, sustainable forest management, CF = -1',
            'hardwood forestry, mixed species, sustainable forest management, CF = -1',
            'Hydrogen, gaseous, 25 bar, from dual fluidised bed gasification of woody biomass with CCS, at gasification plant',
            'market for wood chips, wet, measured as dry mass, CF = -1',
            'Hydrogen, gaseous, 700 bar, from dual fluidised bed gasification of woody biomass with CCS, at fuelling station',
            'SMR BM, HT+LT, + CCS (MDEA), 98 (average), digestate incineration, 26 bar',
            'Hydrogen, gaseous, 700 bar, from SMR of biogas, at fuelling station',
            'SMR NG + CCS (MDEA), 98 (average), 25 bar',
            'SMR BM, HT+LT, with digestate incineration, 26 bar',
            'Hydrogen, gaseous, 700 bar, from dual fluidised bed gasification of woody biomass, at fuelling station',
            'Hydrogen, gaseous, 700 bar, from SMR of biogas with CCS, at fuelling station',
            'SMR NG + CCS (MDEA), 98 (average), 700 bar',
            'Hydrogen, gaseous, 25 bar, from dual fluidised bed gasification of woody biomass, at gasification plant',
            #'Methanol Synthesis',
            #'Diesel production, synthetic, Fischer Tropsch process',
            #'Gasoline production, synthetic, from methanol',
            'Crude vegetable oil | oil mill: extraction of vegetable oil from rapeseed | Alloc Rec, U',
            'biomethane from biogas upgrading - biowaste - amine scrubbing, best - with biogenic carbon uptake, lower bound C sequestration, digestate incineration',
            'Plant oil from crude oil | refining of vegetable oil from oil palm|',
            'ATR BM, with digestate incineration, 25 bar',
            'Hydrogen, gaseous, 25 bar, from gasification of woody biomass in oxy-fired entrained flow gasifier, at gasification plant',
            'Ethanol from wheat grains',
            'Wheat grain cultivation, drying and storage {RER} | wheat grain production Europe | Alloc Rec, U',
            'Hydrogen, gaseous, 700 bar, ATR of NG, with CCS, at fuelling station',
            'ATR BM + CCS (MDEA), 98 (average), with digestate incineration, 25 bar',
            'SMR NG, 25 bar',
            'Plant oil production | refining of crude vegetable oil from rapeseed| Alloc Rec, U',
            'ethanol without biogas',
            'Biodiesel from rapeseed oil',
            'Hydrogen, gaseous, 700 bar, ATR of NG, at fuelling station',
            'Waste Cooking Oil',
            'Carbon monoxide, from RWGS',
            'Oil Palm Tree Cultivation {RER} | Fresh Fruit Bunches (FFBs) production | Alloc Rec, U',
            'ATR NG, 25 bar',
            'Rapeseed cultivation {RER} | rapeseed production Europe | Alloc Rec, U',
            'Hydrogen, gaseous, 700 bar, from ATR of biogas with CCS, at fuelling station',
            'Hydrogen, gaseous, 25 bar, from gasification of woody biomass in oxy-fired entrained flow gasifier, with CCS, at gasification plant',
            'SMR NG, 700 bar',
            'Hydrogen, gaseous, 700 bar, from SMR of NG, with CCS, at fuelling station',
            'ATR NG + CCS (MDEA), 98 (average), 25 bar',
            'Hydrogen, gaseous, 700 bar, from ATR of biogas, at fuelling station',
            'Biodiesel from palm oil',
            'Hydrogen, gaseous, 700 bar, from SMR of NG, at fuelling station',
            'treatment of biowaste by anaerobic digestion, with biogenic carbon uptake, lower bound C sequestration, digestate incineration',
            'Crude Palm Oil extraction from FFBs {RER} |oil mill|',
            'Hydrogen, gaseous, 700 bar, from gasification of woody biomass in oxy-fired entrained flow gasifier, at fuelling station',
            'hardwood forestry, birch, sustainable forest management_CF = -1',
            'Hydrogen, gaseous, 700 bar, from electrolysis, at fuelling station',
            'Diesel production, synthetic, Fischer Tropsch process',
            'softwood forestry, pine, sustainable forest management_CF = -1',
            'softwood forestry, spruce, sustainable forest management_CF = -1',
            'Hydrogen, gaseous, 25 bar, from electrolysis',
            'Methanol distillation',
            'softwood forestry, spruce, sustainable forest management_CF = -1',
            'hardwood forestry, oak, sustainable forest management_CF = -1',
            'softwood forestry, pine, sustainable forest management_CF = -1',
            'Gasoline production, synthetic, from methanol',
            'Methanol Synthesis',
            'Hydrogen, gaseous, 700 bar, from gasification of woody biomass in oxy-fired entrained flow gasifier, with CCS, at fuelling station',
            'hardwood forestry, beech, sustainable forest management_CF = -1',
            'Syngas, RWGS, Production'
        ]

        if isinstance(forbidden_activities, list):
            activities_to_be_removed.extend(forbidden_activities)

        uvek_activities_to_remove = [
            "market for activated carbon, granular",
            "market for iodine",
            "market for manganese sulfate",
            "market for molybdenum trioxide",
            "market for nickel sulfate",
            "market for soda ash, light, crystalline, heptahydrate",
        ]

        ei35_activities_to_remove = [
            "latex production"
        ]

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

                # If ecoinvent_compatibility==False and the activity name is part of the list
                if tuple_output[0] in activities_to_be_removed and not ecoinvent_compatibility:
                    break

                if ecoinvent_compatibility == False:
                    tuple_output = self.map_ecoinvent_remind.get(
                        tuple_output, tuple_output
                    )
                    tuple_input = self.map_ecoinvent_remind.get(
                        tuple_input, tuple_input
                    )

                if ecoinvent_compatibility == True:

                    tuple_output = self.map_remind_ecoinvent.get(
                        tuple_output, tuple_output
                    )
                    tuple_input = self.map_remind_ecoinvent.get(
                        tuple_input, tuple_input
                    )

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
                        else:
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
                    if presamples == True:
                        # Generate pre-sampled values
                        amount = np.median(self.array[:, row, col])  * mult_factor
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
                    # else:
                    #    # Generate uncertainty distribution parameters
                    #    amount = np.median(self.array[:, row, col])
                    #    uncertainty = self.best_fit_distribution(self.array[:, row, col] * -1)

                # Look for a tag, if any
                tag = [self.tags[t] for t in list(self.tags.keys()) if t in tuple_input[0]]
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
                            "tag":tag
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
                            "tag":tag
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
                            "tag":tag
                        }
                    )
                    list_exc[-1].update(uncertainty)

            else:

                # Look for a tag, if any
                tag = [self.tags[t] for t in list(self.tags.keys()) if t in tuple_output[0]]
                if len(tag) > 0:
                    tag = tag[0]
                else:
                    tag = "other"

                if "transport, passenger car" in tuple_output[0]:
                    source = self.references["transport, passenger car"]["source"]
                    description = self.references["transport, passenger car"]["description"]
                    special_remark = self.references["transport, passenger car"]["special remark"]

                elif "Passenger car" in tuple_output[0]:
                    source = self.references["passenger car"]["source"]
                    description = self.references["passenger car"]["description"]
                    special_remark = self.references["passenger car"]["special remark"]

                else:
                    source = self.references[tuple_output[0]]["source"] if tuple_output[0] in self.references else ""
                    description = self.references[tuple_output[0]]["description"] if tuple_output[0] in self.references else ""
                    special_remark = self.references[tuple_output[0]]["special remark"] if tuple_output[0] in self.references else ""

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
                    }
                )
        if presamples:
            return (list_act, presamples_matrix)
        else:
            return list_act

    def write_lci_to_excel(
        self,
        directory,
        ecoinvent_compatibility,
        ecoinvent_version,
        software_compatibility,
        filename=None,
        forbidden_activities=None
    ):
        """
        Export an Excel file that can be consumed by the software defined in `software_compatibility`.

        :param directory: str. path to export the file to.
        :param ecoinvent_compatibility: bool. If True, the inventory is compatible with ecoinvent. If False, the inventory is compatible with REMIND-ecoinvent.
        :param ecoinvent_version: str. "3.5", "3.6" or "uvek"
        :param software_compatibility: str. "brightway2" or "simapro"
        :returns: returns the file path of the exported inventory.
        :rtype: str.
        """

        if software_compatibility == "brightway2":
            if filename is None:
                safe_name = safe_filename(
                    "carculator_inventory_export_{}_brightway2".format(
                        str(datetime.date.today())
                    ),
                    False,
                ) + ".xlsx"
            else:
                safe_name = safe_filename(
                    filename,
                    False,
                ) + ".xlsx"
        else:
            safe_name = safe_filename(
                "carculator_inventory_export_{}_simapro".format(
                    str(datetime.date.today())
                ),
                False,
            ) + ".csv"

        if directory is None:
            filepath_export = safe_name
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)
            filepath_export = os.path.join(directory, safe_name)

        list_act = self.write_lci(False, ecoinvent_compatibility, ecoinvent_version, forbidden_activities)

        if software_compatibility == "brightway2":
            data = []
            data.extend((["Database", self.db_name], ("format", "Excel spreadsheet")))
            data.append([])

            for k in list_act:
                if k.get("exchanges"):
                    data.extend(
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
                                "tag"
                            ],
                        )
                    )

                    for e in k["exchanges"]:
                        data.append(
                            [
                                e["name"],
                                float(e["amount"]),
                                e["database"],
                                e.get("location", "None"),
                                e["unit"],
                                "::".join(e.get("categories", ())),
                                e["type"],
                                e.get("reference product"),
                                e.get("tag", "other")
                            ]
                        )
                else:
                    data.extend(
                        (
                            ["Activity", k["name"]],
                            ("type", "biosphere"),
                            ("unit", k["unit"]),
                            ("worksheet name", "None"),
                        )
                    )
                data.append([])

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
            print("Inventories exported to {}.".format(filepath_export))
            workbook.close()

        else:

            # Load the matching dictionary between ecoinvent and Simapro biosphere flows
            filename = "simapro-biosphere.json"
            filepath = DATA_DIR / filename
            if not filepath.is_file():
                raise FileNotFoundError(
                    "The dictionary of biosphere flow match between ecoinvent and Simapro could not be found."
                )
            with open(filepath) as json_file:
                data = json.load(json_file)
            dict_bio = {}
            for d in data:
                dict_bio[d[2]] = d[1]

            # Load the matching dictionary between ecoinvent and Simapro product flows
            filename = "simapro-technosphere-3.5.csv"
            filepath = DATA_DIR / filename
            with open(filepath) as f:
                csv_list = [
                    [val.strip() for val in r.split(";")] for r in f.readlines()
                ]
            (_, _, *header), *data = csv_list

            dict_tech = {}
            for row in data:
                name, location, simapro_name = row
                simapro_name = simapro_name.split("|")[:2]
                dict_tech[(name, location)] = ("|").join(simapro_name)

            headers = [

                "{SimaPro 9.1.1.1}",
                "{processes}",
                "{Project: carculator import" + f"{datetime.datetime.today():%d.%m.%Y}" + "}",
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
                "Comment",
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

            with open(filepath_export, "w", newline="") as csvFile:
                writer = csv.writer(csvFile, delimiter=";")
                for item in headers:
                    writer.writerow([item])
                writer.writerow([])

                list_own_datasets = []

                database = "carculator"

                for a in list_act:
                    list_own_datasets.append(
                        a["name"].capitalize()
                        + " {"
                        + a.get("location", "GLO")
                        + "}"
                    )

                for a in list_act:
                    for item in fields:

                        if any(i in a["name"].lower() for i in ["disposal", "treatment"]) \
                                and item == "Products":
                            continue

                        if not any(i in a["name"].lower() for i in ["disposal", "treatment"]) \
                                and item == "Waste treatment":
                            continue

                        writer.writerow([item])

                        if item == "Process name":

                            if ecoinvent_version in ("3.5", "3.6"):
                                name = (
                                        a["name"].capitalize()
                                        + " {"
                                        + a.get("location", "GLO")
                                        + "}"
                                        + "| Cut-off, U"
                                )

                            if ecoinvent_version == "uvek":
                                name = a["name"] + "/" + a["location"] + " U"

                            writer.writerow([name])

                        if item == "Type":
                            writer.writerow(["Unit process"])

                        if item == "Comment":
                            if a["name"] in self.references:

                                string = "Originally published in: "
                                string += self.references[a["name"]]["source"]

                                if self.references[a["name"]]["description"] != "":
                                    string += " Description: "
                                    string += self.references[a["name"]]["description"]

                                if self.references[a["name"]]["special remark"] != "":
                                    string += " Special remark(s): "
                                    string += self.references[a["name"]]["special remark"]

                                writer.writerow([string])
                                string = ""

                            else:

                                if "transport, passenger car" in a["name"]:
                                    key = "transport, passenger car"

                                if "Passenger car" in a["name"]:
                                    key = "passenger car"

                                if "fuel supply for" in a["name"]:
                                    key = "fuel supply"
                                if "electricity supply for" in a["name"]:
                                    key = "electricity supply"
                                if "electricity market for fuel" in a["name"]:
                                    key = "electricity market for fuel"
                                if "electricity market for energy storage" in a["name"]:
                                    key = "electricity market for energy storage"

                                if self.references[key]["source"] != "":
                                    string = "Originally published in: "
                                    string += self.references[key]["source"]

                                if self.references[key]["description"] != "":
                                    string += " Description: "
                                    string += self.references[key]["description"]

                                if self.references[key]["special remark"] != "":
                                    string += " Special remark(s): "
                                    string += self.references[key]["special remark"]

                                writer.writerow([string])
                                string, key = ("", "")

                        if item == "Category type":
                            if any(i in a["name"].lower() for i in ["disposal", "treatment"]):
                                name = "waste treatment"
                            else:
                                name = "transport"
                            writer.writerow([name])

                        if item == "Generator":
                            writer.writerow(["carculator " + str(__version__)])

                        if item == "Geography":
                            writer.writerow([a["location"]])

                        if item == "Time Period":
                            writer.writerow(
                                ["Between 2010 and 2020. Extrapolated to the selected years."]
                            )

                        if item == "Date":
                            writer.writerow([f"{datetime.datetime.today():%d.%m.%Y}"])

                        if item in ("Cut off rules",
                                    "Capital goods",
                                    "Technology",
                                    "Representativeness",
                                    "Boundary with nature"):
                            writer.writerow(["Unspecified"])

                        if item == "Infrastructure":
                            writer.writerow(["Yes"])

                        if item == "External documents":
                            writer.writerow(["https://carculator.psi.ch"])

                        if item in ("System description"):
                            writer.writerow(["carculator"])

                        if item in ("Allocation rules"):
                            writer.writerow(
                                ["In the instance of joint-production, allocation of process burden based on"
                                 "economic relative revenue of each co-product."])

                        if item == "Literature references":
                            writer.writerow(["Sacchi et al. 2020"])

                        if item == "Collection method":
                            writer.writerow(
                                [
                                    "Modeling and assumptions: https://carculator.readthedocs.io/en/latest/modeling.html"
                                ]
                            )

                        if item == "Verification":
                            writer.writerow(["In review. Susceptible to change."])

                        if item == "Waste treatment":
                            if ecoinvent_version in ("3.5", "3.6"):
                                writer.writerow(
                                    [
                                        dict_tech.get(
                                            (a["name"], a["location"]), name
                                        ),
                                        simapro_units[a["unit"]],
                                        1.0,
                                        "not defined",
                                        category,
                                    ]
                                )

                            if ecoinvent_version == "uvek":
                                writer.writerow(
                                    [
                                        a["name"] + "/" + a["location"] + " U",
                                        simapro_units[a["unit"]],
                                        1.0,
                                        "not defined",
                                        category,
                                    ])

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

                                    category = database

                                    if a["name"] in self.references:
                                        if self.references[e["name"]]["category 1"] != "":
                                            category += r"\ ".strip() + self.references[e["name"]]["category 1"]

                                        if self.references[e["name"]]["category 2"] != "":
                                            category += r"\ ".strip() + self.references[e["name"]]["category 2"]
                                    else:
                                        if "transport, " in a["name"] and "kilometer" in a["unit"]:
                                            category += r"\ ".strip() + "passenger cars" + r"\ ".strip() + "transport"

                                        if "Passenger car" in a["name"]:
                                            category += r"\ ".strip() + "passenger cars" + r"\ ".strip() + "vehicles"

                                        if "ICEV" in a["name"]:
                                            category += r"\ ".strip() + "combustion"

                                        if " HEV" in a["name"]:
                                            category += r"\ ".strip() + "hybrid"

                                        if "PHEV" in a["name"]:
                                            category += r"\ ".strip() + "hybrid plugin"

                                        if any(i in a["name"] for i in ("BEV", "FCEV")):
                                            category += r"\ ".strip() + "electric"

                                        if any(i in a["name"].lower() for i in ("fuel supply for",
                                                                                "electricity supply for",
                                                                                "electricity market for")):
                                            category += r"\ ".strip() + "energy mix and fuel blends"

                                    if ecoinvent_version in ("3.5", "3.6"):
                                        writer.writerow(
                                            [
                                                dict_tech.get(
                                                    (a["name"], a["location"]), name
                                                ),
                                                simapro_units[a["unit"]],
                                                1.0,
                                                "100%",
                                                "not defined",
                                                category,
                                            ]
                                        )

                                    if ecoinvent_version == "uvek":
                                        writer.writerow(
                                            [
                                                a["name"] + "/" + a["location"] + " U",
                                                simapro_units[a["unit"]],
                                                1.0,
                                                "100%",
                                                "not defined",
                                                category,
                                            ])

                        if item == "Materials/fuels":
                            for e in a["exchanges"]:
                                if e["type"] == "technosphere":
                                    if ecoinvent_version in ("3.5", "3.6"):
                                        if not any(i.lower() in e["name"].lower()
                                                   for i in ("waste", "emissions", "treatment", "scrap",
                                                             "used powertrain", "disposal")) \
                                                or any(i in e["name"]
                                                       for i in ["from municipal waste incineration",
                                                                 ]):

                                            if ecoinvent_version == "3.6":
                                                (e["name"], e["location"], e["unit"],
                                                 e["reference product"]) = self.map_37_to_36.get(
                                                    (e["name"], e["location"], e["unit"], e["reference product"]),
                                                    (e["name"], e["location"], e["unit"], e["reference product"])
                                                )
                                            if ecoinvent_version == "3.5":
                                                (e["name"], e["location"], e["unit"],
                                                 e["reference product"]) = self.map_37_to_35.get(
                                                    (e["name"], e["location"], e["unit"], e["reference product"]),
                                                    (e["name"], e["location"], e["unit"], e["reference product"])
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
                                                    name += "| market for " + e["reference product"].lower() + " "
                                                if "market group" in e["name"]:
                                                    name += "| market group for " + e["reference product"].lower() + " "

                                                if "production" in e["name"]:
                                                    if len(e["reference product"].split(", ")) > 1:
                                                        name += ("| " + e["reference product"].split(", ")[
                                                            0] + " production, "
                                                                 + e["reference product"].split(", ")[1] + " ")

                                            writer.writerow(
                                                [
                                                    dict_tech.get(
                                                        (e["name"], e["location"]), name
                                                    ) + "| Cut-off, U",
                                                    simapro_units[e["unit"]],
                                                    "{:.3E}".format(e["amount"]),
                                                    "undefined",
                                                    0,
                                                    0,
                                                    0,
                                                ]
                                            )

                                    if ecoinvent_version == "uvek":
                                        if not any(i.lower() in e["name"].lower()
                                                   for i in ("waste", "emissions", "treatment", "scrap",
                                                             "used powertrain", "disposal", "used passenger car",
                                                             "used electric passenger car")) \
                                                or any(i in e["name"]
                                                       for i in ["from municipal waste incineration",
                                                                 "aluminium scrap, new",
                                                                 "brake wear emissions",
                                                                 "tyre wear emissions",
                                                                 "road wear emissions",
                                                                 "used powertrain from electric passenger car"
                                                                 ]):

                                            if e["name"] not in [i["name"] for i in list_act]:

                                                name = self.map_36_to_uvek_for_simapro[
                                                    e["name"], e["location"], e["unit"], e["reference product"]
                                                ]

                                            else:
                                                name = e["name"] + "/" + e["location"] + " U"

                                            uvek_multiplication_factors = {
                                                "market for heat, from steam, in chemical industry": 1 / 2.257,
                                                "steam production, as energy carrier, in chemical industry": 1 / 2.257,
                                                "market group for natural gas, high pressure": 0.842,
                                                "market for natural gas, high pressure": 0.842,
                                                "market for natural gas, high pressure, vehicle grade": 0.842,
                                                "market for chemical factory": 1 / 12.6e6,
                                                "market for used powertrain from electric passenger car, manual dismantling": -1
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

                                            writer.writerow(
                                                [
                                                    name,
                                                    simapro_units[e["unit"]],
                                                    "{:.3E}".format(e["amount"] * factor),
                                                    "undefined",
                                                    0,
                                                    0,
                                                    0,
                                                ])

                        if item == "Resources":
                            for e in a["exchanges"]:
                                if (
                                        e["type"] == "biosphere"
                                        and e["categories"][0] == "natural resource"
                                ):
                                    writer.writerow(
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
                                        e["type"] == "biosphere"
                                        and e["categories"][0] == "air"
                                ):
                                    writer.writerow(
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
                                if (
                                        e["type"] == "biosphere"
                                        and e["categories"][0] == "water"
                                ):
                                    if e["name"].lower() == "water":
                                        e["unit"] = "kilogram"
                                        e["amount"] /= 1000

                                    writer.writerow(
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
                                        e["type"] == "biosphere"
                                        and e["categories"][0] == "soil"
                                ):
                                    writer.writerow(
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
                                if e["type"] == "technosphere":
                                    if any(i.lower() in e["name"].lower()
                                           for i in (" waste ",
                                                     "emissions",
                                                     "treatment",
                                                     "scrap",
                                                     "used powertrain",
                                                     "used passenger car",
                                                     "used electric passenger car",
                                                     "municipal solid waste",
                                                     "disposal")
                                           ) \
                                            and not any(i.lower() in e["name"].lower()
                                                        for i in ("anaerobic",
                                                                  "cooking",
                                                                  "heat",
                                                                  "manual dismantling"
                                                                  )):

                                        name = ""

                                        if ecoinvent_version in ("3.5", "3.6"):

                                            if ecoinvent_version == "3.6":
                                                (e["name"], e["location"], e["unit"],
                                                 e["reference product"]) = self.map_37_to_36.get(
                                                    (e["name"], e["location"], e["unit"], e["reference product"]),
                                                    (e["name"], e["location"], e["unit"], e["reference product"])
                                                )
                                            if ecoinvent_version == "3.5":
                                                (e["name"], e["location"], e["unit"],
                                                 e["reference product"]) = self.map_37_to_35.get(
                                                    (e["name"], e["location"], e["unit"], e["reference product"]),
                                                    (e["name"], e["location"], e["unit"], e["reference product"])
                                                )

                                            name = dict_tech.get(
                                                (e["name"], e["location"]),
                                                e["name"] + " {" + e["location"] + "}"
                                            )

                                            writer.writerow(
                                                [
                                                    name + "| Cut-off, U",
                                                    simapro_units[e["unit"]],
                                                    "{:.3E}".format(e["amount"] * -1),
                                                    "undefined",
                                                    0,
                                                    0,
                                                    0,
                                                ])

                                        if ecoinvent_version == "uvek":

                                            if not any(i in e["name"].lower()
                                                       for i in [
                                                           "brake wear",
                                                           "tyre wear",
                                                           "road wear",
                                                           "aluminium scrap, new",
                                                           "used powertrain from electric passenger car",
                                                       ]):

                                                uvek_multiplication_factors = {
                                                    "market for manual dismantling of used electric passenger car": 1 / 1200,
                                                    "manual dismantling of used passenger car with internal combustion engine": 1 / 1200,
                                                    "market for manual dismantling of used passenger car with internal combustion engine": 1 / 1200,

                                                }

                                                if e["name"] in uvek_multiplication_factors:
                                                    factor = uvek_multiplication_factors[e["name"]]
                                                else:
                                                    factor = 1

                                                if e["name"] not in [i["name"] for i in list_act]:
                                                    name = self.map_36_to_uvek_for_simapro[
                                                        e["name"], e["location"], e["unit"], e["reference product"]
                                                    ]

                                                else:
                                                    name = e["name"] + "/" + e["location"] + " U"

                                                writer.writerow(
                                                    [
                                                        name,
                                                        simapro_units[e["unit"]],
                                                        "{:.3E}".format(e["amount"] * factor),
                                                        "undefined",
                                                        0,
                                                        0,
                                                        0,
                                                    ])

                        writer.writerow([])

                #System description
                writer.writerow(["System description"])
                writer.writerow([])
                writer.writerow(["Name"])
                writer.writerow(["carculator"])
                writer.writerow([])
                writer.writerow(["Category"])
                writer.writerow(["transport"])
                writer.writerow([])
                writer.writerow(["Description"])
                writer.writerow(["Prospective life cycle assessment model for passenger cars developed by PSI"])
                writer.writerow([])
                writer.writerow(["Cut-off rules"])
                writer.writerow(["All environmentally-relevant flows are included, as far as the authors knowledge permits."
                                 "Also, residual material (e.g., biomass residue) and energy (e.g., waste heat) "
                                 "come free of burden, except for the necessary steps to make it reusable"
                                 " (transport, conditioning, etc.)."
                                 ])
                writer.writerow([])
                writer.writerow(["Energy model"])
                writer.writerow(["The energy consumption of vehicles calculated based on a physics model, including "
                                 "inertia, rolling resistance, aerodynamic drag, road gradient, etc."])
                writer.writerow([])
                writer.writerow(["Transport model"])
                writer.writerow(["Based on Sacchi et al. 2020 (in review)"])
                writer.writerow([])
                writer.writerow(["Allocation rules"])
                writer.writerow(["The system modeling is attributional. In the instance of joint-production, the allocation of "
                                 "burden between co-products is generally based on the relative economic revenue of "
                                 "each product, to align with the underlying database ecoinvent cut-off."])
                writer.writerow(["End"])
                writer.writerow([])

                # Literature reference
                writer.writerow(["Literature reference"])
                writer.writerow([])
                writer.writerow(["Name"])
                writer.writerow(["Sacchi et al. 2020"])
                writer.writerow([])
                writer.writerow(["Documentation link"])
                writer.writerow(["https://www.psi.ch/en/ta/preprint"])
                writer.writerow([])
                writer.writerow(["Comment"])
                writer.writerow(["Pre-print available at: https://www.psi.ch/en/media/57994/download"])
                writer.writerow([])
                writer.writerow(["Category"])
                writer.writerow(["carculator"])
                writer.writerow([])
                writer.writerow(["Description"])
                description = "carculator: an open-source tool for prospective environmental and " \
                              "economic life cycle assessment of vehicles. When, Where and How can battery-electric " \
                              "vehicles help reduce greenhouse gas emissions?\n"
                description += "Romain Sacchi, Christian Bauer and Brian L. Cox\n"
                description += "Submitted to Environmental Science and Technology on November 17th, 2020"

                writer.writerow([description])


            csvFile.close()

        return filepath_export

    def write_lci_to_bw(self, presamples, ecoinvent_compatibility, ecoinvent_version, forbidden_activities):
        """
        Return a LCIImporter object with the inventory as `data` attribute.

        :return: LCIImporter object to be imported in a Brightway2 project
        :rtype: bw2io.base_lci.LCIImporter
        """
        if presamples == True:
            data, array = self.write_lci(
                presamples, ecoinvent_compatibility, ecoinvent_version, forbidden_activities
            )
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return (i, array)
        else:
            data = self.write_lci(
                presamples, ecoinvent_compatibility, ecoinvent_version, forbidden_activities
            )
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return i

    def best_fit_distribution(self, data, bins=200, ax=None):
        import scipy.stats as st
        import warnings
        import pandas as pd

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
                        end
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
        if self.uncertainty_ID[best_distribution.name] == 2:
            mu, std = st.norm.fit(data)
            return [("uncertainty type", 2), ("scale", std), ("loc", mu)]

        # Normal distribution
        if self.uncertainty_ID[best_distribution.name] == 3:
            return [
                ("uncertainty type", 3),
                ("loc", best_params[0]),
                ("scale", best_params[1]),
            ]

        # Uniform distribution
        if self.uncertainty_ID[best_distribution.name] == 4:
            return [
                ("uncertainty type", 4),
                ("minimum", best_params[0]),
                ("maximum", (best_params[0] + best_params[1])),
            ]

        # Triangular distribution
        if self.uncertainty_ID[best_distribution.name] == 5:
            return [
                ("uncertainty type", 5),
                ("loc", best_params[1]),
                ("minimum", np.min(data)),
                ("maximum", np.max(data)),
            ]

        # Gamma distribution
        if self.uncertainty_ID[best_distribution.name] == 9:
            return [
                ("uncertainty type", 9),
                ("shape", best_params[0]),
                ("scale", best_params[2]),
                ("loc", best_params[1]),
            ]

        # Beta distribution
        if self.uncertainty_ID[best_distribution.name] == 10:
            return [
                ("uncertainty type", 10),
                ("loc", best_params[0]),
                ("shape", best_params[1]),
                ("scale", best_params[3]),
            ]

        # Weibull distribution
        if self.uncertainty_ID[best_distribution.name] == 8:
            return [
                ("uncertainty type", 8),
                ("shape", best_params[0]),
                ("loc", best_params[1]),
                ("scale", best_params[2]),
            ]

        # Student's T distribution
        if self.uncertainty_ID[best_distribution.name] == 12:
            return [
                ("uncertainty type", 12),
                ("shape", best_params[0]),
                ("loc", best_params[1]),
                ("scale", best_params[2]),
            ]

        def make_pdf(self, dist, params, size=10000):
            """Generate distributions's Probability Distribution Function """
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

    def fetch_comment(self, name):

        source, description, special_remarks = d_comment[name]

        source_formatted = "Originally published in: " + source
        description_formatted = "This dataset describes " + description

        if len(special_remarks)>0:
            special_remarks_formatted = "it is important to note the following: " + special_remarks
        else:
            special_remarks_formatted = ""

        return list(source_formatted + description_formatted + special_remarks_formatted)

        return
