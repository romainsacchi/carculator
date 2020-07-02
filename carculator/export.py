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


class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices, db_name="carculator export"):
        self.array = array
        self.indices = indices
        self.db_name = db_name
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

        self.map_36_to_35 = {
            ("market for nylon 6", "RoW", "kilogram", "nylon 6"): (
                "market for nylon 6",
                "GLO",
                "kilogram",
                "nylon 6",
            ),
            (
                "market for aluminium oxide, metallurgical",
                "IAI Area, EU27 & EFTA",
                "kilogram",
                "aluminium oxide, metallurgical",
            ): ("market for aluminium oxide", "GLO", "kilogram", "aluminium oxide"),
            (
                "market for flat glass, coated",
                "RER",
                "kilogram",
                "flat glass, coated",
            ): (
                "market for flat glass, coated",
                "GLO",
                "kilogram",
                "flat glass, coated",
            ),
            (
                "market for water, decarbonised",
                "RoW",
                "kilogram",
                "water, decarbonised",
            ): (
                "market for water, decarbonised, at user",
                "GLO",
                "kilogram",
                "water, decarbonised, at user",
            ),
            (
                "market for water, decarbonised",
                "DE",
                "kilogram",
                "water, decarbonised",
            ): (
                "market for water, decarbonised, at user",
                "GLO",
                "kilogram",
                "water, decarbonised, at user",
            ),
            (
                "market for transport, freight, sea, tanker for petroleum",
                "GLO",
                "ton kilometer",
                "transport, freight, sea, tanker for petroleum",
            ): (
                "market for transport, freight, sea, transoceanic tanker",
                "GLO",
                "ton kilometer",
                "transport, freight, sea, transoceanic tanker",
            ),
            (
                "market for transport, freight, sea, tanker for liquid goods other than petroleum and liquefied natural gas",
                "GLO",
                "ton kilometer",
                "transport, freight, sea, tanker for liquid goods other than petroleum and liquefied natural gas",
            ): (
                "market for transport, freight, sea, transoceanic tanker",
                "GLO",
                "ton kilometer",
                "transport, freight, sea, transoceanic tanker",
            ),
            ("market for water, deionised", "CH", "kilogram", "water, deionised"): (
                "market for water, deionised, from tap water, at user",
                "CH",
                "kilogram",
                "water, deionised, from tap water, at user",
            ),
            (
                "market for styrene butadiene rubber (SBR)",
                "RER",
                "kilogram",
                "styrene butadiene rubber (SBR)",
            ): ("latex production", "RER", "kilogram", "latex"),
            (
                "market for water, deionised",
                "Europe without Switzerland",
                "kilogram",
                "water, deionised",
            ): (
                "market for water, deionised, from tap water, at user",
                "Europe without Switzerland",
                "kilogram",
                "water, deionised, from tap water, at user",
            ),
            ("market for water, deionised", "RoW", "kilogram", "water, deionised",): (
                "market for water, deionised, from tap water, at user",
                "RoW",
                "kilogram",
                "water, deionised, from tap water, at user",
            ),
            (
                "market for flat glass, uncoated",
                "RER",
                "kilogram",
                "flat glass, uncoated",
            ): (
                "market for flat glass, uncoated",
                "GLO",
                "kilogram",
                "flat glass, uncoated",
            ),
            ("market for water, ultrapure", "RoW", "kilogram", "water, ultrapure"): (
                "market for water, ultrapure",
                "GLO",
                "kilogram",
                "water, ultrapure",
            ),
            ("market for concrete block", "DE", "kilogram", "concrete block"): (
                "market for concrete block",
                "GLO",
                "kilogram",
                "concrete block",
            ),
        }

        self.map_36_to_uvek = self.load_mapping_36_to_uvek()

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
            name, ref_prod, unit, location, uvek_name, uvek_ref_prod, uvek_unit, uvek_loc = row
            dict_uvek[(name, ref_prod, unit, location)] = (uvek_name, uvek_ref_prod, uvek_unit, uvek_loc)

        return dict_uvek


    def write_lci(self, presamples, ecoinvent_compatibility, ecoinvent_version):
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
        # They should not appear in the exported inventories, otherwise they will be duplicate
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
            "CO2 storage/natural gas, post, 200km pipeline, storage 1000m/2025",
            "drilling, deep borehole/m",
            "Sugar beet cultivation {RER} | sugar beet production Europe | Alloc Rec, U",
            "Refined Waste Cooking Oil {RER} | Refining of waste cooking oil Europe | Alloc Rec, U",
            "Ethanol from forest residues",
            "Ethanol from sugarbeet",
            "pipeline, supercritical CO2/km",
            "Biodiesel from algae",
            "Maize cultivation, drying and storage {RER} | Maize production Europe | Alloc Rec, U",
            "Fischer Tropsch reactor and upgrading plant, construction",
            #"Methanol Synthesis",
            "Walls and foundations, for hydrogen refuelling station",
            "container, with pipes and fittings, for diaphragm compressor",
            "RWGS tank construction",
            "storage module, high pressure, at fuelling station",
            "pumps, carbon dioxide capture process",
            "PEM electrolyzer, Operation and Maintenance",
            "heat exchanger, carbon dioxide capture process",
            "biogas upgrading - sewage sludge - amine scrubbing - best",
            "Hydrogen refuelling station, SMR",
            "Hydrogen, gaseous, 700 bar, from SMR NG w/o CCS, at H2 fuelling station",
            "transformer and rectifier unit, for electrolyzer",
            "PEM electrolyzer, ACDC Converter",
            #"Hydrogen, gaseous, 25 bar, from electrolysis",
            "carbon dioxide, captured from atmosphere",
            "PEM electrolyzer, Balance of Plant",
            #"Hydrogen, gaseous, 700 bar, from electrolysis, at H2 fuelling station",
            "Sabatier reaction methanation unit",
            #"Diesel production, synthetic, Fischer Tropsch process",
            #"Gasoline production, synthetic, from methanol",
            #"Syngas, RWGS, Production",
            "PEM electrolyzer, Stack",
            "hot water tank, carbon dioxide capture process",
            "cooling unit, carbon dioxide capture process",
            #"Methane production, synthetic, from electrochemical methanation",
            "diaphragm compressor module, high pressure",
            "carbon dioxide capture system",
            "Hydrogen dispenser, for gaseous hydrogen",
            "diaphragms, for diaphragm compressor",
            "MTG production facility, construction",
            "Disposal, hydrogen fuelling station",
            "production of 2 wt-% potassium iodide solution",
            "production of nickle-based catalyst for methanation",
            #"Methanol distillation",
            "wiring and tubing, carbon dioxide capture process",
            "control panel, carbon dioxide capture process",
            "adsorption and desorption unit, carbon dioxide capture process",
            "Buffer tank",
            "frequency converter, for diaphragm compressor",
            'Hydrogen, gaseous, 30 bar, from hard coal gasification and reforming, at coal gasification plant'
        ]

        uvek_activities_to_remove = [
            "market for activated carbon, granular",
            "market for iodine",
            "market for manganese sulfate",
            "market for molybdenum trioxide",
            "market for nickel sulfate",
            "market for soda ash, light, crystalline, heptahydrate",
        ]

        uvek_multiplication_factors = {
            "Steam, for chemical processes, at plant": 1/2.257, # 2.257 MJ/kg steam @ ambient pressure
            "Natural gas, from high pressure network (1-5 bar), at service station": 0.842,
            "Disposal, passenger car": 1/1600
        }

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
                if (
                    ecoinvent_compatibility == False
                    and tuple_output[0] in activities_to_be_removed
                ):
                    break

                if ecoinvent_compatibility == False:
                    tuple_output = self.map_ecoinvent_remind.get(
                        tuple_output, tuple_output
                    )
                    tuple_input = self.map_ecoinvent_remind.get(
                        tuple_input, tuple_input
                    )

                #if (ecoinvent_compatibility == False
                #    and tuple_output[0].startswith("fuel supply")
                #    and tuple_input[0].startswith("electricity market")):
                #    continue


                if ecoinvent_compatibility == True:

                    tuple_output = self.map_remind_ecoinvent.get(
                        tuple_output, tuple_output
                    )
                    tuple_input = self.map_remind_ecoinvent.get(
                        tuple_input, tuple_input
                    )

                    if ecoinvent_version == "3.5":
                        tuple_output = self.map_36_to_35.get(tuple_output, tuple_output)
                        tuple_input = self.map_36_to_35.get(tuple_input, tuple_input)

                    if ecoinvent_version == "uvek":

                        tuple_output = self.map_36_to_uvek.get(tuple_output, tuple_output)

                        if tuple_input[0] in uvek_activities_to_remove:
                            continue
                        else:
                            tuple_input = self.map_36_to_uvek.get(tuple_input, tuple_input)

                        #print(tuple_input[0])
                        if tuple_input[0] in uvek_multiplication_factors:
                            mult_factor = uvek_multiplication_factors[tuple_input[0]]

                if len(self.array[:, row, col]) == 1:
                    # No uncertainty, only one value
                    amount = self.array[0, row, col] * mult_factor
                    uncertainty = [("uncertainty type", 1)]

                elif np.all(
                    np.isclose(self.array[:, row, col], self.array[0, row, col])
                ):
                    # Several values, but all the same, so no uncertainty
                    amount = self.array[0, row, col]  * mult_factor
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
                        }
                    )
                    list_exc[-1].update(uncertainty)
            else:
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
        filename=None
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

        list_act = self.write_lci(False, ecoinvent_compatibility, ecoinvent_version)

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
                dict_tech[(name, location)] = simapro_name

            headers = [
                "{CSV separator: Semicolon}",
                "{CSV Format version: 7.0.0}",
                "{Decimal separator: .}",
                "{Date separator: /}",
                "{Short date format: dd/MM/yyyy}",
            ]

            fields = [
                "Process",
                "Category type",
                "Time Period",
                "Geography",
                "Technology",
                "Representativeness",
                "Multiple output allocation",
                "Substitution allocation",
                "Cut off rules",
                "Capital goods",
                "Date",
                "Boundary with nature",
                "Record",
                "Generator",
                "Literature references",
                "External documents",
                "Collection method",
                "Data treatment",
                "Verification",
                "Products",
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
                "megajoule": "mj",
                "unit": "unit",
                "square meter": "m2",
                "kilowatt": "kW",
                "hour": "h",
                "square meter-year": "m2a",
                "meter": "m",
                "vehicle-kilometer": "vkm",
                "meter-year": "ma",
            }

            with open(filepath_export, "w", newline="") as csvFile:
                writer = csv.writer(csvFile, delimiter=";")
                for item in headers:
                    writer.writerow([item])
                writer.writerow([])

                for a in list_act:
                    for item in fields:
                        writer.writerow([item])

                        if item == "Process":
                            name = (
                                a["name"].capitalize()
                                + " {"
                                + a.get("location", "GLO")
                                + "}"
                                + "| Cut-off, U"
                            )
                            writer.writerow([name])

                        if item == "Generator":
                            writer.writerow(["carculator " + str(__version__)])

                        if item == "Geography":
                            writer.writerow([a["location"]])

                        if item == "Time Period":
                            writer.writerow(
                                [
                                    "Between 2010 and 2020. Extrapolated to the selected years."
                                ]
                            )

                        if item == "Date":
                            writer.writerow([str(datetime.date.today())])

                        if item == "Cut off rules":
                            writer.writerow(["100:0 - polluter pays-principle."])

                        if item == "Multiple output allocation":
                            writer.writerow(["No"])

                        if item == "Substitution allocation":
                            writer.writerow(["No"])

                        if item == "Capital goods":
                            writer.writerow(
                                [
                                    "Included when relevant (e.g., factory and machinery.)"
                                ]
                            )

                        if item == "Literature references":
                            writer.writerow(
                                [
                                    "Sacchi, R. et al., 2020, Renewable and Sustainable Energy Reviews (in review), https://www.psi.ch/en/ta/preprint"
                                ]
                            )

                        if item == "External documents":
                            writer.writerow(["https://carculator.psi.ch"])

                        if item == "Collection method":
                            writer.writerow(
                                [
                                    "Modeling and assumptions: https://carculator.readthedocs.io/en/latest/modeling.html"
                                ]
                            )

                        if item == "Verification":
                            writer.writerow(["In review. Susceptible to change."])

                        if item == "Products":
                            for e in a["exchanges"]:
                                if e["type"] == "production":
                                    name = (
                                        e["reference product"].capitalize()
                                        + " {"
                                        + e.get("location", "GLO")
                                        + "}"
                                        + "| Cut-off, U"
                                    )

                                    writer.writerow(
                                        [
                                            dict_tech.get(
                                                (a["name"], a["location"]), name
                                            ),
                                            simapro_units[a["unit"]],
                                            1.0,
                                            "100%",
                                            "not defined",
                                            a["database"],
                                        ]
                                    )

                        if item == "Materials/fuels":
                            for e in a["exchanges"]:
                                if (
                                    e["type"] == "technosphere"
                                    and "waste" not in e["name"]
                                ):
                                    name = (
                                        e["reference product"].capitalize()
                                        + " {"
                                        + e.get("location", "GLO")
                                        + "}"
                                        + "| Cut-off, U"
                                    )

                                    writer.writerow(
                                        [
                                            dict_tech.get(
                                                (e["name"], e["location"]), name
                                            ),
                                            simapro_units[e["unit"]],
                                            e["amount"],
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
                                    writer.writerow(
                                        [
                                            dict_bio.get(e["name"]),
                                            simapro_units[e["unit"]],
                                            e["amount"],
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
                                            simapro_units[e["unit"]],
                                            e["amount"],
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
                                    writer.writerow(
                                        [
                                            dict_bio.get(e["name"], e["name"]),
                                            simapro_units[e["unit"]],
                                            e["amount"],
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
                                            simapro_units[e["unit"]],
                                            e["amount"],
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                        if item == "Final waste flows":
                            for e in a["exchanges"]:
                                if e["type"] == "technosphere" and "waste" in e["name"]:
                                    writer.writerow(
                                        [
                                            dict_bio.get(e["name"], e["name"]),
                                            simapro_units[e["unit"]],
                                            e["amount"],
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )

                        writer.writerow([])

            csvFile.close()

        return filepath_export

    def write_lci_to_bw(self, presamples, ecoinvent_compatibility, ecoinvent_version):
        """
        Return a LCIImporter object with the inventory as `data` attribute.

        :return: LCIImporter object to be imported in a Brightway2 project
        :rtype: bw2io.base_lci.LCIImporter
        """
        if presamples == True:
            data, array = self.write_lci(
                presamples, ecoinvent_compatibility, ecoinvent_version
            )
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return (i, array)
        else:
            data = self.write_lci(
                presamples, ecoinvent_compatibility, ecoinvent_version
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
