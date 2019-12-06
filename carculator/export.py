"""
.. module: export.py

"""
import numpy as np
from pathlib import Path
from inspect import currentframe, getframeinfo
import csv
import xarray as xr
from . import DATA_DIR
from .background_systems import BackgroundSystemModel
import itertools
from bw2io.export.excel import (
    safe_filename,
    xlsxwriter,
    create_valid_worksheet_name,
)

class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices):
        self.array = array
        self.indices = indices

    def write_lci_excel_to_bw(self):
        list_act = []
        for col in range(0, self.array.shape[1]):
            count = 0

            if len(self.indices[col]) > 3:
                activity_name = self.indices[col][0]
                activity_loc = self.indices[col][1]
                activity_unit = self.indices[col][2]
                activity_ref = self.indices[col][3]

                list_exc = []

                for row in range(0, self.array.shape[0]):
                    if self.array[row, col] != 0:
                        if len(self.indices[row]) > 3:
                            input_activity_name = self.indices[row][0]
                            input_activity_loc = self.indices[row][1]
                            input_activity_unit = self.indices[row][2]
                            input_activity_ref = self.indices[row][3]
                            amount = self.array[row, col]

                            if input_activity_ref == activity_ref:
                                list_exc.append(
                                    {
                                        "name": input_activity_name,
                                        "database": 'test',
                                        "amount": amount,
                                        "unit": input_activity_unit,
                                        "type": 'production',
                                        "location": input_activity_loc,
                                        "reference product": input_activity_ref,
                                        "uncertainty type": 0,
                                    }
                                )
                            else:

                                list_exc.append(
                                    {
                                        "name": input_activity_name,
                                        "database": 'test',
                                        "amount": amount * -1,
                                        "unit": input_activity_unit,
                                        "type": 'technosphere',
                                        "location": input_activity_loc,
                                        "reference product": input_activity_ref,
                                        "uncertainty type": 0,
                                    }
                                )
                                count += 1

                        else:
                            input_bio_name = self.indices[row][0]
                            input_bio_cat = self.indices[row][1]
                            input_bio_unit = self.indices[row][2]
                            amount = self.array[row, col]

                            list_exc.append(
                                {
                                    "name": input_bio_name,
                                    "database": 'test_bio',
                                    "amount": amount * -1,
                                    "unit": input_bio_unit,
                                    "type": 'biosphere',
                                    "categories": input_bio_cat,
                                    "uncertainty type": 0,
                                }
                            )
                            count += 1

                if count > 0:
                    list_act.append(
                        {
                            "production amount": 1,
                            "database": 'test',
                            "name": activity_name,
                            "unit": activity_unit,
                            "location": activity_loc,
                            "exchanges": list_exc,
                            "reference product": activity_ref,
                            "type": "process",
                        }
                    )
                    count = 0

        data = []

        data.extend((["Database", 'test'], ("format", "Excel spreadsheet")))
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

        safe_name = safe_filename('test', False)

        filepath = "lci-" + safe_name + ".xlsx"

        workbook = xlsxwriter.Workbook(filepath)
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

        sheet = workbook.add_worksheet(create_valid_worksheet_name('test'))

        for row_index, row in enumerate(data):
            for col_index, value in enumerate(row):
                if value is None:
                    continue
                elif isinstance(value, float):
                    sheet.write_number(row_index, col_index, value, frmt(value))
                else:
                    sheet.write_string(row_index, col_index, value, frmt(value))
        print("Inventories exported to {}.".format(filepath))
        workbook.close()


