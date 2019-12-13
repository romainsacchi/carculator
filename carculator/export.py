"""
.. module: export.py

"""
from bw2io.export.excel import (
    safe_filename,
    xlsxwriter,
    create_valid_worksheet_name,
)
import bw2io
import uuid
import numpy

class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices, db_name = "carculator export"):
        self.array = array
        self.indices = indices
        self.db_name = db_name

    def write_lci(self):
        """
        Return the inventory as a dictionary

        :return: a dictionary that contains all the exchanges
        :rtype: dict
        """
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
                                        "database": self.db_name,
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
                                        "database": self.db_name,
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
                                    "database": 'biosphere3',
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
                            "database": self.db_name,
                            "name": activity_name,
                            "unit": activity_unit,
                            "location": activity_loc,
                            "exchanges": list_exc,
                            "reference product": activity_ref,
                            "type": "process",
                            "code": str(uuid.uuid1())
                        }
                    )
                    count = 0
        return list_act

    def write_lci_to_excel(self):
        """
        Export an Excel file that can be consumed by Brightway2.

        :return: returns the file path of the exported inventory
        :rtype: str
        """

        list_act = self.write_lci()
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
        return filepath

    def write_lci_to_bw(self):
        """
        Return a LCIImporter object with the inventory as `data` attribute.
        :return: LCIImporter object to be imported in a Brightway2 project
        :rtype: bw2io.base_lci.LCIImporter
        """

        data = self.write_lci()
        i = bw2io.importers.base_lci.LCIImporter(self.db_name)
        i.data = data
        return i

    #Create models from data
    def best_fit_distribution(self, data, bins=200, ax=None):
     import scipy.stats as st
     import warnings
     import pandas as pd

     """Model data by finding best fit distribution to data"""
      #Get histogram of original data
     y, x = np.histogram(data, bins=bins, density=True)
     x = (x + np.roll(x, -1))[:-1] / 2.0

      #Distributions to check
     DISTRIBUTIONS = [
         st.beta,
         st.gamma,
         st.lognorm,
         st.norm,
         st.t,
         st.triang,
         st.uniform,
         st.weibull_min,
     ]

      #Best holders
     best_distribution = st.norm
     best_params = (0.0, 1.0)
     best_sse = np.inf

      #Estimate distribution parameters from data
     for distribution in DISTRIBUTIONS:

          #Try to fit the distribution
         try:
              #Ignore warnings from data that can't be fit
             with warnings.catch_warnings():
                 warnings.filterwarnings("ignore")

                  #fit dist to data
                 params = distribution.fit(data)

                  #Separate parts of parameters
                 arg = params[:-2]
                 loc = params[-2]
                 scale = params[-1]

                  #Calculate fitted PDF and error with fit in distribution
                 pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                 sse = np.sum(np.power(y - pdf, 2.0))

                 # if axis pass in add to plot
                 try:
                     if ax:
                         pd.Series(pdf, x).plot(ax=ax)
                     end
                 except Exception:
                     pass

                  #identify if this distribution is better
                 if best_sse > sse > 0:
                     best_distribution = distribution
                     best_params = params
                     best_sse = sse

         except Exception:
             pass

     return (
         best_distribution.name,
         getattr(st, best_distribution.name),
         best_params,
     )

    def make_pdf(self, dist, params, size=10000):
         """Generate distributions's Probability Distribution Function """
         import pandas as pd

         # Separate parts of parameters
         arg = params[:-2]
         loc = params[-2]
         scale = params[-1]

         # Get sane start and end points of distribution
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



