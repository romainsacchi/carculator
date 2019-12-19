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
import numpy as np
import pyprind
import os

class ExportInventory:
    """
    Export the inventory to various formats

    """

    def __init__(self, array, indices, db_name = "carculator export"):
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

    def write_lci(self, presamples=True):
        """
        Return the inventory as a dictionary
        If if there several values for one exchange, uncertainty information is generated.
        If `presamples` is True, returns the inventory as well as a `presamples` matrix.
        If `presamples` is False, returns the inventory with characterized uncertainty information.

        :returns: a dictionary that contains all the exchanges
        :rtype: dict
        """
        list_act = []

        if presamples:
            presamples_matrix = []

        # List of coordinates for non-zero values
        non_zeroes = np.nonzero(self.array[0,:,:])
        # List of coordinates where activities present more than once (to filter out "empty" activities, that is,
        # activities with only one reference product exchange)
        u, c = np.unique(non_zeroes[1], return_counts=True)
        dup = u[c > 1]

        # Filter out coordinates of "empty" activities
        coords = np.column_stack((
            non_zeroes[0][np.isin(non_zeroes[1], dup)],
            non_zeroes[1][np.isin(non_zeroes[1], dup)]
        ))

        # Iterate through activities
        bar = pyprind.ProgBar(len(dup))
        for d in dup:
            bar.update(item_id=d)
            list_exc = []
            for row, col in coords[coords[:,1] == d]:
                tuple_output = self.indices[col]
                tuple_input = self.indices[row]

                if len(self.array[:, row, col]) == 1:
                    # No uncertainty, only one value
                    amount = self.array[0, row, col]
                    uncertainty = [('uncertainty type', 1)]

                elif np.all(np.isclose(self.array[:, row, col], self.array[0, row, col])):
                    # Several values, but all the same, so no uncertainty
                    amount = self.array[0, row, col]
                    uncertainty = [('uncertainty type', 1)]
                else:
                    # Uncertainty
                    if presamples == True:
                        # Generate pre-sampled values
                        amount = np.median(self.array[:, row, col])
                        uncertainty = [('uncertainty type', 1)]
                        if len(tuple_input)>3:
                            type_exc = "technosphere"
                        else:
                            type_exc = "biosphere"

                        presamples_matrix.append(
                                                    (
                                                     self.array[:, row, col] * -1,
                                                     [(tuple_input, tuple_output, type_exc)], type_exc
                                                    )
                                                )
                    #else:
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
                                        "type": 'production',
                                        "location": tuple_output[1],
                                        "reference product": tuple_output[3]
                                    }
                        )
                    list_exc[-1].update(uncertainty)
                # If not, if input is technosphere exchange
                elif len(tuple_input)>3:
                    list_exc.append(
                                    {
                                        "name": tuple_input[0],
                                        "database": self.db_name,
                                        "amount": amount * -1,
                                        "unit": tuple_input[2],
                                        "type": 'technosphere',
                                        "location": tuple_input[1],
                                        "reference product": tuple_input[3]
                                    }
                            )
                    list_exc[-1].update(uncertainty)
                # If not, then input is biosphere exchange
                else:
                    list_exc.append(
                                    {
                                        "name": tuple_input[0],
                                        "database": 'biosphere3',
                                        "amount": amount * -1,
                                        "unit": tuple_input[2],
                                        "type": 'biosphere',
                                        "categories": tuple_input[1]
                                    }
                                )
                    list_exc[-1].update(uncertainty)

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
                    "code": str(uuid.uuid1())
                }
            )
        if presamples:
            return (
                list_act,
                presamples_matrix
            )
        else:
            return list_act



    def write_lci_to_excel(self, directory=None):
        """
        Export an Excel file that can be consumed by Brightway2.

        :returns: returns the file path of the exported inventory.
        :rtype: str.
        """

        list_act = self.write_lci(presamples = False)
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

        if directory is None:
            filepath = "lci-" + safe_name + ".xlsx"
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)
            filepath = os.path.join(directory, "lci-" + safe_name + ".xlsx")
            print(filepath)

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

    def write_lci_to_bw(self, presamples):
        """
        Return a LCIImporter object with the inventory as `data` attribute.

        :return: LCIImporter object to be imported in a Brightway2 project
        :rtype: bw2io.base_lci.LCIImporter
        """
        if presamples == True:
            data, array = self.write_lci(presamples)
            i = bw2io.importers.base_lci.LCIImporter(self.db_name)
            i.data = data
            return (i, array)
        else:
            data = self.write_lci(presamples)
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
        #Get histogram of original data
        y, x = np.histogram(data, bins=bins, density=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0

        #Distributions to check
        DISTRIBUTIONS = [
         #st.beta,
         #st.gamma,
         #st.lognorm,
         st.norm,
         #st.t,
         #st.triang,
         #st.uniform,
         #st.weibull_min,
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


        # Lognormal distribution
        if self.uncertainty_ID[best_distribution.name] == 2:
            mu, std = st.norm.fit(data)
            return [
                    ("uncertainty type", 2),
                    ("scale", std),
                    ("loc", mu)
            ]

        # Normal distribution
        if self.uncertainty_ID[best_distribution.name] == 3:
            return [
                    ("uncertainty type", 3),
                    ("loc", best_params[0]),
                    ("scale", best_params[1])
            ]

        # Uniform distribution
        if self.uncertainty_ID[best_distribution.name] == 4:
            return [
                    ("uncertainty type", 4),
                    ("minimum", best_params[0]),
                    ("maximum", (
                                 best_params[0] + best_params[1]))
                    ]

        # Triangular distribution
        if self.uncertainty_ID[best_distribution.name] == 5:
            return [
                    ("uncertainty type", 5),
                    ("loc", best_params[1]),
                    ("minimum", np.min(data)),
                    ("maximum", np.max(data))
                    ]


        # Gamma distribution
        if self.uncertainty_ID[best_distribution.name] == 9:
            return [
                    ("uncertainty type", 9),
                    ("shape", best_params[0]),
                    ("scale", best_params[2]),
                    ("loc", best_params[1])
                    ]

        # Beta distribution
        if self.uncertainty_ID[best_distribution.name] == 10:
            return [
                    ("uncertainty type", 10),
                    ("loc", best_params[0]),
                    ("shape", best_params[1]),
                    ("scale", best_params[3])
            ]

        # Weibull distribution
        if self.uncertainty_ID[best_distribution.name] == 8:
            return [
                    ("uncertainty type", 8),
                    ("shape", best_params[0]),
                    ("loc", best_params[1]),
                    ("scale", best_params[2])
            ]

        # Student's T distribution
        if self.uncertainty_ID[best_distribution.name] == 12:
            return [
                    ("uncertainty type", 12),
                    ("shape", best_params[0]),
                    ("loc", best_params[1]),
                    ("scale", best_params[2])
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



