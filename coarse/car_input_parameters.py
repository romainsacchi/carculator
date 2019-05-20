from . import DATA_DIR
from klausen import NamedParameters
from stats_arrays import uncertainty_choices
import itertools
import numpy as np
import pandas as pd


UNCERTAINTY_MAPPING = {
    dist.description.replace(" uncertainty", "").lower().strip(): dist.id
    for dist in uncertainty_choices
}


def load_excel_parameters(filepath=None, worksheet=None):
    if filepath is None:
        filepath = DATA_DIR / "car_parameters.xlsx"
        worksheet = "Car parameters"
    return pd.read_excel(filepath, sheet_name=worksheet, header=[0])


def to_float(x):
    if "%" in x:
        return float(x.replace("%", "") / 100)
    return float(x)


class CarInputParameters(NamedParameters):
    def __init__(self):
        """Create a `klausen <https://github.com/cmutel/klausen>`__ model with the car input parameters.

        The parameter names are not unique, so we append a simple counter to their labels."""
        super().__init__(None)
        self.df = load_excel_parameters()

        self.sizes = self.get_unique_labels(self.df['size'])
        self.powertrains = self.get_unique_labels(self.df['powertrain'])
        # self.parameters = sorted(self.df['parameter'].unique().tolist())
        self.years = [2017, 2040]

        self.add_car_parameters()

    def get_unique_labels(self, series):
        return sorted({
            x.strip()
            for o in series.unique()
            for x in o.split(',')
            if x != 'all' and x.strip()
        })

    def add_car_parameters(self):
        params = {}
        count = itertools.count()

        for i, row in self.df.iterrows():
            for year in self.years:
                if np.isnan(row["{} base".format(year)]):
                    continue

                label = "{}-{}-{}".format(next(count), year, row['parameter'])
                params[label] = {
                    'metadata': {
                        'unit': row['unit'],
                        'source': row['source'],
                        'comment': row['comment'],
                        'sizes': self.sizes if row['size'] == 'all' else row['size'].split(", "),
                        'powertrain': self.powertrains if row['powertrain'] == 'all' else row['powertrain'].split(", "),
                        'category': row['category'],
                    },
                    'kind': 'distribution',
                    'uncertainty_type': UNCERTAINTY_MAPPING[row['uncertainty distribution']],
                    'amount': row["{} base".format(year)],
                    'loc': row["{} base".format(year)],
                    'minimum': row["{} low".format(year)],
                    'maximum': row["{} high".format(year)],
                }

        self.add_parameters(params)
