from . import DATA_DIR
import pandas as pd
import xarray as xr
import numpy as np


def load_excel_parameters(filepath=None, worksheet=None):
    if filepath is None:
        filepath = DATA_DIR / "car_parameters.xlsx"
        worksheet = "Car parameters"
    return pd.read_excel(filepath, sheet_name=worksheet, header=[0])


def create_builtin_xarray_inputs_static():
    """Create an ``xarray`` labeled array for input parameters.

    Dimensions:

        0: Vehicle size, e.g. "small", "medium". str.
        1: Powertrain, e.g. "ICE-p", "BEV". str.
        2: Parameter name, e.g. "energy battery mass". str.
        3. Triangular distribution values, i.e. "minimum", "mode", "maximum". str.
        4: Scenario label, e.g. "base". str. Optional.
        5. Year, e.g. 2020. int. Optional.

    """
    df = load_excel_parameters()

    def get_unique_labels(series):
        return {
            x.strip()
            for o in df['size'].unique()
            for x in o.split(',')
            if x != 'all' and x.strip()
        }

    sizes = sorted(get_unique_labels(df['size']))
    powertrains = sorted(get_unique_labels(df['powertrain']))
    parameters = sorted(df['parameter'].unique().tolist())
    scenarios = ['low', 'base', 'high']
    years = [2018, 2040]
    triangular = ["minimum", "mode", "maximum"]

    array = xr.DataArray(
        np.zeros((len(sizes), len(powertrains), len(parameters),
                 len(triangular), len(scenarios), len(years))),
        coords=[sizes, powertrains, parameters,
                triangular, scenarios, years],
        dims=['size', 'powertrain', 'parameter',
              'values', 'scenario', 'year']
    )

    compartment = None

    for i in enumerate(len(df)):
        row = df.loc[i]
        if row['Unnamed: 0']:
            compartment = row['Unnamed: 0']
        array.attrs[i] = {
            'compartment': compartment,
            'source': row['Source'],
            'comment': row['comment'],
            'unit': row['unit']
        }
        for s in (sizes if row['size'] == 'all' else )

    return array
