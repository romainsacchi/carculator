from . import DATA_DIR
import pandas as pd
import xarray as xr
import numpy as np


def load_excel_parameters(filepath=None, worksheet=None):
    if filepath is None:
        filepath = DATA_DIR / "car_parameters.xlsx"
        worksheet = "Car parameters"
    return pd.read_excel(filepath, sheet_name=worksheet, header=[0])


def create_builtin_xarray_inputs(iterations=1000):
    """Create an ``xarray`` labeled array for input parameters.

    Dimensions:

        0: Vehicle size, e.g. "small", "medium". str.
        1: Powertrain, e.g. "ICE-p", "BEV". str.
        2: Parameter name, e.g. "energy battery mass". str.
        3. Values (static or stochastic). float.
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

    array = xr.DataArray(
        np.zeros((len(sizes), len(powertrains), len(parameters),
                 iterations, len(scenarios), len(years))),
        coords=[sizes, powertrains, parameters,
                np.arange(iterations), scenarios, years],
        dims=['size', 'powertrain', 'parameter',
              'values', 'scenario', 'year']
    )

    return array
