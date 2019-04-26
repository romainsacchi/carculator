from . import DATA_DIR
import pandas as pd
import xarray as xr


def load_excel_parameters(filepath=None, worksheet=None):
    if filepath is None:
        filepath = DATA_DIR / "car_parameters.xlsx"
        worksheet = "Car parameters"
    return pd.read_excel(filepath, sheet_name=worksheet, header=[0])
