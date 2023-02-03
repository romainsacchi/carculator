from pathlib import Path

import numpy as np
import pandas as pd

from carculator import *

DATA = Path(__file__, "..").resolve() / "fixtures" / "cars_values.xlsx"
OUTPUT = Path(__file__, "..").resolve() / "fixtures" / "test_model_results.xlsx"
ref = pd.read_excel(DATA, index_col=0)

cip = CarInputParameters()
cip.static()
dcts, arr = fill_xarray_from_input_parameters(cip)
cm = CarModel(arr, cycle="WLTC")
cm.set_all()


def test_model_results():
    list_powertrains = [
        "ICEV-p",
        "ICEV-d",
        "PHEV-p",
        "PHEV-d",
        "BEV",
        "ICEV-g",
        "HEV-d",
        "HEV-p",
    ]
    list_sizes = [
        # "Small",
        # "Lower medium",
        "Medium",
        # "Large"
    ]
    list_years = [
        2020,
        # 2030,
        # 2040,
        # 2050
    ]

    l_res = []

    for pwt in list_powertrains:
        for size in list_sizes:
            for year in list_years:
                for param in cm.array.parameter.values:
                    val = float(
                        cm.array.sel(
                            powertrain=pwt,
                            size=size,
                            year=year,
                            parameter=param,
                            value=0,
                        ).values
                    )

                    try:
                        ref_val = (
                            ref.loc[
                                (ref["powertrain"] == pwt)
                                & (ref["size"] == size)
                                & (ref["parameter"] == param),
                                year,
                            ]
                            .values.astype(float)
                            .item(0)
                        )
                    except:
                        ref_val = 1

                    _ = lambda x: np.where(ref_val == 0, 1, ref_val)
                    diff = val / _(ref_val)
                    l_res.append([pwt, size, year, param, val, ref_val, diff])

    pd.DataFrame(
        l_res,
        columns=["powertrain", "size", "year", "parameter", "val", "ref_val", "diff"],
    ).to_excel(OUTPUT)


def test_setting_batt_cap():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]},
    )

    batt_cap = {
        "capacity": {
            ("BEV", "Medium", 2020): 50,
        }
    }

    cm = CarModel(arr, cycle="WLTC", energy_storage=batt_cap)
    cm.set_all()

    assert (
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="electric energy stored",
            value=0,
        ).values
        == 50
    )


def test_setting_battery_chemistry():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]},
    )

    batt_chem = {
        "electric": {
            ("BEV", "Medium", 2020): "LFP",
        }
    }

    cm = CarModel(arr, cycle="WLTC", energy_storage=batt_chem)
    cm.set_all()

    assert cm.array.sel(
        powertrain="BEV",
        size="Medium",
        year=2020,
        parameter="battery cell energy density",
        value=0,
    ).values == np.array(0.15, dtype=np.float32)


def test_setting_range():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]},
    )

    range = {
        ("BEV", "Medium", 2020): 100,
    }

    cm = CarModel(arr, cycle="WLTC", target_range=range)
    cm.set_all()

    assert (
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="electric energy stored",
            value=0,
        ).values
        <= 25
    )

    assert np.isclose(
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="range",
            value=0,
        ).values,
        np.array(100, dtype=np.float32),
        rtol=0.01,
    )


def test_setting_mass():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]},
    )

    mass = {
        ("BEV", "Medium", 2020): 2000,
    }

    cm = CarModel(arr, cycle="WLTC", target_mass=mass)
    cm.set_all()

    assert (
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="curb mass",
            value=0,
        ).values
        == 2000
    )


def test_setting_ttw_energy():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV", "ICEV-p"], "year": [2020]},
    )

    ttw_energy = {
        ("BEV", "Medium", 2020): 1000,
        ("ICEV-p", "Medium", 2020): 2500,
    }

    cm = CarModel(arr, cycle="WLTC", energy_consumption=ttw_energy)
    cm.set_all()

    assert (
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="TtW energy",
            value=0,
        ).values
        == 1000
    )

    assert np.isclose(
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="electricity consumption",
            value=0,
        ).values,
        (1000 / 3600) * 1.17,
        rtol=0.01,
    )

    assert (
        cm.array.sel(
            powertrain="ICEV-p",
            size="Medium",
            year=2020,
            parameter="TtW energy",
            value=0,
        ).values
        == 2500
    )

    _ = lambda x: np.where(x == 0, 1, x)
    assert np.array_equal(
        cm["fuel consumption"],
        (cm["fuel mass"] / (cm["range"]) / _(cm["fuel density per kg"])),
    )


def test_setting_power():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV", "ICEV-p"], "year": [2020]},
    )

    power = {
        ("BEV", "Medium", 2020): 100,
        ("ICEV-p", "Medium", 2020): 200,
    }

    cm = CarModel(arr, cycle="WLTC", power=power)
    cm.set_all()

    assert (
        cm.array.sel(
            powertrain="BEV",
            size="Medium",
            year=2020,
            parameter="electric power",
            value=0,
        ).values
        == 100
    )

    assert (
        cm.array.sel(
            powertrain="ICEV-p",
            size="Medium",
            year=2020,
            parameter="combustion power",
            value=0,
        ).values
        == 200
    )

    assert (
        cm.array.sel(
            powertrain="ICEV-p",
            size="Medium",
            year=2020,
            parameter="power",
            value=0,
        ).values
        == 200
    )

    assert np.array_equal(
        cm["combustion engine mass"],
        (
            cm["combustion power"] * cm["combustion mass per power"]
            + cm["combustion fixed mass"]
        ),
    )


def test_setting_battery_origin():
    cip = CarInputParameters()
    cip.static()
    dcts, arr = fill_xarray_from_input_parameters(
        cip,
        scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]},
    )

    battery_origin = {"origin": "FR"}

    cm = CarModel(arr, cycle="WLTC", energy_storage=battery_origin)
    cm.set_all()

    assert cm.energy_storage["origin"] == "FR"
