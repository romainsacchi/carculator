import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from carculator import *
from carculator.model import CarModel


class TestCarModel(unittest.TestCase):

    DATA = Path(__file__, "..").resolve() / "fixtures" / "cars_values.xlsx"
    OUTPUT = Path(__file__, "..").resolve() / "fixtures" / "test_model_results.xlsx"
    ref = pd.read_excel(DATA, index_col=0)

    def setUp(self):
        cip = CarInputParameters()
        cip.static()
        dcts, arr = fill_xarray_from_input_parameters(cip)
        self.cm = CarModel(arr, cycle="WLTC")
        self.cm.set_all()

    def test_model_results(self):
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
            "Medium",
        ]
        list_years = [
            2020,
        ]

        l_res = []

        for pwt in list_powertrains:
            for size in list_sizes:
                for year in list_years:
                    for param in self.cm.array.parameter.values:
                        val = float(
                            self.cm.array.sel(
                                powertrain=pwt,
                                size=size,
                                year=year,
                                parameter=param,
                                value=0,
                            ).values
                        )

                        try:
                            ref_val = (
                                self.ref.loc[
                                    (self.ref["powertrain"] == pwt)
                                    & (self.ref["size"] == size)
                                    & (self.ref["parameter"] == param),
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
            columns=[
                "powertrain",
                "size",
                "year",
                "parameter",
                "val",
                "ref_val",
                "diff",
            ],
        ).to_excel(self.OUTPUT)

    def test_setting_batt_cap(self):
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

        assert np.isclose(
            cm.array.sel(
                powertrain="BEV",
                size="Medium",
                year=2020,
                parameter="electric energy stored",
                value=0,
            ).values,
            50,
            rtol=0.01,
        )

    def test_setting_battery_chemistry(self):
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

        assert np.isclose(
            cm.array.sel(
                powertrain="BEV",
                size="Medium",
                year=2020,
                parameter="battery cell energy density",
                value=0,
            ).values,
            0.16,
            rtol=0.01,
        )

    def test_setting_range(self):
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

    def test_setting_mass(self):
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

    def test_setting_ttw_energy(self):
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

    def test_setting_power(self):
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

    def test_setting_battery_origin(self):
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

    # New tests start here

    def test_set_all(self):
        # Test that set_all completes without errors
        self.assertTrue(hasattr(self.cm, "ecm"))
        self.assertGreater(len(self.cm.array), 0)

    def test_set_battery_chemistry(self):
        self.cm.set_battery_chemistry()
        # Check if the energy_storage dictionary has been populated with expected values
        self.assertIn("electric", self.cm.energy_storage)
        self.assertEqual(self.cm.energy_storage["origin"], "CN")

    def test_adjust_cost(self):
        self.cm.adjust_cost()
        # Check if the costs have been adjusted as expected
        assert self.cm["total cost per km"].sum() > 0

    def test_set_vehicle_mass(self):
        # Mock necessary attributes for testing
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="glider base mass",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([800, 900, 1000]).T
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="lightweighting",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([0.1, 0.1, 0.1])
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="fuel mass",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([50, 60, 70])
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="average passengers",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([2, 2, 2])
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="average passenger mass",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([70, 70, 70])
        self.cm.array.loc[
            dict(
                powertrain=["ICEV-p", "ICEV-d", "ICEV-g"],
                parameter="cargo mass",
                size="Medium",
                year=2020,
                value=0,
            )
        ] = np.array([100, 120, 140])

        # Perform the mass calculation
        self.cm.set_vehicle_mass()

        # Verify the mass has been set correctly
        self.assertIn("driving mass", self.cm.array.coords["parameter"].values)
        self.assertGreater(self.cm["driving mass"].sum(), 0)

    def test_set_electric_utility_factor(self):

        # Check if the electric utility factor is set within expected limits
        self.assertGreaterEqual(self.cm["electric utility factor"].min(), 0)
        self.assertLessEqual(self.cm["electric utility factor"].max(), 0.75)

    def test_remove_energy_consumption_from_unavailable_vehicles(self):

        self.cm.remove_energy_consumption_from_unavailable_vehicles()

        # Check that energy consumption is set to 0 for vehicles that should be unavailable
        self.assertEqual(
            self.cm["TtW energy"].sel(year=2010, powertrain="BEV").sum(), 0
        )


TestCarModel()
