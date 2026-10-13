from itertools import product

import numpy as np
import xarray as xr
import yaml
from carculator_utils.energy_consumption import EnergyConsumptionModel
from carculator_utils.model import VehicleModel
from carculator_utils.model_run import repeatable_run

from . import DATA_DIR


class CarModel(VehicleModel):
    """
    This class represents the entirety of the vehicles considered, with useful attributes, such as an array that stores
    all the vehicles parameters.

    :ivar array: multi-dimensional numpy-like array that contains parameters' value(s)
    :ivar cycle: name of a driving cycle, or custom driving cycle
    :ivar gradient: series of gradients, for each second of the driving cycle
    :ivar energy_storage: dictionary with selection of battery chemistry for each powertrain

    Example usage:

    1. **Basic Setup and Simulation**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize car input parameters
            cip = CarInputParameters()
            cip.static()  # Static input parameters

            # Fill the xarray from input parameters
            dcts, arr = fill_xarray_from_input_parameters(cip)

            # Create a CarModel instance and run the simulation
            cm = CarModel(arr, cycle="WLTC")
            cm.set_all()

            # Access simulation results
            print(cm.array.sel(parameter="electricity consumption", year=2020).values)

    2. **Setting Battery Capacity Manually**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize parameters
            cip = CarInputParameters()
            cip.static()

            # Define specific battery capacity
            batt_cap = {"capacity": {("BEV", "Medium", 2020): 50}}

            # Fill xarray
            dcts, arr = fill_xarray_from_input_parameters(cip, scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]})

            # Create CarModel with specific battery capacity
            cm = CarModel(arr, cycle="WLTC", energy_storage=batt_cap)
            cm.set_all()

            # Check the battery capacity setting
            print(cm.array.sel(powertrain="BEV", size="Medium", year=2020, parameter="electric energy stored", value=0).values)

    3. **Customizing Battery Chemistry**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize parameters
            cip = CarInputParameters()
            cip.static()

            # the following chemistries are available
            # NMC-111
            # NMC-532
            # NMC-622
            # NMC-811
            # NMC-955
            # LFP
            # NCA
            # Li-S
            # Li-O2
            # SiB

            # by default, the battery chemistry is:
            # NMC-111 for 2000-2015
            # NMC-622 for 2020
            # NMC-811 for 2025
            # NMC-955 for 2030


            # Define battery chemistry
            batt_chem = {"electric": {("BEV", "Medium", 2020): "LFP"}}

            # Fill xarray
            dcts, arr = fill_xarray_from_input_parameters(cip, scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]})

            # Create CarModel with specific battery chemistry
            cm = CarModel(arr, cycle="WLTC", energy_storage=batt_chem)
            cm.set_all()

            # Check the battery chemistry setting
            print(cm.array.sel(powertrain="BEV", size="Medium", year=2020, parameter="battery cell energy density", value=0).values)

    4. **Setting Target Vehicle Mass**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize parameters
            cip = CarInputParameters()
            cip.static()

            # Define a custom target mass
            mass = {("BEV", "Medium", 2020): 2000}

            # Fill xarray
            dcts, arr = fill_xarray_from_input_parameters(cip, scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]})

            # Create CarModel with target mass
            cm = CarModel(arr, cycle="WLTC", target_mass=mass)
            cm.set_all()

            # Check if the target mass was set
            print(cm.array.sel(powertrain="BEV", size="Medium", year=2020, parameter="curb mass", value=0).values)

    5. **Setting Target Driving Range**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize parameters
            cip = CarInputParameters()
            cip.static()

            # Define a custom target range
            range = {("BEV", "Medium", 2020): 100}

            # Fill xarray
            dcts, arr = fill_xarray_from_input_parameters(cip, scope={"size": ["Medium"], "powertrain": ["BEV"], "year": [2020]})

            # Create CarModel with target range
            cm = CarModel(arr, cycle="WLTC", target_range=range)
            cm.set_all()

            # Check if the target range was achieved
            print(cm.array.sel(powertrain="BEV", size="Medium", year=2020, parameter="range", value=0).values)

    6. **Customizing Power Values**

        .. code-block:: python

            from carculator import CarModel, CarInputParameters, fill_xarray_from_input_parameters

            # Initialize parameters
            cip = CarInputParameters()
            cip.static()

            # Define custom power values
            power = {("BEV", "Medium", 2020): 100, ("ICEV-p", "Medium", 2020): 200}

            # Fill xarray
            dcts, arr = fill_xarray_from_input_parameters(cip, scope={"size": ["Medium"], "powertrain": ["BEV", "ICEV-p"], "year": [2020]})

            # Create CarModel with custom power settings
            cm = CarModel(arr, cycle="WLTC", power=power)
            cm.set_all()

            # Check if the power values were set correctly
            print(cm.array.sel(powertrain="BEV", size="Medium", year=2020, parameter="electric power", value=0).values)
            print(cm.array.sel(powertrain="ICEV-p", size="Medium", year=2020, parameter="combustion power", value=0).values)

    """

    @repeatable_run
    def set_all(self):
        """
        This method runs a series of other methods to obtain the tank-to-wheel energy requirement, efficiency
        of the car, costs, etc.

        :meth:`set_component_masses()`, :meth:`set_car_masses()` and :meth:`set_power_parameters()` are interdependent.
        `powertrain_mass` depends on `power`, `curb_mass` is affected by changes in `powertrain_mass`,
        `combustion engine mass` and `electric engine mass`, and `power` is a function of `curb_mass`.
        Sizing converges per vehicle and sample at a relative tolerance of 1e-5.
        With a BEV target range, battery mass and driving mass converge together:
        each step recalculates energy demand and the capacity needed for that range.

        PHEV utility factors and component retention use the constructor options
        ``electric_utility_factor`` and ``drop_hybrids``.

        """

        self.ecm = EnergyConsumptionModel(
            vehicle_type="car",
            vehicle_size=list(self.array.coords["size"].values),
            powertrains=list(self.array.coords["powertrain"].values),
            cycle=self.cycle,
            gradient=self.gradient,
            country=self.country,
        )

        size_for_range = any(
            key[0] == "BEV" and target is not None
            for key, target in (self.target_range or {}).items()
        )
        sizing_state = (
            ["driving mass", "energy battery mass"]
            if size_for_range
            else "driving mass"
        )
        for iteration in self.iterate_sizing(sizing_state, rtol=1e-05):
            self.set_vehicle_masses()
            self.set_power_parameters()
            self.set_component_masses()
            self.set_auxiliaries()
            self.set_power_battery_properties()
            self.set_battery_preferences()
            self.set_battery_properties()
            self.set_energy_stored_properties()
            # Keep the input recuperation assumption fixed. Energy calculations
            # write a cycle-average transmission efficiency back to the array;
            # feeding that output into later iterations would change the inputs.
            if iteration == 1:
                self.set_recuperation()
            self.set_battery_preferences()

            if "FCEV" in self.array.powertrain.values:
                self.set_fuel_cell_power()
                self.set_fuel_cell_mass()

            # if user-provided values are passed,
            # they override the default values
            if "capacity" in self.energy_storage:
                self.override_battery_capacity()

            if size_for_range:
                self.calculate_ttw_energy()
                # Range retains precedence over a capacity override for the
                # same BEV; other vehicles retain their specified capacities.
                self.override_range()

        if not size_for_range:
            self.calculate_ttw_energy()
        self.set_ttw_efficiency()
        self.set_range()

        self.set_share_recuperated_energy()
        self.set_battery_fuel_cell_replacements()
        self.adjust_cost()

        self.set_electric_utility_factor()
        self.set_electricity_consumption()
        self.set_costs()
        self.set_hot_emissions()
        self.set_particulates_emission()
        self.set_noise_emissions()
        self.set_vehicle_masses()
        self.create_PHEV()
        if self.drop_hybrids:
            self.drop_hybrid()

        self.remove_energy_consumption_from_unavailable_vehicles()

    def set_battery_chemistry(self):
        # override default values for batteries
        # if provided by the user
        if "electric" not in self.energy_storage:
            self.energy_storage["electric"] = {}

        default_chemistries = {
            2000: "NMC-111",
            2005: "NMC-111",
            2010: "NMC-111",
            2015: "NMC-111",
            2020: "NMC-622",
            2025: "NMC-811",
            2030: "NMC-955",
        }

        powertrains = self.array.coords["powertrain"].values
        sizes = self.array.coords["size"].values
        years = self.array.coords["year"].values

        for x in product(
            list(powertrains) if isinstance(powertrains, str) else powertrains,
            list(sizes) if isinstance(sizes, str) else sizes,
            list(years) if isinstance(years, int) else years,
        ):
            if x not in self.energy_storage["electric"]:
                if x[-1] in default_chemistries:
                    self.energy_storage["electric"][x] = default_chemistries[x[-1]]
                elif x[-1] < min(default_chemistries.keys()):
                    self.energy_storage["electric"][x] = "NMC-111"
                else:
                    self.energy_storage["electric"][x] = "NMC-955"

        if "origin" not in self.energy_storage:
            self.energy_storage.update({"origin": "CN"})

    def adjust_cost(self) -> None:
        """
        This method adjusts costs of energy storage over time, to correct for the overly optimistic linear
        interpolation between years.

        """

        cost_factor, _ = self._get_cost_factors()
        years = self.array.year

        # Correction of hydrogen tank cost, per kg
        # Correction of fuel cell stack cost, per kW
        if "FCEV" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel tank cost per kg")
            ] = (1.078e58 * np.exp(-6.32e-2 * years) + 3.43e2) * cost_factor

            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell cost per kW")
            ] = (3.15e66 * np.exp(-7.35e-2 * years) + 2.39e1) * cost_factor

        # Correction of energy battery system cost, per kWh
        list_batt = [
            i
            for i in ["BEV", "PHEV-e", "PHEV-c-p", "PHEV-c-d"]
            if i in self.array.powertrain
        ]
        if len(list_batt) > 0:
            self.array.loc[
                dict(powertrain=list_batt, parameter="energy battery cost per kWh")
            ] = (2.75e86 * np.exp(-9.61e-2 * years) + 5.059e1) * cost_factor

        # Correction of power battery system cost, per kW
        list_pwt = [
            i
            for i in [
                "ICEV-p",
                "ICEV-d",
                "ICEV-g",
                "PHEV-c-p",
                "PHEV-c-d",
                "FCEV",
                "HEV-p",
                "HEV-d",
            ]
            if i in self.array.powertrain
        ]

        if len(list_pwt) > 0:
            self.array.loc[
                dict(powertrain=list_pwt, parameter="power battery cost per kW")
            ] = (8.337e40 * np.exp(-4.49e-2 * years) + 11.17) * cost_factor

        # Correction of combustion powertrain cost for ICEV-g
        if "ICEV-g" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="ICEV-g", parameter="combustion powertrain cost per kW")
            ] = np.clip(
                ((5.92e160 * np.exp(-0.1819 * years) + 26.76) * cost_factor),
                None,
                100,
            )

        self.apply_battery_cost_inputs(projected=True)

    def calculate_ttw_energy(self) -> None:
        """
        This method calculates the energy required to operate auxiliary services as well
        as to move the car. The sum is stored under the parameter label "TtW energy" in :attr:`self.array`.

        """

        self.energy = self.ecm.motive_energy_per_km(
            combustion_controls=self.combustion_controls,
            engine_efficiency=self.get_energy_efficiency_override("engine efficiency"),
            transmission_efficiency=self.get_energy_efficiency_override(
                "transmission efficiency"
            ),
            driving_mass=self["driving mass"],
            rr_coef=self["rolling resistance coefficient"],
            drag_coef=self["aerodynamic drag coefficient"],
            frontal_area=self["frontal area"],
            electric_motor_power=self["electric power"],
            engine_power=self["power"],
            combustion_engine_power=self["combustion power"],
            recuperation_efficiency=self["recuperation efficiency"],
            electric_motor_efficiency=self.get_electric_motor_efficiency(),
            aux_power=self["auxiliary power demand"],
            battery_charge_eff=self["battery charge efficiency"],
            battery_discharge_eff=self["battery discharge efficiency"],
            fuel_cell_system_efficiency=self["fuel cell system efficiency"],
        )

        self.energy = self.energy.assign_coords(
            {
                "powertrain": self.array.powertrain,
                "year": self.array.year,
                "size": self.array.coords["size"],
                "value": self.array.coords["value"],
            }
        )

        if self.energy_consumption:
            self.override_ttw_energy()

        distance = self.energy.sel(parameter="velocity").sum(dim="second") / 1000

        self["transmission efficiency"] = (
            np.ma.array(
                self.energy.loc[dict(parameter="transmission efficiency")],
                mask=self.energy.loc[dict(parameter="power load")] == 0,
            )
            .mean(axis=0)
            .T
        )

        self["engine efficiency"] = (
            (
                (self.energy.sel(parameter="motive energy at wheels").sum(dim="second"))
                / (self.energy.sel(parameter="motive energy").sum(dim="second"))
            ).T
            / self["transmission efficiency"]
            / np.where(
                self["fuel cell system efficiency"] > 0,
                self["fuel cell system efficiency"],
                1,
            )
        )

        self["TtW energy"] = (
            self.energy.sel(
                parameter=[
                    "motive energy",
                    "auxiliary energy",
                ]
                + (
                    ["combustion control energy"]
                    if "combustion control energy" in self.energy.parameter
                    else []
                )
            ).sum(dim=["second", "parameter"])
            / distance
        ).T

        self["TtW energy"] += self.get_regeneration_credit()
        self.set_battery_energy_balance()

        self["TtW energy, combustion mode"] = self["TtW energy"] * (
            self["combustion power share"] > 0
        )
        self["TtW energy, electric mode"] = self["TtW energy"] * (
            self["combustion power share"] == 0
        )

        self["auxiliary energy"] = (
            self.energy.sel(parameter="auxiliary energy").sum(dim="second") / distance
        ).T

    def set_vehicle_masses(self) -> None:
        """
        Define ``curb mass``, ``driving mass``, and ``total cargo mass``.

            * `curb mass <https://en.wikipedia.org/wiki/Curb_weight>`__ is the mass of the vehicle and fuel, without people or cargo.
            * ``total cargo mass`` is the mass of the cargo and passengers.
            * ``driving mass`` is the ``curb mass`` plus ``total cargo mass``.

        .. note::
            driving mass = total cargo mass + driving mass

        """

        self["curb mass"] = self["glider base mass"] * (1 - self["lightweighting"])

        curb_mass_includes = [
            "fuel mass",
            "charger mass",
            "converter mass",
            "inverter mass",
            "power distribution unit mass",
            # Updates with set_components_mass
            "combustion engine mass",
            # Updates with set_components_mass
            "electric engine mass",
            # Updates with set_components_mass
            "powertrain mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
            "battery cell mass",
            "battery BoP mass",
            "fuel tank mass",
        ]

        self["curb mass"] += self[curb_mass_includes].sum(axis=2)

        if self.target_mass:
            self.override_vehicle_mass()

        self["total cargo mass"] = (
            self["average passengers"] * self["average passenger mass"]
            + self["cargo mass"]
        )
        self["driving mass"] = self["curb mass"] + self["total cargo mass"]

    def set_electric_utility_factor(self) -> None:
        """Set the electric utility factor according to a sampled values in Germany (ICTT 2022)
        https://theicct.org/wp-content/uploads/2022/06/real-world-phev-use-jun22-1.pdf

        Real-world range in simulation 20 km 30 km 40 km 50 km 60 km 70 km 80 km
        Observed UF for Germany (Sample-size weighted regression ± 2 standard errors)
        Observed UF private (in %) 30±2 41±2 50±3 58±3 65±3 71±3 75±3

        which correlated the share of km driven in electric-mode to
        the capacity of the battery
        (the range that can be driven in battery-depleting mode).

        The argument `uf` is used to override this relation, if needed.
        `uf` must be a ratio between 0 and .75, for each."""

        if "PHEV-e" in self.array.coords["powertrain"].values:
            if self.electric_utility_factor is None:
                self.array.loc[
                    dict(powertrain="PHEV-e", parameter="electric utility factor")
                ] = np.clip(
                    np.interp(
                        self.array.loc[dict(powertrain="PHEV-e", parameter="range")],
                        [0, 20, 30, 40, 50, 60, 70, 80, 100, 120, 200],
                        [0, 0.13, 0.18, 0.23, 0.28, 0.30, 0.35, 0.40, 0.45, 0.5, 0.75],
                    ),
                    0,
                    0.75,
                )
            else:
                for key, val in self.electric_utility_factor.items():
                    if (
                        "PHEV-e" in self.array.powertrain.values
                        and key in self.array.year.values
                    ):
                        self.array.loc[
                            dict(
                                powertrain="PHEV-e",
                                parameter="electric utility factor",
                                year=key,
                            )
                        ] = val

    def set_costs(self) -> None:
        """
        Calculate the different cost types.
        :return:
        """
        self["glider cost"] = (
            self["glider base mass"] * self["glider cost slope"]
            + self["glider cost intercept"]
        )
        self["lightweighting cost"] = (
            self["glider base mass"]
            * self["lightweighting"]
            * self["glider lightweighting cost per kg"]
        )
        self["electric powertrain cost"] = (
            self["electric powertrain cost per kW"] * self["electric power"]
        )
        self["combustion powertrain cost"] = (
            self["combustion power"] * self["combustion powertrain cost per kW"]
        )
        self["fuel cell cost"] = self["fuel cell power"] * self["fuel cell cost per kW"]
        self["power battery cost"] = (
            self["battery power"] * self["power battery cost per kW"]
        )
        self["energy battery cost"] = (
            self["energy battery cost per kWh"] * self.battery_capacity_for_purchase()
        )
        self["fuel tank cost"] = self["fuel tank cost per kg"] * self["fuel mass"]
        # Per km
        self["energy cost"] = self["energy cost per kWh"] * self["TtW energy"] / 3600

        # Retain the fuel-cost convention; electric modes use grid purchases.
        _ = lambda x: np.where(x == 0, 1, x)
        self["energy cost"] /= _(self["battery charge efficiency"])
        self.set_electricity_costs(self["energy cost per kWh"])

        self["component replacement cost"] = (
            self["energy battery cost"] * self["battery lifetime replacements"]
            + self["fuel cell cost"] * self["fuel cell lifetime replacements"]
        )

        with open(DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
            to_markup = yaml.safe_load(stream)["markup"]

        self[to_markup] *= self["markup factor"]

        # calculate costs per km:
        self["lifetime"] = self["lifetime kilometers"] / self["kilometers per year"]

        with open(DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
            purchase_cost_params = yaml.safe_load(stream)["purchase"]

        self["purchase cost"] = self[purchase_cost_params].sum(axis=2)
        # per km
        valid_lifetime = self["lifetime"] > 0
        safe_lifetime = self["lifetime"].where(valid_lifetime, 1)
        safe_kilometers_per_year = self["kilometers per year"].where(
            self["kilometers per year"] > 0, 1
        )
        safe_interest_rate = self["interest rate"].where(self["interest rate"] != 0, 1)
        amortisation_factor = np.where(
            self["interest rate"] == 0,
            np.array(1) / safe_lifetime,
            safe_interest_rate
            + (
                safe_interest_rate
                / ((np.array(1) + safe_interest_rate) ** safe_lifetime - np.array(1))
            ),
        )
        amortisation_factor = np.where(valid_lifetime, amortisation_factor, 0)
        self["amortised purchase cost"] = np.where(
            valid_lifetime,
            self["purchase cost"] * amortisation_factor / safe_kilometers_per_year,
            0,
        )

        # per km
        self["maintenance cost"] = (
            self["maintenance cost per glider cost"]
            * self["glider cost"]
            / self["kilometers per year"]
        )

        # simple assumption that component replacement
        # occurs at half of life.
        self["amortised component replacement cost"] = (
            (
                self["component replacement cost"]
                * np.where(
                    self["interest rate"] == 0,
                    1,
                    (np.array(1) + safe_interest_rate) ** (-safe_lifetime / 2),
                )
            )
            * amortisation_factor
            / safe_kilometers_per_year
        )

        self["total cost per km"] = (
            self["energy cost"]
            + self["amortised purchase cost"]
            + self["maintenance cost"]
            + self["amortised component replacement cost"]
        )

    def remove_energy_consumption_from_unavailable_vehicles(self):
        """Apply the existing technology/year policy to all energy outputs."""
        pwt = self.array.powertrain
        year = self.array.year
        size = self.array.coords["size"]
        available = xr.ones_like(self["TtW energy"], dtype=bool)
        available &= ~(
            pwt.isin(
                [
                    "BEV",
                    "PHEV-e",
                    "PHEV-c-p",
                    "PHEV-c-d",
                    "FCEV",
                    "PHEV-p",
                    "PHEV-d",
                    "HEV-d",
                    "HEV-p",
                ]
            )
            & (year < 2013)
        )
        available &= ~((size == "Micro") & (pwt != "BEV"))
        available &= ~((pwt == "BEV") & (year <= 2010))
        self.mask_energy_outputs(available)
        if "Micro" in size.values:
            self.array.loc[dict(size="Micro")] = self.array.loc[
                dict(size="Micro")
            ].fillna(0)
