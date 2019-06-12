"""
.. module: eda.py

"""

import numpy as np
from .energy_consumption import EnergyConsumptionModel


DEFAULT_MAPPINGS = {
    "electric": {"BEV", "PHEV-c", "PHEV-e", "FCEV"},
    "combustion": {"ICEV-p", "HEV-p", "PHEV-c", "ICEV-g", "ICEV-d"},
    "pure_combustion": {"ICEV-p", "ICEV-g", "ICEV-d"},
    "petrol": {"ICEV-p", "HEV-p", "PHEV-c"},
    "cng": {"ICEV-g"},
    "fuel_cell": {"FCEV"},
    "hybrid": {"PHEV-c", "PHEV-e"},
    "combustion_hybrid": {"PHEV-c"},
    "electric_hybrid": {"PHEV-e"},
    "diesel": {"ICEV-d"},
    "battery": {"BEV"},
}


def finite(array, mask_value=0):
    return np.where(np.isfinite(array), array, mask_value)


class CarModel:

    def __init__(self, array, mappings=None, cycle="WLTC"):
        self.array = array
        self.mappings = mappings or DEFAULT_MAPPINGS
        self.ecm = EnergyConsumptionModel(cycle)

    # Set up this class as a context manager, so we can have some nice syntax
    # e.g. the following:
    # with class('some powertrain') as cpm:
    #     cpm['something']. <- Will be filtered for the correct powertrain
    # On with block exit, this filter is cleared
    # https://stackoverflow.com/a/10252925/164864
    def __call__(self, key):
        self.__cache = self.array
        self.array = self.array.sel(powertrain=key)
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.array = self.__cache
        del self.__cache

    # Make class['foo'] automatically filter for the parameter 'foo'
    # Makes the model code much cleaner
    def __getitem__(self, key):
        return self.array.sel(parameter=key)

    def __setitem__(self, key, value):
        self.array.loc[{"parameter": key}] = value

    # Make it easier/more flexible to filter by powertrain types
    def __getattr__(self, key):
        if key in self.mappings:
            return self.mappings[key]
        else:
            return super().__getattr__(key)

    def set_all(self):
        self.set_recuperation()
        self.set_auxiliaries()
        self.set_component_masses()
        self.set_car_masses()
        self.set_power_parameters()
        self.set_battery_fuel_cell_replacements()
        self.set_energy_stored_properties()
        self.set_battery_properties()
        self.set_fuel_cell_parameters()
        self.set_ttw_efficiency()
        self.calculate_ttw_energy()
        self.set_costs()

    def calculate_ttw_energy(self):
        aux_energy = self.ecm.aux_energy_per_km(self["auxiliary power demand"])

        for pt in self.pure_combustion:
            with self(pt):
                aux_energy.loc[{"powertrain": pt}] /= self['engine efficiency']
        for pt in self.fuel_cell:
            with self(pt):
                aux_energy.loc[{"powertrain": pt}] /= self['fuel cell system efficiency']

        motive_energy = self.ecm.motive_energy_per_km(
            driving_mass=self["driving mass"],
            rr_coef=self["rolling resistance coefficient"],
            drag_coef=self["aerodynamic drag coefficient"],
            frontal_area=self["frontal area"],
            ttw_efficiency=self["TtW efficiency"],
            recuperation_efficiency=self["recuperation efficiency"],
            motor_power=self["electric power"],
        ).sum(axis=-1)

        self.motive_energy = motive_energy

        self["TtW energy"] = aux_energy + motive_energy

    def set_fuel_cell_parameters(self):
        """Specific setup for fuel cells, which are mild hybrids.

        Must be called after ``.set_power_parameters``."""
        for pt in self.fuel_cell:
            with self(pt):
                self["fuel cell system efficiency"] = (
                    self["fuel cell stack efficiency"]
                    / self["fuel cell own consumption"]
                )

                self["fuel cell power"] = (
                    self["power"]
                    * self["fuel cell power share"]
                    * self["fuel cell own consumption"]
                )
                # our basic fuel cell mass is based on a car fuel cell with 800 mW/cm2 and 0.51 kg/kW
                self["fuel cell stack mass"] = (
                    0.51
                    * self["fuel cell power"]
                    * 800
                    / self["fuel cell power area density"]
                )
                self["fuel cell ancillary BoP mass"] = (
                    self["fuel cell power"]
                    * self["fuel cell ancillary BoP mass per power"]
                )
                self["fuel cell essential BoP mass"] = (
                    self["fuel cell power"]
                    * self["fuel cell essential BoP mass per power"]
                )

                self["battery power"] = self["fuel cell power"] * (
                    1 - self["fuel cell power share"]
                )
                self["battery cell mass"] = (
                    self["battery power"] / self["battery cell power density"]
                )
                self["battery BoP mass"] = self["battery cell mass"] * (
                    1 - self["battery cell mass share"]
                )

                self["oxidation energy stored"] = self["fuel mass"] * 120 / 3.6  # kWh
                self["fuel tank mass"] = (
                    self["oxidation energy stored"] * self["H2 tank mass per energy"]
                )

    def set_auxiliaries(self):
        self["auxiliary power demand"] = (
            self["auxilliary power base demand"]
            + self["heating thermal demand"] * self["heating energy consumption"]
            + self["cooling thermal demand"] * self["cooling energy consumption"]
        )

    def set_recuperation(self):
        self["recuperation efficiency"] = (
            self["drivetrain efficiency"] * self["battery charge efficiency"]
        )

    def set_battery_fuel_cell_replacements(self):
        # Here we assume that we can use fractions of a battery/fuel cell
        # (averaged across the fleet)
        self['battery lifetime replacements'] = finite(np.clip(
            (self['lifetime kilometers'] / self['battery lifetime kilometers']) - 1,
            0,
            None
        ))
        self['fuel cell lifetime replacements'] = finite(np.clip(
            (self['lifetime kilometers'] / self['fuel cell lifetime kilometers']) - 1,
            0,
            None
        ))

    def set_car_masses(self):
        """Define ``curb mass``, ``driving mass``, and ``total cargo mass``.

        * `curb mass <https://en.wikipedia.org/wiki/Curb_weight>`__ is the mass of the vehicle and fuel, without people or cargo.
        * ``total cargo mass`` is the mass of the cargo and passengers.
        * ``driving mass`` is the ``curb mass`` plus ``total cargo mass``.

        """
        self["curb mass"] = self["glider base mass"] * (1 - self["lightweighting"])
        curb_mass_includes = [
            "fuel mass",
            "charger mass",
            "converter mass",
            "glider base mass",
            "inverter mass",
            "power distribution unit mass",
            "combustion engine mass",
            "electric engine mass",
            "powertrain mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
            "battery cell mass",
            "battery BoP mass",
            "fuel tank mass",
        ]
        for elem in curb_mass_includes:
            self["curb mass"] += self[elem]

        self["total cargo mass"] = (
            self["average passengers"] * self["average passenger mass"]
            + self["cargo mass"]
        )
        self["driving mass"] = self["curb mass"] + self["total cargo mass"]

    def set_power_parameters(self):
        """Set electric and combustion motor powers based on input parameter ``power to mass ratio``."""
        # Convert from W/kg to kW
        self["power"] = self["power to mass ratio"] * self["curb mass"] / 1000
        self["combustion power"] = self["power"] * self["combustion power share"]
        self["electric power"] = self["power"] * (1 - self["combustion power share"])

    def set_component_masses(self):
        self["combustion engine mass"] = (
            self["combustion power"] * self["combustion mass per power"]
            + self["combustion fixed mass"]
        )
        self["electric engine mass"] = (
            self["electric power"] * self["electric mass per power"]
            + self["electric fixed mass"]
        )
        self["powertrain mass"] = (
            self["power"] * self["powertrain mass per power"]
            + self["powertrain fixed mass"]
        )

    def set_battery_properties(self):
        for pt in self.combustion:
            with self(pt):
                self["battery power"] = self["electric power"]
                self["battery cell mass"] = (
                    self["battery power"] / self["battery cell power density"]
                )
                self["battery BoP mass"] = self["battery cell mass"] * (
                    1 - self["battery cell mass share"]
                )
        for pt in self.electric:
            with self(pt):
                self["battery cell mass"] = (
                    self["energy battery mass"] * self["battery cell mass share"]
                )
                self["battery BoP mass"] = self["energy battery mass"] * (
                    1 - self["battery cell mass share"]
                )

    def set_energy_stored_properties(self):
        for pt in self.petrol:
            with self(pt):
                # Assume 42.4 MJ/kg of gasoline, convert to kWh
                self["oxidation energy stored"] = self["fuel mass"] * 42.4 / 3.6
                self["fuel tank mass"] = (
                    self["oxidation energy stored"] * self["fuel tank mass per energy"]
                )
        for pt in self.diesel:
            with self(pt):
                # Assume 48 MJ/kg of gasoline, convert to kWh
                self["oxidation energy stored"] = self["fuel mass"] * 48 / 3.6
                self["fuel tank mass"] = (
                    self["oxidation energy stored"] * self["fuel tank mass per energy"]
                )
        for pt in self.cng:
            with self(pt):
                # Assume 55.5 MJ/kg of gasoline, convert to kWh
                self["oxidation energy stored"] = self["fuel mass"] * 55.5 / 3.6
                self["fuel tank mass"] = (
                    self["oxidation energy stored"] * self["CNG tank mass slope"]
                    + self["CNG tank mass intercept"]
                )
        for pt in self.battery:
            with self(pt):
                self["electric energy stored"] = (
                    self["battery cell mass"] * self["battery cell energy density"]
                )
        for pt in self.electric_hybrid:
            with self(pt):
                self["electric energy stored"] = (
                    self["battery cell mass"] * self["battery cell energy density"]
                )
                # Assume 42.4 MJ/kg of gasoline
                self["fuel tank mass"] = (
                    self["fuel mass"] * 42.4 / 3.6 * self["fuel tank mass per energy"]
                )

        self["battery cell production electricity"] = (
            self["battery cell production energy"]
            * self["battery cell production energy electricity share"]
        )
        self["battery cell production heat"] = (
            self["battery cell production energy"]
            - self["battery cell production electricity"]
        )

    def set_costs(self):
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
            self["energy battery cost per kWh"]
            * self["battery cell mass"]
            * self["battery cell energy density"]
        )
        self["fuel tank cost"] = self["fuel tank cost per kg"] * self["fuel tank mass"]
        # Per km
        self["energy cost"] = self["energy cost per kWh"] * self["TtW energy"] / 3600

        # For battery, need to divide cost of electricity in battery by efficiency of charging
        for pt in self.battery:
            with self(pt):
                self["energy cost"] /= self["battery charge efficiency"]

        self["component replacement cost"] = (
            self["energy battery cost"] * self["battery lifetime replacements"]
            + self["fuel cell cost"] * self["fuel cell lifetime replacements"]
        )

        # calculate costs per km:
        self["lifetime"] = self["lifetime kilometers"] / self["kilometers per year"]
        i = self["interest rate"]
        amortisation_factor = i + (i / ((1 + i) ** self["lifetime"] - 1))

        purchase_cost_list = [
            "battery onboard charging infrastructure cost",
            "combustion exhaust treatment cost",
            "combustion powertrain cost",
            "electric powertrain cost",
            "energy battery cost",
            "fuel cell cost",
            "fuel tank cost",
            "glider cost",
            "heat pump cost",
            "lightweighting cost",
            "power battery cost",
        ]
        for item in purchase_cost_list:
            self["purchase cost"] += self[item]

        # per km
        self["amortised purchase cost"] = (
            self["purchase cost"] * amortisation_factor / self["kilometers per year"]
        )
        # per km
        self["maintenance cost"] = (
            self["maintenance cost per glider cost"]
            * self["glider cost"]
            / self["kilometers per year"]
        )

        # simple assumption that component replacement occurs at half of life.
        self["amortised component replacement cost"] = (
            self["component replacement cost"]
            * ((1 - self["interest rate"]) ** self["lifetime"] / 2)
            * amortisation_factor
            / self["kilometers per year"]
        )

        to_markup = [
            "combustion powertrain cost",
            "component replacement cost",
            "electric powertrain cost",
            "energy battery cost",
            "fuel cell cost",
            "fuel tank cost",
            "glider cost",
            "lightweighting cost",
            "power battery cost",
        ]
        for item in to_markup:
            self[item] *= self["markup factor"]

        self["total cost per km"] = (
            self["energy cost"]
            + self["amortised purchase cost"]
            + self["maintenance cost"]
            + self["amortised component replacement cost"]
        )

    def set_ttw_efficiency(self):
        _ = lambda array: np.where(array == 0, 1, array)

        self["TtW efficiency"] = (
            _(self["battery discharge efficiency"])
            * _(self["fuel cell system efficiency"])
            * self["drivetrain efficiency"]
            * self["engine efficiency"]
        )
