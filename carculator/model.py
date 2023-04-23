from itertools import product

import numpy as np
import yaml
from carculator_utils.energy_consumption import EnergyConsumptionModel
from carculator_utils.model import VehicleModel

from . import DATA_DIR


class CarModel(VehicleModel):
    """
    This class represents the entirety of the vehicles considered, with useful attributes, such as an array that stores
    all the vehicles parameters.

    :ivar array: multi-dimensional numpy-like array that contains parameters' value(s)
    :ivar cycle: name of a driving cycle, or custom driving cycle
    :ivar gradient: series of gradients, for each second of the driving cycle
    :ivar energy_storage: dictionary with selection of battery chemistry for each powertrain

    """

    def set_all(self):
        """
        This method runs a series of other methods to obtain the tank-to-wheel energy requirement, efficiency
        of the car, costs, etc.

        :meth:`set_component_masses()`, :meth:`set_car_masses()` and :meth:`set_power_parameters()` are interdependent.
        `powertrain_mass` depends on `power`, `curb_mass` is affected by changes in `powertrain_mass`,
        `combustion engine mass` and `electric engine mass`, and `power` is a function of `curb_mass`.
        The current solution is to loop through the methods until the increment in driving mass is
        inferior to 0.1%.

        :param drop_hybrids: boolean. True by default. If False, the underlying vehicles used to build plugin-hybrid
                vehicles remain present in the array.
        :param electric_utility_factor: array. If an array is passed, its values are used to override the
                electric utility factor for plugin hybrid vehicles. If not, this factor is calculated using a relation
                described in `set_electric_utility_factor()`

        :returns: Does not return anything. Modifies ``self.array`` in place.

        """

        self.ecm = EnergyConsumptionModel(
            vehicle_type="car",
            vehicle_size=list(self.array.coords["size"].values),
            powertrains=list(self.array.coords["powertrain"].values),
            cycle=self.cycle,
            gradient=self.gradient,
            country=self.country,
        )

        diff = 1.0

        while diff > 0.0001:
            old_driving_mass = self["driving mass"].sum().values
            self.set_vehicle_mass()
            self.set_power_parameters()
            self.set_component_masses()
            self.set_auxiliaries()
            self.set_power_battery_properties()
            self.set_battery_properties()
            self.set_energy_stored_properties()
            self.set_recuperation()

            if "FCEV" in self.array.powertrain.values:
                self.set_fuel_cell_power()
                self.set_fuel_cell_mass()

            # if user-provided values are passed,
            # they override the default values
            if "capacity" in self.energy_storage:
                self.override_battery_capacity()

            diff = (self["driving mass"].sum().values - old_driving_mass) / self[
                "driving mass"
            ].sum()

        self.set_ttw_efficiency()
        self.calculate_ttw_energy()

        self.set_range()

        if self.target_range:
            self.override_range()

        self.set_share_recuperated_energy()
        self.set_battery_fuel_cell_replacements()
        self.adjust_cost()

        self.set_electric_utility_factor()
        self.set_electricity_consumption()
        self.set_costs()
        self.set_hot_emissions()
        self.set_particulates_emission()
        self.set_noise_emissions()
        self.create_PHEV()
        if self.drop_hybrids:
            self.drop_hybrid()

        self.remove_energy_consumption_from_unavailable_vehicles()

    def set_battery_chemistry(self):
        # override default values for batteries
        # if provided by the user
        if "electric" not in self.energy_storage:
            self.energy_storage.update(
                {
                    "electric": {
                        x: "NMC-622"
                        for x in product(
                            ["BEV", "PHEV-e", "HEV-d", "HEV-p"],
                            self.array.coords["size"].values,
                            self.array.year.values,
                        )
                    },
                }
            )
        if "origin" not in self.energy_storage:
            self.energy_storage.update({"origin": "CN"})

    def adjust_cost(self) -> None:
        """
        This method adjusts costs of energy storage over time, to correct for the overly optimistic linear
        interpolation between years.

        """

        n_iterations = self.array.shape[-1]
        n_year = len(self.array.year.values)

        # If uncertainty is not considered, the cost factor equals 1.
        # Otherwise, a variability of +/-30% is added.

        if n_iterations == 1:
            cost_factor = 1
        else:
            if "reference" in self.array.value.values.tolist():
                cost_factor = np.ones((n_iterations, 1))
            else:
                cost_factor = np.random.triangular(0.7, 1, 1.3, (n_iterations, 1))

        # Correction of hydrogen tank cost, per kg
        # Correction of fuel cell stack cost, per kW
        if "FCEV" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel tank cost per kg")
            ] = np.reshape(
                (1.078e58 * np.exp(-6.32e-2 * self.array.year.values) + 3.43e2)
                * cost_factor,
                (1, n_year, n_iterations),
            )

            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel tank cost per kg")
            ] = np.reshape(
                (3.15e66 * np.exp(-7.35e-2 * self.array.year.values) + 2.39e1)
                * cost_factor,
                (1, n_year, n_iterations),
            )

        # Correction of energy battery system cost, per kWh
        list_batt = [
            i
            for i in ["BEV", "PHEV-e", "PHEV-c-p", "PHEV-c-d"]
            if i in self.array.powertrain
        ]
        if len(list_batt) > 0:
            self.array.loc[
                dict(powertrain=list_batt, parameter="energy battery cost per kWh")
            ] = np.reshape(
                (2.75e86 * np.exp(-9.61e-2 * self.array.year.values) + 5.059e1)
                * cost_factor,
                (1, 1, n_year, n_iterations),
            )

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
            ] = np.reshape(
                (8.337e40 * np.exp(-4.49e-2 * self.array.year.values) + 11.17)
                * cost_factor,
                (1, 1, n_year, n_iterations),
            )

        # Correction of combustion powertrain cost for ICEV-g
        if "ICEV-g" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="ICEV-g", parameter="combustion powertrain cost per kW")
            ] = np.clip(
                np.reshape(
                    (5.92e160 * np.exp(-0.1819 * self.array.year.values) + 26.76)
                    * cost_factor,
                    (1, n_year, n_iterations),
                ),
                None,
                100,
            )

    def calculate_ttw_energy(self) -> None:
        """
        This method calculates the energy required to operate auxiliary services as well
        as to move the car. The sum is stored under the parameter label "TtW energy" in :attr:`self.array`.

        """

        self.energy = self.ecm.motive_energy_per_km(
            driving_mass=self["driving mass"],
            rr_coef=self["rolling resistance coefficient"],
            drag_coef=self["aerodynamic drag coefficient"],
            frontal_area=self["frontal area"],
            electric_motor_power=self["electric power"],
            engine_power=self["power"],
            recuperation_efficiency=self["recuperation efficiency"],
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

        self["TtW energy"] = (
            self.energy.sel(
                parameter=["motive energy", "auxiliary energy", "recuperated energy"]
            ).sum(dim=["second", "parameter"])
            / distance
        ).T

        self["engine efficiency"] = (
            np.ma.array(
                self.energy.loc[dict(parameter="engine efficiency")],
                mask=self.energy.loc[dict(parameter="power load")] == 0,
            )
            .mean(axis=0)
            .T
        )

        self["transmission efficiency"] = (
            np.ma.array(
                self.energy.loc[dict(parameter="transmission efficiency")],
                mask=self.energy.loc[dict(parameter="power load")] == 0,
            )
            .mean(axis=0)
            .T
        )

        self["TtW energy, combustion mode"] = self["TtW energy"] * (
            self["combustion power share"] > 0
        )
        self["TtW energy, electric mode"] = self["TtW energy"] * (
            self["combustion power share"] == 0
        )

        self["auxiliary energy"] = (
            self.energy.sel(parameter="auxiliary energy").sum(dim="second") / distance
        ).T

    def set_vehicle_mass(self) -> None:
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
            self["energy battery cost per kWh"] * self["electric energy stored"]
        )
        self["fuel tank cost"] = self["fuel tank cost per kg"] * self["fuel mass"]
        # Per km
        self["energy cost"] = self["energy cost per kWh"] * self["TtW energy"] / 3600

        # For battery, need to divide cost of electricity
        # at battery by efficiency of charging
        # to get costs at the "wall socket".

        _ = lambda x: np.where(x == 0, 1, x)
        self["energy cost"] /= _(self["battery charge efficiency"])

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
        amortisation_factor = self["interest rate"] + (
            self["interest rate"]
            / (
                (np.array(1) + self["interest rate"]) ** self["lifetime kilometers"]
                - np.array(1)
            )
        )
        self["amortised purchase cost"] = (
            self["purchase cost"] * amortisation_factor / self["kilometers per year"]
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
                * (
                    (np.array(1) - self["interest rate"]) ** self["lifetime kilometers"]
                    / 2
                )
            )
            * amortisation_factor
            / self["kilometers per year"]
        )

        self["total cost per km"] = (
            self["energy cost"]
            + self["amortised purchase cost"]
            + self["maintenance cost"]
            + self["amortised component replacement cost"]
        )

    def remove_energy_consumption_from_unavailable_vehicles(self):
        """
        This method sets the energy consumption of vehicles that are not available to zero.
        """

        # we flag cars that have a range inferior to 100 km
        # and also BEVs, PHEVs and FCEVs from before 2013
        self["TtW energy"] = np.where((self["range"] < 100), 0, self["TtW energy"])

        pwts = [
            pt
            for pt in [
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
            if pt in self.array.coords["powertrain"].values
        ]

        years = [y for y in self.array.year.values if y < 2013]

        if years:
            self.array.loc[
                dict(
                    parameter="TtW energy",
                    powertrain=pwts,
                    year=years,
                )
            ] = 0

        # and also Micro cars other than BEVs
        if "Micro" in self.array.coords["size"].values:
            self.array.loc[
                dict(
                    parameter="TtW energy",
                    powertrain=[
                        pt
                        for pt in [
                            "PHEV-e",
                            "PHEV-c-p",
                            "PHEV-c-d",
                            "FCEV",
                            "PHEV-p",
                            "PHEV-d",
                            "HEV-d",
                            "HEV-p",
                            "ICEV-p",
                            "ICEV-d",
                            "ICEV-g",
                        ]
                        if pt in self.array.coords["powertrain"].values
                    ],
                    size="Micro",
                )
            ] = 0

            # replace Nans with zeros
            self.array.loc[dict(size="Micro")] = self.array.loc[
                dict(size="Micro")
            ].fillna(0)

        if "BEV" in self.array.coords["powertrain"].values:
            # set the `TtW energy` of BEV vehicles before 2010 to zero
            self.array.loc[
                dict(
                    powertrain="BEV",
                    year=slice(None, 2010),
                    parameter="TtW energy",
                )
            ] = 0
