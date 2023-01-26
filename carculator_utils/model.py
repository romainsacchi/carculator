from itertools import product
from pathlib import Path
from typing import Dict, List, Union

import numexpr as ne
import numpy as np
import xarray as xr
import yaml

from .background_systems import BackgroundSystemModel
from .driving_cycles import detect_vehicle_type
from .energy_consumption import get_default_driving_cycle_name
from .hot_emissions import HotEmissionsModel
from .noise_emissions import NoiseEmissionsModel
from .particulates_emissions import ParticulatesEmissionsModel


def finite(array, mask_value=0):
    return np.where(np.isfinite(array), array, mask_value)


class VehicleModel:

    """
    This class represents the entirety of the vehicles considered,
    with useful attributes, such as an array that stores
    all the vehicles parameters.

    :ivar array: multi-dimensional numpy-like array that contains parameters' value(s)
    :ivar cycle: name of a driving cycle, or custom driving cycle
    :ivar gradient: series of gradients, for each second of the driving cycle
    :ivar energy_storage: dictionary with selection of battery chemistry for each powertrain

    """

    DATA_DIR = Path(__file__).resolve().parent / "data"

    def __init__(
        self,
        array: xr.DataArray,
        country="CH",
        cycle: Union[None, str, np.ndarray] = None,
        gradient: Union[None, np.ndarray] = None,
        energy_storage: Union[None, Dict] = None,
        electric_utility_factor: float = None,
        drop_hybrids: bool = True,
        payload=None,
        energy_target=None,
        energy_consumption: dict = None,
        target_range: dict = None,
        target_mass: dict = None,
        power: dict = None,
        fuel_blend: dict = None,
    ) -> None:

        self.array = array
        self.country = country

        self.vehicle_type = detect_vehicle_type(list(self.array.coords["size"].values))
        self.cycle = (
            cycle
            if isinstance(cycle, str)
            else get_default_driving_cycle_name(self.vehicle_type)
        )

        self.gradient = gradient
        self.energy_storage = energy_storage or {}
        self.energy_target = energy_target or {2025: 0.85, 2030: 0.7, 2050: 0.6}
        self.payload = payload or {}
        self.fuel_blend = fuel_blend

        self.set_battery_chemistry()
        self.set_battery_preferences()
        self.energy = None
        self.electric_utility_factor = electric_utility_factor
        self.drop_hybrids = drop_hybrids
        self.energy_consumption = energy_consumption or None
        # a range to reach can be defined by the user
        self.target_range = target_range
        self.override_range()
        # a curb mass to reach can be defined by the user
        self.target_mass = target_mass
        # overrides the engine/motor power
        self.power = power

        self.bs = BackgroundSystemModel()
        self.fuel_blend = fuel_blend or self.bs.define_fuel_blends(
            self.array.powertrain.values, self.country, self.array.year.values
        )

    def __call__(self, key: Union[str, List]):

        """
        This method fixes a dimension of the `array` attribute given
        a powertrain technology selected.
        Set up this class as a context manager,
        so we can have some nice syntax

        .. code-block:: python

            with class('some powertrain') as cpm:
                cpm['something']. # Will be filtered for the correct powertrain

        On with block exit, this filter is cleared
        https://stackoverflow.com/a/10252925/164864

        :param key: A powertrain type, e.g., "FCEV"
        :return: An instance of `array` filtered after the powertrain selected.

        """
        if isinstance(key, str):
            key = [key]

        self.__cache = self.array
        self.array = self.array.loc[
            dict(powertrain=[k for k in key if k in self.array.powertrain])
        ]
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.array = self.__cache
        del self.__cache

    def __getitem__(self, key: Union[str, List]) -> xr.DataArray:
        """
        Make class['foo'] automatically filter for the parameter 'foo'
        Makes the model code much cleaner

        :param key: Parameter name
        :return: `array` filtered after the parameter selected
        """

        return self.array.loc[dict(parameter=key)]

    def __setitem__(self, key, value):
        self.array.loc[{"parameter": key}] = value

    def set_all(self):
        pass

    def set_battery_chemistry(self):
        pass

    def set_battery_preferences(self):

        l_parameters = [
            p
            for p in [
                "battery cell energy density",
                "battery cell mass share",
                "battery cycle life",
            ]
            if p in self.array.parameter.values
        ]

        for key, val in self.energy_storage["electric"].items():
            pwt, size, year = key

            if (
                (val is not None)
                & (pwt in self.array.powertrain.values)
                & (year in self.array.year.values)
                & (size in self.array["size"].values)
            ):

                cell_params = self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter=[
                            f"{p}, {val.split('-')[0].strip()}" for p in l_parameters
                        ],
                    )
                ]

                self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter=l_parameters,
                    )
                ] = cell_params.values

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

    def drop_hybrid(self) -> None:
        """
        This method drops the powertrains `PHEV-c-p`, `PHEV-c-d` and `PHEV-e` as they were only used to create the
        `PHEV` powertrain.
        :returns: Does not return anything. Modifies ``self.array`` in place.
        """
        self.array = self.array.sel(
            powertrain=[
                pt
                for pt in self.array.coords["powertrain"].values
                if pt not in ["PHEV-e", "PHEV-c-p", "PHEV-c-d"]
            ]
        )

    def set_electricity_consumption(self) -> None:
        """
        This method calculates the total electricity consumption for BEV
        and plugin-hybrid vehicles
        :returns: Does not return anything. Modifies ``self.array`` in place.
        """
        _ = lambda x: np.where(x == 0, 1, x)

        self["electricity consumption"] = (
            self["TtW energy"]
            / _(self["battery charge efficiency"])
            / _(self["charger efficiency"])
            / 3600
            * (self["charger mass"] > 0)
        )

        var = (
            "range"
            if "range" in self.array.coords["parameter"].values
            else (
                "target range"
                if "target range" in self.array.coords["parameter"].values
                else "daily distance"
            )
        )

        self["fuel consumption"] = (
            self["fuel mass"] / _(self[var]) / _(self["fuel density per kg"])
        )

    def override_ttw_energy(self):

        # override of TtW energy, provided by the user
        if self.energy_consumption:
            for key, val in self.energy_consumption.items():
                pwt, size, year = key
                if val is not None:
                    print(
                        f"Overriding TtW energy for {pwt} {size} {year} "
                        f"with {val} kj/km"
                    )

                    self.array.loc[
                        dict(
                            powertrain=pwt, size=size, year=year, parameter="TtW energy"
                        )
                    ] = val

                    self.energy.loc[
                        dict(
                            powertrain=pwt,
                            size=size,
                            year=year,
                            parameter="motive energy",
                        )
                    ] = (
                        val  # kj/km
                        * self.energy.distance  # km
                        / self.energy.shape[-1]  # seconds
                    ).reshape(
                        1, -1
                    )

    def calculate_ttw_energy(self) -> None:
        """
        This method calculates the energy required to operate auxiliary
        services as well as to move the car. The sum is stored under the
        parameter label "TtW energy" in :attr:`self.array`.

        """

        pass

    def set_fuel_cell_mass(self):

        """
        Specific setup for fuel cells, which are mild hybrids.
        Must be called after :meth:`.set_power_parameters`.
        """

        # our basic fuel cell mass is based
        # on a car fuel cell with 800 mW/cm2
        # the cell power density is adapted for truck or bus use
        # it is decreased comparatively to that of a passenger car
        # to reflect increased durability

        _ = lambda x: np.where(x == 0, 1, x)

        self["fuel cell stack mass"] = (
            self["fuel cell power density"]
            * self["fuel cell power"]
            * (800 / _(self["fuel cell power area density"]))
        )
        self["fuel cell ancillary BoP mass"] = (
            self["fuel cell power"] * self["fuel cell ancillary BoP mass per power"]
        )
        self["fuel cell essential BoP mass"] = (
            self["fuel cell power"] * self["fuel cell essential BoP mass per power"]
        )

        if "FCEV" in self.array.powertrain.values:
            self.array.loc[
                dict(parameter="battery power", powertrain="FCEV")
            ] = self.array.loc[dict(parameter="fuel cell power", powertrain="FCEV")] * (
                np.array(1)
                - self.array.loc[
                    dict(parameter="fuel cell power share", powertrain="FCEV")
                ]
            )

            self.array.loc[dict(parameter="battery cell mass", powertrain="FCEV")] = (
                self.array.loc[dict(parameter="battery power", powertrain="FCEV")]
                / self.array.loc[
                    dict(parameter="battery cell power density", powertrain="FCEV")
                ]
            )

            self.array.loc[
                dict(parameter="battery BoP mass", powertrain="FCEV")
            ] = self.array.loc[
                dict(parameter="battery cell mass", powertrain="FCEV")
            ] * (
                np.array(1)
                - self.array.loc[
                    dict(parameter="battery cell mass share", powertrain="FCEV")
                ]
            )

    def set_fuel_cell_power(self) -> None:
        """
        Specific setup for fuel cells, which are mild hybrids.
        Must be called after :meth:`.set_power_parameters`.
        """

        _ = lambda x: np.where(x == 0, 1, x)

        self["fuel cell system efficiency"] = (
            self["fuel cell stack efficiency"]
            / _(self["fuel cell own consumption"])
            * (self["fuel cell own consumption"] > 0)
        )

        self["fuel cell power"] = (
            self["power"]
            * self["fuel cell power share"]
            * self["fuel cell own consumption"]
        )

    def set_auxiliaries(self) -> None:
        """
        Calculates the power needed to operate the auxiliary services
        of the vehicle (heating, cooling).

        The demand for heat and cold are expressed as a fraction of the
        heating and cooling capacities

        .. note:

            Auxiliary power demand (W) = Base auxiliary power (W) +
            (Heating demand (dimensionless, between 0 and 1) * Heating power (W)) +
            (Cooling demand (dimensionless, between 0 and 1) * Cooling power (W))

        """
        self["auxiliary power demand"] = (
            self["auxilliary power base demand"]
            + self["heating thermal demand"] * self["heating energy consumption"]
            + self["cooling thermal demand"] * self["cooling energy consumption"]
        )

    def set_recuperation(self):

        _ = lambda x: np.where(x == 0, 1, x)
        self["recuperation efficiency"] = _(
            self["transmission efficiency"] * (self["combustion power share"] < 1)
        )

    def set_battery_fuel_cell_replacements(self) -> None:
        """
        Calculates the fraction of the replacement battery
        needed to match the vehicle lifetime.

        .. note::
            if ``car lifetime`` = 200000 (km) and
            ``battery lifetime`` = 190000 (km)
            then ``replacement battery`` = 0.05

        .. note::
            It is debatable whether this is realistic or not.
            Car owners may not decide to invest in a new
            battery if the remaining lifetime of the car is
            only 10000 km. Also, a battery lifetime may be expressed
            in other terms, e.g., charging cycles.

        """
        # Here we assume that we can use fractions of a battery
        # (averaged across the fleet)
        self["battery lifetime replacements"] = finite(
            np.clip(
                (self["lifetime kilometers"] / self["battery lifetime kilometers"]) - 1,
                0,
                None,
            )
        )

        # The number of fuel cell replacements is based on the
        # average distance driven with a set of fuel cells given
        # their lifetime expressed in hours of use.
        # The number of replacement is rounded *up* as we assume
        # no allocation of burden with a second life

        average_speed = (
            np.nanmean(
                np.where(
                    self.energy.sel(parameter="velocity") > 0,
                    self.energy.sel(parameter="velocity"),
                    np.nan,
                ),
                0,
            )
            * 3.6
        )

        _ = lambda array: np.where(array == 0, 1, array)

        self["fuel cell lifetime replacements"] = np.ceil(
            np.clip(
                self["lifetime kilometers"]
                / (average_speed.T * _(self["fuel cell lifetime hours"]))
                - 1,
                0,
                5,
            )
        ) * (self["fuel cell lifetime hours"] > 0)

    def override_vehicle_mass(self):

        for key, target_mass in ().items():
            pwt, size, year = key

            if target_mass:
                current_curb_mass = self.array.loc[
                    dict(powertrain=pwt, size=size, year=year, parameter="curb mass")
                ]
                mass_difference = target_mass - current_curb_mass

                lightweighting = self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter="lightweighting",
                    )
                ]

                self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter="glider base mass",
                    )
                ] += mass_difference / (np.array(1) - lightweighting)

        self.set_vehicle_masses()

    def set_vehicle_masses(self) -> None:
        """
        Define ``curb mass``, ``driving mass``, and ``total cargo mass``.

            * `curb mass <https://en.wikipedia.org/wiki/Curb_weight>`__
            is the mass of the vehicle and fuel, without people or cargo.
            * ``total cargo mass`` is the mass of the cargo and passengers.
            * ``driving mass`` is the ``curb mass`` plus ``total cargo mass``.

        .. note::
            driving mass = total cargo mass + driving mass

        """

        pass

    def override_power(self):

        if self.power:
            for key, power in self.power.items():
                pwt, size, year = key
                if power:
                    self.array.loc[
                        dict(powertrain=pwt, size=size, year=year, parameter="power")
                    ] = power

    def set_power_parameters(self) -> None:
        """
        Set electric and combustion motor powers
        based on input parameter ``power to mass ratio``.
        """
        # Convert from W/kg to kW
        self["power"] = self["power to mass ratio"] * self["curb mass"] / 1000

        self["combustion power share"] = self["combustion power share"].clip(
            min=0, max=1
        )
        self["combustion power"] = self["power"] * self["combustion power share"]
        self["electric power"] = self["power"] * (
            np.array(1) - self["combustion power share"]
        )

    def set_component_masses(self) -> None:

        self["combustion engine mass"] = (
            self["combustion power"] * self["combustion mass per power"]
            + self["combustion fixed mass"]
        )

        self["electric engine mass"] = (
            self["electric power"] * self["electric mass per power"]
            + self["electric fixed mass"]
        ) * (self["electric power"] > 0)

        self["powertrain mass"] = (
            self["power"] * self["powertrain mass per power"]
            + self["powertrain fixed mass"]
        )

    def set_share_recuperated_energy(self) -> None:
        """
        Calculate the share of recuperated energy,
        over the total negative motive energy.
        """

        _ = lambda x: np.where(x == 0, 1, x)

        self["share recuperated energy"] = (
            self.energy.sel(parameter="recuperated energy").sum(dim="second")
            / _(self.energy.sel(parameter="negative motive energy").sum(dim="second"))
            * (self["combustion power share"] < 1)
        ).T

        if "PHEV-d" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="PHEV-c-d", parameter="share recuperated energy")
            ] = self.array.loc[
                dict(powertrain="PHEV-e", parameter="share recuperated energy")
            ]

        if "PHEV-p" in self.array.powertrain:
            self.array.loc[
                dict(powertrain="PHEV-c-p", parameter="share recuperated energy")
            ] = self.array.loc[
                dict(powertrain="PHEV-e", parameter="share recuperated energy")
            ]

    def set_electric_utility_factor(self) -> None:
        pass

    def create_PHEV(self):
        """
        PHEV-p/d is the range-weighted average
        between PHEV-c-p/PHEV-c-d and PHEV-e.
        """
        _ = lambda array: np.where(array == 0, 1, array)

        for pwt, pwtc in (("PHEV-d", "PHEV-c-d"), ("PHEV-p", "PHEV-c-p")):
            if pwt in self.array.coords["powertrain"].values:

                self.array.loc[:, pwt] = (
                    self.array.loc[:, "PHEV-e"]
                    * self.array.loc[:, "PHEV-e", "electric utility factor"]
                ) + (
                    self.array.loc[:, pwtc]
                    * (
                        np.array(1)
                        - self.array.loc[:, "PHEV-e", "electric utility factor"]
                    )
                )

                self.array.loc[:, pwt, "electric utility factor"] = self.array.loc[
                    :, "PHEV-e", "electric utility factor"
                ]

                self.energy.loc[
                    dict(
                        powertrain=pwt,
                    )
                ] = self.energy.loc[dict(powertrain="PHEV-e")]

                self.energy.loc[dict(powertrain=pwt,)] *= self.array.loc[
                    dict(parameter="electric utility factor", powertrain="PHEV-e")
                ].T.values[None, ..., None]

                self.energy.loc[dict(powertrain=pwt,)] += (
                    np.array(1)
                    - self.array.loc[
                        dict(parameter="electric utility factor", powertrain="PHEV-e")
                    ].T.values[None, ..., None]
                ) * self.energy.loc[dict(powertrain=pwtc)]

                # We need to preserve the fuel mass and fuel tank mass
                self.array.loc[
                    dict(
                        parameter=[
                            # "fuel mass",
                            "fuel tank mass",
                            "oxidation energy stored",
                            "LHV fuel MJ per kg",
                            "fuel density per kg",
                        ],
                        powertrain=pwt,
                    )
                ] = self.array.loc[
                    dict(
                        parameter=[
                            # "fuel mass",
                            "fuel tank mass",
                            "oxidation energy stored",
                            "LHV fuel MJ per kg",
                            "fuel density per kg",
                        ],
                        powertrain=pwtc,
                    )
                ]

                # we need to preserve the battery mass from PHEV-e as well
                self.array.loc[
                    dict(
                        parameter=[
                            "energy battery mass",
                            "battery BoP mass",
                            "battery cell mass",
                            "battery DoD",
                            "battery cell energy density",
                            "battery charge efficiency",
                            "battery discharge efficiency",
                            "battery lifetime kilometers",
                            "charger efficiency",
                            "recuperation efficiency",
                        ],
                        powertrain=pwt,
                    )
                ] = self.array.loc[
                    dict(
                        parameter=[
                            "energy battery mass",
                            "battery BoP mass",
                            "battery cell mass",
                            "battery DoD",
                            "battery cell energy density",
                            "battery charge efficiency",
                            "battery discharge efficiency",
                            "battery lifetime kilometers",
                            "charger efficiency",
                            "recuperation efficiency",
                        ],
                        powertrain="PHEV-e",
                    )
                ]

                # We store the tank-to-wheel energy consumption
                # in combustion and electric mode separately
                self.array.loc[
                    dict(parameter="TtW energy", powertrain=pwt)
                ] = self.array.loc[dict(parameter="TtW energy", powertrain=pwtc)] * (
                    1
                    - self.array.loc[
                        dict(parameter="electric utility factor", powertrain="PHEV-e")
                    ]
                )
                self.array.loc[
                    dict(parameter="TtW energy", powertrain=pwt)
                ] += self.array.loc[
                    dict(parameter="TtW energy", powertrain="PHEV-e")
                ] * (
                    self.array.loc[
                        dict(parameter="electric utility factor", powertrain="PHEV-e")
                    ]
                )

                self.array.loc[
                    dict(parameter="TtW energy, combustion mode", powertrain=pwt)
                ] = self.array.loc[dict(parameter="TtW energy", powertrain=pwtc)]

                self.array.loc[
                    dict(parameter="TtW energy, electric mode", powertrain=pwt)
                ] = self.array.loc[dict(parameter="TtW energy", powertrain="PHEV-e")]

                # We need to recalculate the range as well

                var = (
                    "range"
                    if "range" in self.array.coords["parameter"].values
                    else (
                        "target range"
                        if "target range" in self.array.coords["parameter"].values
                        else "daily distance"
                    )
                )

                self.array.loc[dict(parameter=var, powertrain=pwt)] = (
                    self.array.loc[
                        dict(parameter="oxidation energy stored", powertrain=pwt)
                    ]
                    * 3600
                    / self.array.loc[
                        dict(parameter="TtW energy, combustion mode", powertrain=pwt)
                    ]
                )

                self.array.loc[dict(parameter=var, powertrain=pwt)] += (
                    self.array.loc[
                        dict(parameter="electric energy stored", powertrain=pwt)
                    ]
                    * 3600
                    / self.array.loc[
                        dict(parameter="TtW energy, electric mode", powertrain="PHEV-e")
                    ]
                )

    def set_battery_properties(self) -> None:
        """
        Calculate mass and power of batteries.
        :return:
        """

        self["battery cell mass"] = (
            self["energy battery mass"] * self["battery cell mass share"]
        )

        self["battery BoP mass"] = self["energy battery mass"] * (
            np.array(1.0) - self["battery cell mass share"]
        )

    def override_battery_capacity(self) -> None:
        """
        Override battery capacity.
        :return:
        """

        for key, val in self.energy_storage["capacity"].items():
            pwt, size, year = key
            if val:
                self.array.loc[
                    dict(
                        parameter="energy battery mass",
                        powertrain=pwt,
                        size=size,
                        year=year,
                    )
                ] = (
                    val
                    / self.array.loc[
                        dict(
                            parameter="battery cell energy density",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ]
                    / self.array.loc[
                        dict(
                            parameter="battery cell mass share",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ]
                )

        self.set_battery_properties()

    def override_range(self):
        """
        Set storage size or range for each powertrain.
        :return:
        """

        if self.target_range:
            for key, val in self.target_range.items():
                pwt, size, year = key

                if pwt == "BEV" and val is not None:
                    battery_DoD = self.array.loc[
                        dict(
                            powertrain=pwt,
                            size=size,
                            year=year,
                            parameter="battery DoD",
                        )
                    ]  # maximum depth of discharge allowed (80%)
                    TtW = self.array.loc[
                        dict(
                            powertrain=pwt, size=size, year=year, parameter="TtW energy"
                        )
                    ]  # kj/km

                    energy_stored = val * (TtW / battery_DoD / 3600)

                    self.array.loc[
                        dict(
                            parameter="energy battery mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        np.array(energy_stored)
                        / self.array.loc[
                            dict(
                                parameter="battery cell energy density",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                        / self.array.loc[
                            dict(
                                parameter="battery cell mass share",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                    )

            self.set_battery_properties()

    def set_range(self) -> None:
        """
        Calculate range autonomy of vehicles
        :return:
        """

        self["range"] = (
            self["fuel mass"] * self["LHV fuel MJ per kg"] * np.array(1000)
        ) / self["TtW energy"]

        self["range"] += (
            self["electric energy stored"]
            * self["battery DoD"]
            * np.array(3600)
            / self["TtW energy"]
        )

    def set_average_lhv(self) -> None:
        """
        Calculate average LHV of fuel.
        :return:
        """

        d_map_fuel = {
            "ICEV-p": "petrol",
            "ICEV-d": "diesel",
            "HEV-d": "diesel",
            "HEV-p": "petrol",
            "PHEV-c-d": "diesel",
            "PHEV-c-p": "petrol",
            "ICEV-g": "cng",
            "FCEV": "hydrogen",
        }

        for pt in [
            pwt
            for pwt in [
                "ICEV-p",
                "ICEV-d",
                "HEV-d",
                "HEV-p",
                "PHEV-c-d",
                "PHEV-c-p",
                "ICEV-g",
                "FCEV",
            ]
            if pwt in self.array.coords["powertrain"].values
        ]:
            # calculate the average LHV based on fuel blend
            fuel_type = d_map_fuel[pt]
            primary_fuel_share = self.fuel_blend[fuel_type]["primary"]["share"]
            primary_fuel_lhv = self.fuel_blend[fuel_type]["primary"]["lhv"]
            primary_fuel_density = self.fuel_blend[fuel_type]["primary"]["density"]
            secondary_fuel_share = self.fuel_blend[fuel_type]["secondary"]["share"]
            secondary_fuel_lhv = self.fuel_blend[fuel_type]["secondary"]["lhv"]
            secondary_fuel_density = self.fuel_blend[fuel_type]["secondary"]["density"]

            self.array.loc[dict(powertrain=pt, parameter="LHV fuel MJ per kg")] = (
                (np.array(primary_fuel_share) * primary_fuel_lhv)
                + (np.array(secondary_fuel_share) * secondary_fuel_lhv)
            ).reshape(1, -1, 1)

            self.array.loc[dict(powertrain=pt, parameter="fuel density per kg")] = (
                (np.array(primary_fuel_share) * primary_fuel_density)
                + (np.array(secondary_fuel_share) * secondary_fuel_density)
            ).reshape(1, -1, 1)

    def set_energy_stored_properties(self) -> None:
        """
        Calculate size and capacity of onboard
        energy storage components.
        :return:
        """

        self.set_average_lhv()
        self["oxidation energy stored"] = (
            self["fuel mass"] * self["LHV fuel MJ per kg"] / 3.6
        )

        self["fuel tank mass"] = (
            self["oxidation energy stored"] * self["fuel tank mass per energy"]
        )

        if "ICEV-g" in self.array.coords["powertrain"].values:
            self["fuel tank mass"] += (
                self["oxidation energy stored"] * self["CNG tank mass slope"]
                + self["CNG tank mass intercept"]
            )

        self["electric energy stored"] = (
            self["battery cell mass"] * self["battery cell energy density"]
        )

    def set_power_battery_properties(self):

        _ = lambda x: np.where(x == 0, 1, x)

        self["battery power"] = self["electric power"] * (
            self["combustion power share"] > 0
        )

        self["battery cell mass"] += (
            self["battery power"]
            / _(self["battery cell power density"])
            * (self["combustion power share"] > 0)
        )

        self["battery BoP mass"] += (
            self["battery cell mass"]
            * (np.array(1) - self["battery cell mass share"])
            * (self["combustion power share"] > 0)
        )

    def set_cargo_mass_and_annual_mileage(self):
        pass

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

        with open(self.DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
            to_markup = yaml.safe_load(stream)["markup"]

        self[to_markup] *= self["markup factor"]

        # calculate costs per km:
        self["lifetime"] = self["lifetime kilometers"] / self["kilometers per year"]

        with open(self.DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
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

    def set_ttw_efficiency(self) -> None:
        """
        Fill in the tank-to-wheel efficiency
        calculated by `calculate_ttw_efficiency`.
        :return:
        """
        _ = lambda array: np.where(array == 0, 1, array)

        if "FCEV" in self.array.coords["powertrain"].values:
            self["TtW efficiency"] = (
                _(self["fuel cell system efficiency"])
                * self["transmission efficiency"]
                * self["engine efficiency"]
            )
        else:
            self["TtW efficiency"] = (
                self["transmission efficiency"] * self["engine efficiency"]
            )

        self["TtW efficiency"] *= np.where(
            self["charger mass"] > 0, self["battery discharge efficiency"], 1
        )

    def set_hot_emissions(self) -> None:
        """
        Calculate hot pollutant emissions based on ``driving cycle``.
        The driving cycle is passed to the :class:`HotEmissionsModel` class
        and :meth:`get_emissions_per_powertrain`
        return emissions per substance per second of driving cycle.
        Those are summed up and divided by
        the distance driven, to obtain emissions, in kg per km.
        :return: Does not return anything. Modifies ``self.array`` in place.
        """

        hem = HotEmissionsModel(
            velocity=self.energy.sel(parameter="velocity"),
            cycle_name=self.cycle,
            vehicle_type=self.vehicle_type,
            powertrains=self.array.coords["powertrain"].values,
            sizes=self.array.coords["size"].values,
        )

        with open(
            self.DATA_DIR / "emission_factors" / "exhaust_flows.yaml", "r"
        ) as stream:
            list_direct_emissions = sorted(yaml.safe_load(stream))

        list_direct_emissions = [
            e + f", {c}"
            for c in ["urban", "suburban", "rural"]
            for e in list_direct_emissions
        ]

        with open(
            self.DATA_DIR / "emission_factors" / "euro_classes.yaml", "r"
        ) as stream:
            euro_classes = yaml.safe_load(stream)[self.vehicle_type]

        list_years = np.clip(
            self.array.coords["year"].values,
            min(euro_classes.keys()),
            max(euro_classes.keys()),
        )

        list_euro_classes = [euro_classes[y] for y in list(list_years)]

        # to calculate emissions and degradation factors
        # we need the vehicle's lifetime, annual mileage
        # as well as its instant fuel consumption

        energy_consumption = self.energy.sel(
            parameter=["motive energy", "auxiliary energy", "recuperated energy"],
            size=self.array.coords["size"].values,
            powertrain=self.array.coords["powertrain"].values,
        ).sum(dim="parameter")

        self.array.loc[
            dict(
                parameter=sorted(list_direct_emissions),
            )
        ] = hem.get_hot_emissions(
            euro_class=list_euro_classes,
            lifetime_km=self["lifetime kilometers"],
            energy_consumption=energy_consumption,
            yearly_km=self["kilometers per year"],
        )

    def set_particulates_emission(self) -> None:
        """
        Calculate the emission of particulates according to
        https://www.eea.europa.eu/ds_resolveuid/6USNA27I4D

        and further disaggregated in:
        https://doi.org/10.1016/j.atmosenv.2020.117886

        for:

        - brake wear
        - tire wear
        - road wear
        - re-suspended road dust

        by considering:

        - vehicle mass
        - driving situation (urban, rural, motorway)

        into the following fractions:

        - PM 2.5
        - PM 10

        Emissions are subdivided in compartments: urban, suburban and rural.

        """

        list_param = [
            "tire wear emissions",
            "brake wear emissions",
            "road wear emissions",
            "road dust emissions",
        ]

        pem = ParticulatesEmissionsModel(
            velocity=self.energy.sel(parameter="velocity"),
            mass=self["driving mass"],
        )

        self[list_param] = pem.get_abrasion_emissions()

        # brake emissions are discounted by
        # the use of regenerative braking
        self["brake wear emissions"] *= np.array(1) - self["share recuperated energy"]

    def set_noise_emissions(self) -> None:
        """
        Calculate noise emissions based on ``driving cycle``.
        The driving cycle is passed to the :class:`NoiseEmissionsModel` class
        and :meth:`get_sound_power_per_compartment`
        returns emissions per compartment type ("rural", "non-urban" and "urban")
        per second of driving cycle.

        Noise emissions are not differentiated by size classes at the moment,
        but only by powertrain "type"
        (e.g., combustion, hybrid and electric)

        :return: Does not return anything. Modifies ``self.array`` in place.
        """
        velocity = self.energy.sel(parameter="velocity")
        nem = NoiseEmissionsModel(velocity, vehicle_type=self.vehicle_type)

        with open(
            self.DATA_DIR / "emission_factors" / "noise_flows.yaml", "r"
        ) as stream:
            list_noise_emissions = yaml.safe_load(stream)

        self.array.loc[
            dict(parameter=list_noise_emissions)
        ] = nem.get_sound_power_per_compartment()

    def calculate_cost_impacts(self, sensitivity=False) -> xr.DataArray:
        """
        This method returns an array with cost values per vehicle-km,
        subdivided into the following groups:

            * Purchase
            * Maintenance
            * Component replacement
            * Energy
            * Total cost of ownership

        :return: A xarray array with cost information per vehicle-km
        """

        list_cost_cat = [
            "purchase",
            "maintenance",
            "component replacement",
            "energy",
            "total",
        ]

        response = self.array.sel(
            parameter=[
                "amortised purchase cost",
                "maintenance cost",
                "amortised component replacement cost",
                "energy cost",
                "total cost per km",
            ],
        )

        response.coords["parameter"] = list_cost_cat

        if not sensitivity:
            return response

        return response / response.sel(value="reference")
