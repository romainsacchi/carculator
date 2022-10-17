from itertools import product
from typing import Dict, List, Union

import numexpr as ne
import numpy as np
import xarray as xr
import yaml

from . import DATA_DIR
from .energy_consumption import EnergyConsumptionModel
from .hot_emissions import HotEmissionsModel
from .noise_emissions import NoiseEmissionsModel
from .particulates_emissions import ParticulatesEmissionsModel


def finite(array, mask_value=0):
    return np.where(np.isfinite(array), array, mask_value)


class CarModel:

    """
    This class represents the entirety of the vehicles considered, with useful attributes, such as an array that stores
    all the vehicles parameters.

    :ivar array: multi-dimensional numpy-like array that contains parameters' value(s)
    :ivar cycle: name of a driving cycle, or custom driving cycle
    :ivar gradient: series of gradients, for each second of the driving cycle
    :ivar energy_storage: dictionary with selection of battery chemistry for each powertrain

    """

    def __init__(
        self,
        array: xr.DataArray,
        country="RER",
        cycle: Union[None, str, np.ndarray] = None,
        gradient: Union[None, np.ndarray] = None,
        energy_storage: Union[None, Dict] = None,
        electric_utility_factor: float = None,
        drop_hybrids: bool = True,
        energy_consumption: dict = None,
        target_range: dict = None,
        target_mass: dict = None,
        power: dict = None,
    ) -> None:

        self.array = array
        self.country = country

        if cycle is None:
            self.ecm = EnergyConsumptionModel("WLTC")
        else:
            self.ecm = EnergyConsumptionModel(cycle=cycle, gradient=gradient)

        # override default values for batteries
        # if provided by the user
        self.energy_storage = {
            "electric": {
                x: "NMC-622"
                for x in product(
                    ["BEV", "PHEV-e", "HEV-d", "HEV-p"],
                    self.array.coords["size"].values,
                    self.array.year.values,
                )
            }
        }
        self.energy_storage["origin"] = "CN"
        self.energy_storage.update(energy_storage or {})
        self.set_battery_preferences()
        self.energy = None
        self.electric_utility_factor = electric_utility_factor
        self.drop_hybrids = drop_hybrids
        self.energy_consumption = energy_consumption or None
        # a range to reach cna be defined by the
        self.target_range = target_range
        # a curb mass to reach can be defined by the user
        self.target_mass = target_mass
        # overrides the engine/motor power
        self.power = power

    def __call__(self, key: Union[str, List]):

        """
        This method fixes a dimension of the `array` attribute given a powertrain technology selected.
        Set up this class as a context manager, so we can have some nice syntax

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

        diff = 1.0

        while diff > 0.001:
            old_driving_mass = self["driving mass"].sum().values

            self.set_car_masses()
            self.set_power_parameters()
            self.set_component_masses()
            self.set_battery_properties()
            self.set_battery_fuel_cell_replacements()
            self.set_recuperation()
            self.set_fuel_cell_parameters()
            self.set_energy_stored_properties()

            diff = (self["driving mass"].sum().values - old_driving_mass) / self[
                "driving mass"
            ].sum()

        self.set_auxiliaries()
        self.set_ttw_efficiency()
        self.calculate_ttw_energy()
        self.set_share_recuperated_energy()
        self.adjust_cost()

        # if user-provided values are passed,
        # they override the default values
        if self.target_range or self.energy_storage:
            self.set_storage_size_or_range()
            self.set_energy_stored_properties()

        self.set_range()
        self.set_electric_utility_factor()
        self.set_electricity_consumption()
        self.set_costs()
        self.set_hot_emissions()
        self.set_particulates_emission()
        self.set_noise_emissions()
        self.create_PHEV()
        if self.drop_hybrids:
            self.drop_hybrid()

        # we flag cars that have a range inferior to 100 km
        # and also BEVs, PHEVs and FCEVs from before 2013
        self.array.loc[dict(parameter="has_low_range")] = (
            self.array.loc[dict(parameter="range")] < 100
        )

        self.array.loc[
            dict(
                parameter="has_low_range",
                powertrain=[
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
                ],
                year=[y for y in self.array.year.values if y < 2013],
            )
        ] = 1

        # and also Micro cars other than BEVs
        if "Micro" in self.array.coords["size"].values:
            self.array.loc[
                dict(
                    parameter="has_low_range",
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
            ] = 1

        # fill in NaNs due to `Micro` only existing for `BEV` powertrain
        if "Micro" in self.array.coords["size"].values:
            if "BEV" in self.array.powertrain:
                self.array.loc[dict(size="Micro")] = self.array.loc[
                    dict(size="Micro")
                ].fillna(1)
                self.array.loc[
                    dict(size="Micro", parameter="lifetime kilometers")
                ] = self.array.loc[
                    dict(
                        size="Micro", parameter="lifetime kilometers", powertrain="BEV"
                    )
                ]
                self.array.loc[dict(size="Micro")] = self.array.loc[
                    dict(size="Micro")
                ].fillna(1)
                self.array.loc[
                    dict(size="Micro", parameter="kilometers per year")
                ] = self.array.loc[
                    dict(
                        size="Micro", parameter="kilometers per year", powertrain="BEV"
                    )
                ]

    def set_battery_preferences(self):

        for key, val in self.energy_storage["electric"].items():
            pwt, size, year = key

            if val is not None and pwt in self.array.powertrain.values:

                cell_params = self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter=[
                            f"battery cell energy density, {val.split('-')[0].strip()}",
                            f"battery cell mass share, {val.split('-')[0].strip()}",
                        ],
                    )
                ]

                self.array.loc[
                    dict(
                        powertrain=pwt,
                        size=size,
                        year=year,
                        parameter=[
                            "battery cell energy density",
                            "battery cell mass share",
                        ],
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

    def adjust_fuel_mass(self) -> None:
        """
        This method adjusts the fuel mass over the years, to correct for the linear
        interpolation between years.

        """

        n_iterations = self.array.shape[-1]
        n_year = len(self.array.year.values)

        # Correction of hydrogen mass
        self.array.loc[dict(powertrain="FCEV", parameter="fuel mass")] = np.reshape(
            (1.078e58 * np.exp(-6.32e-2 * self.array.year.values) + 3.43e2),
            (1, 1, n_year, n_iterations),
        )

        # Correction of CNG mass

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
        This method calculates the total electricity consumption for BEV and plugin-hybrid vehicles
        :returns: Does not return anything. Modifies ``self.array`` in place.
        """

        self["electricity consumption"] = (
            self["TtW energy"]
            / self["battery charge efficiency"]
            / 3600
            * (self["charger mass"] > 0)
        )

    def calculate_ttw_energy(self) -> None:
        """
        This method calculates the energy required to operate auxiliary services as well
        as to move the car. The sum is stored under the parameter label "TtW energy" in :attr:`self.array`.

        """

        self.energy = xr.DataArray(
            np.zeros(
                (
                    len(self.array.coords["size"]),
                    len(self.array.coords["powertrain"]),
                    4,
                    len(self.array.coords["year"]),
                    len(self.array.coords["value"]),
                    self.ecm.cycle.shape[0],
                )
            ).astype("float32"),
            coords=[
                self.array.coords["size"],
                self.array.coords["powertrain"],
                [
                    "auxiliary energy",
                    "motive energy",
                    "recuperated energy",
                    "negative motive energy",
                ],
                self.array.coords["year"],
                self.array.coords["value"],
                np.arange(self.ecm.cycle.shape[0]),
            ],
            dims=["size", "powertrain", "parameter", "year", "value", "second"],
        )

        self.energy.loc[dict(parameter="auxiliary energy")] = (
            self["auxiliary power demand"].values[..., None] / 1000
        )

        self.energy.loc[
            dict(
                parameter="auxiliary energy",
                powertrain=[
                    pt
                    for pt in [
                        "ICEV-d",
                        "ICEV-p",
                        "ICEV-g",
                        "PHEV-c-d",
                        "PHEV-c-p",
                        "HEV-p",
                        "HEV-d",
                    ]
                    if pt in self.array.coords["powertrain"].values
                ],
            )
        ] /= self.array.loc[
            dict(
                parameter="engine efficiency",
                powertrain=[
                    pt
                    for pt in [
                        "ICEV-d",
                        "ICEV-p",
                        "ICEV-g",
                        "PHEV-c-d",
                        "PHEV-c-p",
                        "HEV-p",
                        "HEV-d",
                    ]
                    if pt in self.array.coords["powertrain"].values
                ],
            )
        ].values[
            ..., None
        ]

        self.energy.loc[
            dict(
                parameter="auxiliary energy",
                powertrain=[
                    pt
                    for pt in ["FCEV"]
                    if pt in self.array.coords["powertrain"].values
                ],
            )
        ] /= self.array.sel(
            parameter="fuel cell system efficiency",
            powertrain=[
                pt for pt in ["FCEV"] if pt in self.array.coords["powertrain"].values
            ],
        ).values[
            ..., None
        ]

        motive_energy, recuperated_energy, distance = self.ecm.motive_energy_per_km(
            driving_mass=self["driving mass"],
            rr_coef=self["rolling resistance coefficient"],
            drag_coef=self["aerodynamic drag coefficient"],
            frontal_area=self["frontal area"],
            sizes=self.array.coords["size"].values,
            motor_power=self["electric power"],
        )

        self.energy.loc[dict(parameter="motive energy")] = np.clip(
            motive_energy / 1000, 0, None
        )

        self.energy.loc[dict(parameter="negative motive energy")] = (
            np.clip(motive_energy / 1000, None, 0) * -1
        )

        self.energy.loc[dict(parameter="motive energy")] /= self["TtW efficiency"]

        self.energy.loc[dict(parameter="motive energy")] = np.clip(
            self.energy.loc[dict(parameter="motive energy")],
            0,
            self["power"].values[..., None],
        )

        self.energy.loc[dict(parameter="recuperated energy")] = np.clip(
            recuperated_energy / 1000,
            self["power"].values[..., None] * -1,
            0,
        )

        self.energy.loc[dict(parameter="recuperated energy")] *= self[
            "recuperation efficiency"
        ].values[..., None]

        self.energy = self.energy.fillna(0)
        self.energy *= np.isfinite(self.energy)

        self["auxiliary energy"] = (
            self.energy.sel(parameter="auxiliary energy").sum(dim="second").T / distance
        ).T

        self["TtW energy"] = (
            self.energy.sel(
                parameter=["motive energy", "auxiliary energy", "recuperated energy"]
            )
            .sum(dim=["second", "parameter"])
            .T
            / distance
        ).T

        # override of TtW energy, provided by the user
        if self.energy_consumption:
            for key, val in self.energy_consumption.items():
                pwt, size, year = key
                if val is not None:

                    print(
                        f"Overriding TtW energy for {pwt} {size} {year} with {val} kj/km"
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
                        * distance.mean()  # km
                        / self.energy.shape[-1]  # seconds
                    ).reshape(
                        1, -1
                    )

        self["TtW energy, combustion mode"] = self["TtW energy"] * (
            self["combustion power share"] > 0
        )
        self["TtW energy, electric mode"] = self["TtW energy"] * (
            self["combustion power share"] == 0
        )

    def set_fuel_cell_parameters(self) -> None:
        """
        Specific setup for fuel cells, which are mild hybrids.
        Must be called after :meth:`.set_power_parameters`.
        """
        if "FCEV" in self.array.coords["powertrain"]:

            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell system efficiency")
            ] = (
                self.array.loc[
                    dict(powertrain="FCEV", parameter="fuel cell stack efficiency")
                ]
                / self.array.loc[
                    dict(powertrain="FCEV", parameter="fuel cell own consumption")
                ]
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell power share")
            ] = self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell power share")
            ].clip(
                min=0, max=1
            )
            self.array.loc[dict(powertrain="FCEV", parameter="fuel cell power")] = (
                self.array.loc[dict(powertrain="FCEV", parameter="power")]
                * self.array.loc[
                    dict(powertrain="FCEV", parameter="fuel cell power share")
                ]
                * self.array.loc[
                    dict(powertrain="FCEV", parameter="fuel cell own consumption")
                ]
            )
            # our basic fuel cell mass is based on a car fuel cell with 800 mW/cm2 and 0.51 kg/kW
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell stack mass")
            ] = (
                0.51
                * self.array.loc[dict(powertrain="FCEV", parameter="fuel cell power")]
                * (
                    800
                    / self.array.loc[
                        dict(
                            powertrain="FCEV", parameter="fuel cell power area density"
                        )
                    ]
                )
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell ancillary BoP mass")
            ] = (
                self.array.loc[dict(powertrain="FCEV", parameter="fuel cell power")]
                * self.array.loc[
                    dict(
                        powertrain="FCEV",
                        parameter="fuel cell ancillary BoP mass per power",
                    )
                ]
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell essential BoP mass")
            ] = (
                self.array.loc[dict(powertrain="FCEV", parameter="fuel cell power")]
                * self.array.loc[
                    dict(
                        powertrain="FCEV",
                        parameter="fuel cell essential BoP mass per power",
                    )
                ]
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="battery power")
            ] = self.array.loc[dict(powertrain="FCEV", parameter="fuel cell power")] * (
                1
                - self.array.loc[
                    dict(powertrain="FCEV", parameter="fuel cell power share")
                ]
            )
            self.array.loc[dict(powertrain="FCEV", parameter="battery cell mass")] = (
                self.array.loc[dict(powertrain="FCEV", parameter="battery power")]
                / self.array.loc[
                    dict(powertrain="FCEV", parameter="battery cell power density")
                ]
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="battery BoP mass")
            ] = self.array.loc[
                dict(powertrain="FCEV", parameter="battery cell mass")
            ] * (
                1
                - self.array.loc[
                    dict(powertrain="FCEV", parameter="battery cell mass share")
                ]
            )
            self.array.loc[
                dict(powertrain="FCEV", parameter="oxidation energy stored")
            ] = (
                self.array.loc[dict(powertrain="FCEV", parameter="fuel mass")]
                * 120
                / 3.6
            )  # kWh
            self.array.loc[dict(powertrain="FCEV", parameter="fuel tank mass")] = (
                self.array.loc[
                    dict(powertrain="FCEV", parameter="oxidation energy stored")
                ]
                * self.array.loc[
                    dict(powertrain="FCEV", parameter="H2 tank mass per energy")
                ]
            )

    def set_auxiliaries(self) -> None:
        """
        Calculates the power needed to operate the auxiliary services of the vehicle (heating, cooling).

        The demand for heat and cold are expressed as a fraction of the heating and cooling capacities

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
        self["recuperation efficiency"] = (
            self["drivetrain efficiency"] * self["battery charge efficiency"]
        )

    def set_battery_fuel_cell_replacements(self) -> None:
        """
        This methods calculates the fraction of the replacement battery needed to match the vehicle lifetime.

        .. note::
            if ``car lifetime`` = 200000 (km) and ``battery lifetime`` = 190000 (km) then ``replacement battery`` = 0.05

        .. note::
            It is debatable whether this is realistic or not. Car owners may not decide to invest in a new
            battery if the remaining lifetime of the car is only 10000 km. Also, a battery lifetime may be expressed
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

        # The number of fuel cell replacements is based on the average distance driven
        # with a set of fuel cells given their lifetime expressed in hours of use.
        # The number is replacement is rounded *up* as we assume no allocation of burden
        # with a second life

        if "FCEV" in self.array.coords["powertrain"]:

            self.array.loc[
                dict(powertrain="FCEV", parameter="fuel cell lifetime replacements")
            ] = np.ceil(
                np.clip(
                    (
                        self.array.loc[
                            dict(powertrain="FCEV", parameter="lifetime kilometers")
                        ]
                        / (
                            self.ecm.cycle.sum(axis=0)
                            / self.ecm.cycle.shape[0]
                            * self.array.loc[
                                dict(
                                    powertrain="FCEV",
                                    parameter="fuel cell lifetime hours",
                                )
                            ].T
                        )
                    )
                    - 1,
                    0,
                    5,
                )
            )

    def set_car_masses(self) -> None:
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
            for key, target_mass in self.target_mass.items():
                pwt, size, year = key

                if target_mass:
                    current_curb_mass = self.array.loc[
                        dict(
                            powertrain=pwt, size=size, year=year, parameter="curb mass"
                        )
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
                    ] += mass_difference / (1 - lightweighting)

            self["curb mass"] = self["glider base mass"] * (1 - self["lightweighting"])
            self["curb mass"] += self[curb_mass_includes].sum(axis=2)

        self["total cargo mass"] = (
            self["average passengers"] * self["average passenger mass"]
            + self["cargo mass"]
        )
        self["driving mass"] = self["curb mass"] + self["total cargo mass"]

    def set_power_parameters(self) -> None:
        """Set electric and combustion motor powers based on input parameter ``power to mass ratio``."""
        # Convert from W/kg to kW
        self["power"] = self["power to mass ratio"] * self["curb mass"] / 1000

        if self.power:
            for key, power in self.power.items():
                pwt, size, year = key
                if power:
                    self.array.loc[
                        dict(powertrain=pwt, size=size, year=year, parameter="power")
                    ] = power

        self["combustion power share"] = self["combustion power share"].clip(
            min=0, max=1
        )
        self["combustion power"] = self["power"] * self["combustion power share"]
        self["electric power"] = self["power"] * (1 - self["combustion power share"])

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
        """Calculate the share of recuperated energy, over the total negative motive energy"""

        self["share recuperated energy"] = (
            self.energy.loc[dict(parameter="recuperated energy")].sum(dim="second") * -1
        ) / self.energy.loc[dict(parameter="negative motive energy")].sum(dim="second")

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

    def create_PHEV(self):
        """PHEV-p/d is the range-weighted average between PHEV-c-p/PHEV-c-d and PHEV-e."""

        if "PHEV-d" in self.array.coords["powertrain"].values:

            self.array.loc[:, "PHEV-d", :, :, :] = (
                self.array.loc[:, "PHEV-e", :, :, :]
                * self.array.loc[:, "PHEV-e", "electric utility factor", :, :]
            ) + (
                self.array.loc[:, "PHEV-c-d", :, :, :]
                * (1 - self.array.loc[:, "PHEV-e", "electric utility factor", :, :])
            )

            self.array.loc[
                dict(parameter="electric utility factor", powertrain=["PHEV-d"])
            ] = self.array.loc[
                dict(parameter="electric utility factor", powertrain=["PHEV-e"])
            ].values

            self.energy.loc[
                dict(
                    parameter=[
                        "motive energy",
                        "auxiliary energy",
                        "recuperated energy",
                    ],
                    powertrain=["PHEV-d"],
                )
            ] = (
                self.array.loc[
                    dict(parameter="electric utility factor", powertrain=["PHEV-e"])
                ]
                * self.energy.loc[
                    dict(
                        parameter=[
                            "motive energy",
                            "auxiliary energy",
                            "recuperated energy",
                        ],
                        powertrain=["PHEV-e"],
                    )
                ]
            ).values.transpose(
                0, 1, 4, 2, 3, 5
            ) + (
                (
                    1
                    - self.array.loc[
                        dict(parameter="electric utility factor", powertrain="PHEV-e")
                    ]
                )
                * self.energy.loc[
                    dict(
                        parameter=[
                            "motive energy",
                            "auxiliary energy",
                            "recuperated energy",
                        ],
                        powertrain=["PHEV-c-d"],
                    )
                ]
            ).values.transpose(
                0, 3, 4, 1, 2, 5
            )

            # We need to preserve the fuel mass and fuel tank mass
            self.array.loc[
                dict(
                    parameter=[
                        "fuel mass",
                        "fuel tank mass",
                        "oxidation energy stored",
                    ],
                    powertrain="PHEV-d",
                )
            ] = self.array.loc[
                dict(
                    parameter=[
                        "fuel mass",
                        "fuel tank mass",
                        "oxidation energy stored",
                    ],
                    powertrain="PHEV-c-d",
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
                    ],
                    powertrain="PHEV-d",
                )
            ] = self.array.loc[
                dict(
                    parameter=[
                        "energy battery mass",
                        "battery BoP mass",
                        "battery cell mass",
                        "battery DoD",
                        "battery cell energy density",
                    ],
                    powertrain="PHEV-e",
                )
            ]

            # We store the tank-to-wheel energy consumption
            # in combustion and electric mode separately
            self.array.loc[
                dict(parameter="TtW energy, combustion mode", powertrain="PHEV-d")
            ] = self.array.loc[dict(parameter="TtW energy", powertrain="PHEV-c-d")]

            self.array.loc[
                dict(parameter="TtW energy, electric mode", powertrain="PHEV-d")
            ] = self.array.loc[dict(parameter="TtW energy", powertrain="PHEV-e")]

            # We need to recalculate the range as well
            self.array.loc[dict(parameter="range", powertrain="PHEV-d")] = (
                self.array.loc[
                    dict(parameter="oxidation energy stored", powertrain="PHEV-d")
                ]
                * 3600
                / self.array.loc[
                    dict(parameter="TtW energy, combustion mode", powertrain="PHEV-d")
                ]
            )

            self.array.loc[dict(parameter="range", powertrain="PHEV-d")] += (
                self.array.loc[
                    dict(parameter="electric energy stored", powertrain="PHEV-d")
                ]
                * 3600
                / self.array.loc[
                    dict(parameter="TtW energy, electric mode", powertrain="PHEV-e")
                ]
            )

        if "PHEV-p" in self.array.coords["powertrain"].values:

            self.array.loc[:, "PHEV-p", :, :, :] = (
                self.array.loc[:, "PHEV-e", :, :, :]
                * self.array.loc[:, "PHEV-e", "electric utility factor", :, :]
            ) + (
                self.array.loc[:, "PHEV-c-p", :, :, :]
                * (1 - self.array.loc[:, "PHEV-e", "electric utility factor", :, :])
            )

            self.array.loc[
                dict(parameter="electric utility factor", powertrain=["PHEV-p"])
            ] = self.array.loc[
                dict(parameter="electric utility factor", powertrain=["PHEV-e"])
            ].values

            self.energy.loc[
                dict(
                    parameter=[
                        "motive energy",
                        "auxiliary energy",
                        "recuperated energy",
                    ],
                    powertrain=["PHEV-p"],
                )
            ] = (
                self.array.loc[
                    dict(parameter="electric utility factor", powertrain=["PHEV-e"])
                ]
                * self.energy.loc[
                    dict(
                        parameter=[
                            "motive energy",
                            "auxiliary energy",
                            "recuperated energy",
                        ],
                        powertrain=["PHEV-e"],
                    )
                ]
            ).values.transpose(
                0, 1, 4, 2, 3, 5
            ) + (
                (
                    1
                    - self.array.loc[
                        dict(parameter="electric utility factor", powertrain="PHEV-e")
                    ]
                )
                * self.energy.loc[
                    dict(
                        parameter=[
                            "motive energy",
                            "auxiliary energy",
                            "recuperated energy",
                        ],
                        powertrain=["PHEV-c-p"],
                    )
                ]
            ).values.transpose(
                0, 3, 4, 1, 2, 5
            )

            # We need to preserve the fuel mass and fuel tank mass
            self.array.loc[
                dict(
                    parameter=[
                        "fuel mass",
                        "fuel tank mass",
                        "oxidation energy stored",
                    ],
                    powertrain="PHEV-p",
                )
            ] = self.array.loc[
                dict(
                    parameter=[
                        "fuel mass",
                        "fuel tank mass",
                        "oxidation energy stored",
                    ],
                    powertrain="PHEV-c-p",
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
                    ],
                    powertrain="PHEV-p",
                )
            ] = self.array.loc[
                dict(
                    parameter=[
                        "energy battery mass",
                        "battery BoP mass",
                        "battery cell mass",
                        "battery DoD",
                        "battery cell energy density",
                    ],
                    powertrain="PHEV-e",
                )
            ]

            # We store the tank-to-wheel energy consumption
            # in combustion and electric mode separately
            self.array.loc[
                dict(parameter="TtW energy, combustion mode", powertrain="PHEV-p")
            ] = self.array.loc[dict(parameter="TtW energy", powertrain="PHEV-c-p")]

            self.array.loc[
                dict(parameter="TtW energy, electric mode", powertrain="PHEV-p")
            ] = self.array.loc[dict(parameter="TtW energy", powertrain="PHEV-e")]

            # We need to recalculate the range as well
            self.array.loc[dict(parameter="range", powertrain="PHEV-p")] = (
                self.array.loc[
                    dict(parameter="oxidation energy stored", powertrain="PHEV-p")
                ]
                * 3600
                / self.array.loc[
                    dict(parameter="TtW energy, combustion mode", powertrain="PHEV-p")
                ]
            )

            self.array.loc[dict(parameter="range", powertrain="PHEV-p")] += (
                self.array.loc[
                    dict(parameter="electric energy stored", powertrain="PHEV-p")
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

        self.array.loc[dict(parameter="battery cell mass")] = (
            self.array.loc[dict(parameter="energy battery mass")]
            * self.array.loc[dict(parameter="battery cell mass share")]
        )

        self.array.loc[dict(parameter="battery BoP mass")] = self.array.loc[
            dict(parameter="energy battery mass")
        ] * (1 - self.array.loc[dict(parameter="battery cell mass share")])

    def set_storage_size_or_range(self):
        """
        Set storage size or range for each powertrain.
        :return:
        """

        if "capacity" in self.energy_storage:
            for key, val in self.energy_storage["capacity"].items():
                pwt, size, year = key
                if val:
                    self.array.loc[
                        dict(
                            parameter="battery cell mass",
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
                    )

                    self.array.loc[
                        dict(
                            parameter="energy battery mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        self.array.loc[
                            dict(
                                parameter="battery cell mass",
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

                    self.array.loc[
                        dict(
                            parameter="battery BoP mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        self.array.loc[
                            dict(
                                parameter="energy battery mass",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                        - self.array.loc[
                            dict(
                                parameter="battery cell mass",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                    )

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
                            parameter="battery cell mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        energy_stored
                        / self.array.loc[
                            dict(
                                parameter="battery cell energy density",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                    )

                    self.array.loc[
                        dict(
                            parameter="energy battery mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        self.array.loc[
                            dict(
                                parameter="battery cell mass",
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

                    self.array.loc[
                        dict(
                            parameter="battery BoP mass",
                            powertrain=pwt,
                            size=size,
                            year=year,
                        )
                    ] = (
                        self.array.loc[
                            dict(
                                parameter="energy battery mass",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                        - self.array.loc[
                            dict(
                                parameter="battery cell mass",
                                powertrain=pwt,
                                size=size,
                                year=year,
                            )
                        ]
                    )

    def set_range(self) -> None:
        """
        Calculate range autonomy of vehicles
        :return:
        """

        list_combustion = [
            i
            for i in [
                "ICEV-p",
                "HEV-p",
                "HEV-d",
                "PHEV-c-p",
                "PHEV-c-d",
                "ICEV-d",
                "ICEV-g",
                "FCEV",
            ]
            if i in self.array.powertrain
        ]

        list_electric = [i for i in ["BEV", "PHEV-e"] if i in self.array.powertrain]

        if len(list_combustion) > 0:

            fuel_mass = self.array.loc[
                dict(powertrain=list_combustion, parameter="fuel mass")
            ]  # in kg
            lhv = self.array.loc[
                dict(powertrain=list_combustion, parameter="LHV fuel MJ per kg")
            ]  # MJ/kg
            TtW = self.array.loc[
                dict(powertrain=list_combustion, parameter="TtW energy")
            ]  # kj/km
            self.array.loc[dict(powertrain=list_combustion, parameter="range")] = (
                fuel_mass * lhv * 1000
            ) / TtW  # -> km

        if len(list_electric) > 0:

            energy_stored = self.array.loc[
                dict(powertrain=list_electric, parameter="electric energy stored")
            ]  # in kWh
            battery_DoD = self.array.loc[
                dict(powertrain=list_electric, parameter="battery DoD")
            ]  # maximum depth of discharge allowed (80%)
            TtW = self.array.loc[
                dict(powertrain=list_electric, parameter="TtW energy")
            ]  # kj/km
            self.array.loc[dict(powertrain=list_electric, parameter="range")] = (
                energy_stored * battery_DoD * 3.6 * 1000 / TtW
            )  # -> km

    def set_energy_stored_properties(self) -> None:
        """
        Calculate size and capacity of onboard
        energy storage components.
        :return:
        """

        list_combustion = [
            i
            for i in [
                "ICEV-p",
                "HEV-p",
                "HEV-d",
                "PHEV-c-p",
                "PHEV-c-d",
                "ICEV-d",
            ]
            if i in self.array.powertrain
        ]

        if len(list_combustion) > 0:

            self.array.loc[
                dict(powertrain=list_combustion, parameter="oxidation energy stored")
            ] = (
                self.array.loc[dict(powertrain=list_combustion, parameter="fuel mass")]
                * self.array.loc[
                    dict(powertrain=list_combustion, parameter="LHV fuel MJ per kg")
                ]
                / 3.6
            )

            self.array.loc[
                dict(powertrain=list_combustion, parameter="fuel tank mass")
            ] = (
                self.array.loc[
                    dict(
                        powertrain=list_combustion, parameter="oxidation energy stored"
                    )
                ]
                * self.array.loc[
                    dict(
                        powertrain=list_combustion,
                        parameter="fuel tank mass per energy",
                    )
                ]
            )

        if "ICEV-g" in self.array.coords["powertrain"].values:

            self.array.loc[
                dict(powertrain="ICEV-g", parameter="oxidation energy stored")
            ] = (
                self.array.loc[dict(powertrain="ICEV-g", parameter="fuel mass")]
                * self.array.loc[
                    dict(powertrain="ICEV-g", parameter="LHV fuel MJ per kg")
                ]
                / 3.6
            )

            self.array.loc[dict(powertrain="ICEV-g", parameter="fuel tank mass")] = (
                self.array.loc[
                    dict(powertrain="ICEV-g", parameter="oxidation energy stored")
                ]
                * self.array.loc[
                    dict(powertrain="ICEV-g", parameter="CNG tank mass slope")
                ]
                + self.array.loc[
                    dict(powertrain="ICEV-g", parameter="CNG tank mass intercept")
                ]
            )

        list_electric = [i for i in ["BEV", "PHEV-e"] if i in self.array.powertrain]

        for pwt in list_electric:
            self.array.loc[dict(powertrain=pwt, parameter="electric energy stored")] = (
                self.array.loc[dict(powertrain=pwt, parameter="battery cell mass")]
                * self.array.loc[
                    dict(powertrain=pwt, parameter="battery cell energy density")
                ]
            )

        if "capacity" in self.energy_storage:
            for key, val in self.energy_storage["capacity"].items():
                pwt, size, year = key
                if val:
                    self.array.loc[
                        dict(
                            powertrain=pwt,
                            size=size,
                            year=year,
                            parameter="electric energy stored",
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

        if "BEV" in self.array.coords["powertrain"].values:
            self.array.loc[
                dict(powertrain="BEV", parameter="energy cost")
            ] /= self.array.loc[
                dict(powertrain="BEV", parameter="battery charge efficiency")
            ]

        self["component replacement cost"] = (
            self["energy battery cost"] * self["battery lifetime replacements"]
            + self["fuel cell cost"] * self["fuel cell lifetime replacements"]
        )

        with open(DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
            to_markup = yaml.safe_load(stream)["markup"]

        self[to_markup] *= self["markup factor"]

        # calculate costs per km:
        self["lifetime"] = self["lifetime kilometers"] / self["kilometers per year"]
        i = self["interest rate"]
        lifetime = self["lifetime"]
        amortisation_factor = ne.evaluate("i + (i / ((1 + i) ** lifetime - 1))")

        with open(DATA_DIR / "purchase_cost_params.yaml", "r") as stream:
            purchase_cost_params = yaml.safe_load(stream)["purchase"]

        self["purchase cost"] = self[purchase_cost_params].sum(axis=2)

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
        km_per_year = self["kilometers per year"]
        com_repl_cost = self["component replacement cost"]
        self["amortised component replacement cost"] = ne.evaluate(
            "(com_repl_cost * ((1 - i) ** lifetime / 2) * amortisation_factor / km_per_year)"
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

        self["TtW efficiency"] = (
            _(self["battery discharge efficiency"])
            * _(self["fuel cell system efficiency"])
            * self["drivetrain efficiency"]
            * self["engine efficiency"]
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
        hem = HotEmissionsModel(self.ecm.cycle, self.ecm.cycle_name)

        with open(DATA_DIR / "exhaust_flows.yaml", "r") as stream:
            list_direct_emissions = yaml.safe_load(stream)

        l_y = []
        for year in self.array.year.values:
            # European emission standards function of registration year
            # Here the registration year is synonymous of production year.
            if year < 1993:
                l_y.append(1)
            if 1993 <= year < 1997:
                l_y.append(1)
            if 1997 <= year < 2001:
                l_y.append(2)
            if 2001 <= year < 2006:
                l_y.append(3)
            if 2006 <= year < 2011:
                l_y.append(4)
            if 2011 <= year < 2015:
                l_y.append(5)
            if 2015 <= year < 2017:
                l_y.append(6.0)
            if 2017 <= year < 2019:
                l_y.append(6.1)
            if 2019 <= year < 2021:
                l_y.append(6.2)
            if year >= 2021:
                l_y.append(6.3)

        # to calculate emissions and degradation factors
        # we need the vehicle's lifetime, annual mileage
        # as well as its instant fuel consumption

        list_diesel = [
            i
            for i in ["ICEV-d", "PHEV-c-d", "HEV-d", "PHEV-d"]
            if i in self.array.powertrain
        ]

        if len(list_diesel) > 0:

            energy_consumption = self.energy.sel(
                powertrain=list_diesel,
                parameter=["motive energy", "auxiliary energy", "recuperated energy"],
            ).sum(dim="parameter")

            self.array.loc[
                dict(
                    powertrain=list_diesel,
                    parameter=list_direct_emissions,
                )
            ] = hem.get_hot_emissions(
                powertrain_type="ICEV-d",
                euro_class=l_y,
                lifetime_km=self.array.loc[
                    dict(powertrain=list_diesel, parameter="lifetime kilometers")
                ],
                energy_consumption=energy_consumption,
                yearly_km=self.array.loc[
                    dict(powertrain=list_diesel, parameter="kilometers per year")
                ],
            )

        list_petrol = [
            i
            for i in ["ICEV-p", "PHEV-c-p", "HEV-p", "PHEV-p"]
            if i in self.array.powertrain
        ]

        if len(list_petrol) > 0:

            energy_consumption = self.energy.sel(
                powertrain=list_petrol,
                parameter=["motive energy", "auxiliary energy", "recuperated energy"],
            ).sum(dim="parameter")

            self.array.loc[
                dict(
                    powertrain=list_petrol,
                    parameter=list_direct_emissions,
                )
            ] = hem.get_hot_emissions(
                powertrain_type="ICEV-p",
                euro_class=l_y,
                lifetime_km=self.array.loc[
                    dict(powertrain=list_petrol, parameter="lifetime kilometers")
                ],
                energy_consumption=energy_consumption,
                yearly_km=self.array.loc[
                    dict(powertrain=list_petrol, parameter="kilometers per year")
                ],
            )

        if "ICEV-g" in self.array.powertrain:

            energy_consumption = self.energy.sel(
                powertrain=["ICEV-g"],
                parameter=["motive energy", "auxiliary energy", "recuperated energy"],
            ).sum(dim="parameter")

            self.array.loc[
                dict(
                    powertrain=["ICEV-g"],
                    parameter=list_direct_emissions,
                )
            ] = hem.get_hot_emissions(
                powertrain_type=["ICEV-g"],
                euro_class=l_y,
                lifetime_km=self.array.loc[
                    dict(powertrain=["ICEV-g"], parameter="lifetime kilometers")
                ],
                energy_consumption=energy_consumption,
                yearly_km=self.array.loc[
                    dict(powertrain=["ICEV-g"], parameter="kilometers per year")
                ],
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
            cycle_name=self.ecm.cycle_name,
            cycle=self.ecm.cycle,
            mass=self["driving mass"],
        )

        res = pem.get_abrasion_emissions()
        self[list_param] = res

        # brake emissions are discounted by the use of regenerative braking
        self["brake wear emissions"] *= 1 - self["share recuperated energy"]

    def set_noise_emissions(self) -> None:
        """
        Calculate noise emissions based on ``driving cycle``.
        The driving cycle is passed to the :class:`NoiseEmissionsModel` class and :meth:`get_sound_power_per_compartment`
        returns emissions per compartment type ("rural", "non-urban" and "urban") per second of driving cycle.

        Noise emissions are not differentiated by size classes at the moment, but only by powertrain "type"
        (e.g., combustion, hybrid and electric)

        :return: Does not return anything. Modifies ``self.array`` in place.
        """
        nem = NoiseEmissionsModel(self.ecm.cycle, self.ecm.cycle_name)

        with open(DATA_DIR / "noise_flows.yaml", "r") as stream:
            list_noise_emissions = yaml.safe_load(stream)

        list_pwt = [
            i
            for i in [
                "ICEV-p",
                "PHEV-c-p",
                "ICEV-g",
                "ICEV-d",
                "PHEV-c-d",
            ]
            if i in self.array.powertrain
        ]

        self.array.loc[
            dict(powertrain=list_pwt, parameter=list_noise_emissions)
        ] = nem.get_sound_power_per_compartment("combustion")

        list_pwt = [i for i in ["BEV", "FCEV", "PHEV-e"] if i in self.array.powertrain]

        self.array.loc[
            dict(powertrain=list_pwt, parameter=list_noise_emissions)
        ] = nem.get_sound_power_per_compartment("electric")

        list_pwt = [i for i in ["HEV-p", "HEV-d"] if i in self.array.powertrain]

        self.array.loc[
            dict(powertrain=list_pwt, parameter=list_noise_emissions)
        ] = nem.get_sound_power_per_compartment("hybrid")

    def calculate_cost_impacts(self, sensitivity=False, scope=None) -> xr.DataArray:
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

        if scope is None:
            scope = {
                "size": self.array.coords["size"].values.tolist(),
                "powertrain": self.array.coords["powertrain"].values.tolist(),
                "year": self.array.coords["year"].values.tolist(),
            }

        list_cost_cat = [
            "purchase",
            "maintenance",
            "component replacement",
            "energy",
            "total",
        ]

        response = self.array.sel(
            powertrain=scope["powertrain"],
            size=scope["size"],
            year=scope["year"],
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
