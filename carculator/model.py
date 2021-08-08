import numexpr as ne
import numpy as np
import xarray as xr

from .energy_consumption import EnergyConsumptionModel
from .hot_emissions import HotEmissionsModel
from .noise_emissions import NoiseEmissionsModel

DEFAULT_MAPPINGS = {
    "electric": {"BEV", "PHEV-e"},
    "combustion": {
        "ICEV-p",
        "HEV-p",
        "HEV-d",
        "PHEV-c-p",
        "ICEV-g",
        "ICEV-d",
        "PHEV-c-d",
    },
    "combustion_wo_cng": {"ICEV-p", "HEV-p", "HEV-d", "PHEV-c-p", "ICEV-d", "PHEV-c-d"},
    "pure_combustion": {"ICEV-p", "ICEV-g", "ICEV-d"},
    "petrol": {"ICEV-p", "HEV-p", "PHEV-c-p"},
    "cng": {"ICEV-g"},
    "fuel_cell": {"FCEV"},
    "hybrid": {"PHEV-c-p", "PHEV-e", "PHEV-c-d"},
    "combustion_hybrid": {"PHEV-c-p", "PHEV-c-d"},
    "electric_hybrid": {"PHEV-e"},
    "diesel": {"ICEV-d", "PHEV-c-d", "HEV-d"},
    "battery": {"BEV"},
}


def finite(array, mask_value=0):
    return np.where(np.isfinite(array), array, mask_value)


class CarModel:

    """
    This class represents the entirety of the vehicles considered, with useful attributes, such as an array that stores
    all the vehicles parameters.

    :ivar array: multi-dimensional numpy-like array that contains parameters' value(s)
    :vartype array: xarray.DataArray
    :ivar mappings: Dictionary with names correspondence
    :vartype mappings: dict
    :ivar ecm: instance of :class:`EnergyConsumptionModel` class for a given driving cycle
    :vartype ecm: coarse.energy_consumption.EnergyConsumptionModel

    """

    def __init__(
        self, array, mappings=None, cycle=None, gradient=None, energy_storage=None
    ):

        self.array = array
        self.mappings = mappings or DEFAULT_MAPPINGS

        if cycle is None:
            self.ecm = EnergyConsumptionModel("WLTC")
        else:
            self.ecm = EnergyConsumptionModel(cycle=cycle, gradient=gradient)

        self.energy_storage = energy_storage or {
            "electric": {"BEV": "NMC-111", "PHEV-e": "NMC-111"}
        }

    def __call__(self, key):
        """
        This method fixes a dimension of the `array` attribute given a powertrain technology selected.

        Set up this class as a context manager, so we can have some nice syntax

        .. code-block:: python

            with class('some powertrain') as cpm:
                cpm['something']. # Will be filtered for the correct powertrain

        On with block exit, this filter is cleared
        https://stackoverflow.com/a/10252925/164864

        :param key: A powertrain type, e.g., "FCEV"
        :type key: str
        :return: An instance of `array` filtered after the powertrain selected.

        """
        self.__cache = self.array
        self.array = self.array.sel(powertrain=key)
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.array = self.__cache
        del self.__cache

    def __getitem__(self, key):
        """
        Make class['foo'] automatically filter for the parameter 'foo'
        Makes the model code much cleaner

        :param key: Parameter name
        :type key: str
        :return: `array` filtered after the parameter selected
        """

        return self.array.sel(parameter=key)

    def __setitem__(self, key, value):
        self.array.loc[{"parameter": key}] = value

    # Make it easier/more flexible to filter by powertrain types
    def __getattr__(self, key):
        if key in self.mappings:
            return self.mappings[key]

    def set_all(self, drop_hybrids=True, electric_utility_factor=None):
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
        # TODO: Converging towards a satisfying curb mass is taking too long! Needs to be optimized.

        diff = 1.0

        while diff > 0.01:
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
        self.adjust_cost()
        self.set_range()
        self.set_electric_utility_factor(electric_utility_factor)
        self.set_electricity_consumption()
        self.set_costs()
        self.set_hot_emissions()
        self.set_noise_emissions()
        self.create_PHEV()
        if drop_hybrids:
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

    def adjust_cost(self):
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
        self.array.loc[
            :,
            [pt for pt in ["FCEV"] if pt in self.array.coords["powertrain"].values],
            "fuel tank cost per kg",
            :,
            :,
        ] = np.reshape(
            (1.078e58 * np.exp(-6.32e-2 * self.array.year.values) + 3.43e2)
            * cost_factor,
            (1, 1, n_year, n_iterations),
        )

        # Correction of fuel cell stack cost, per kW
        self.array.loc[
            :,
            [pt for pt in ["FCEV"] if pt in self.array.coords["powertrain"].values],
            "fuel cell cost per kW",
            :,
            :,
        ] = np.reshape(
            (3.15e66 * np.exp(-7.35e-2 * self.array.year.values) + 2.39e1)
            * cost_factor,
            (1, 1, n_year, n_iterations),
        )

        # Correction of energy battery system cost, per kWh
        self.array.loc[
            :,
            [
                pt
                for pt in ["BEV", "PHEV-e", "PHEV-c-p", "PHEV-c-d"]
                if pt in self.array.coords["powertrain"].values
            ],
            "energy battery cost per kWh",
            :,
            :,
        ] = np.reshape(
            (2.75e86 * np.exp(-9.61e-2 * self.array.year.values) + 5.059e1)
            * cost_factor,
            (1, 1, n_year, n_iterations),
        )

        # Correction of power battery system cost, per kW
        self.array.loc[
            :,
            [
                pt
                for pt in [
                    "ICEV-p",
                    "ICEV-d",
                    "ICEV-g",
                    "PHEV-c-p",
                    "PHEV-c-d",
                    "FCEV",
                    "HEV-p",
                    "HEV-d",
                ]
                if pt in self.array.coords["powertrain"].values
            ],
            "power battery cost per kW",
            :,
            :,
        ] = np.reshape(
            (8.337e40 * np.exp(-4.49e-2 * self.array.year.values) + 11.17)
            * cost_factor,
            (1, 1, n_year, n_iterations),
        )

        # Correction of combustion powertrain cost for ICEV-g
        self.array.loc[
            :,
            [pt for pt in ["ICEV-g"] if pt in self.array.coords["powertrain"].values],
            "combustion powertrain cost per kW",
            :,
            :,
        ] = np.clip(
            np.reshape(
                (5.92e160 * np.exp(-0.1819 * self.array.year.values) + 26.76)
                * cost_factor,
                (1, 1, n_year, n_iterations),
            ),
            None,
            100,
        )

    def adjust_fuel_mass(self):
        """
        This method adjusts the fuel mass over the years, to correct for the linear
        interpolation between years.

        """

        n_iterations = self.array.shape[-1]
        n_year = len(self.array.year.values)

        # Correction of hydrogen mass
        self.array.loc[
            :,
            [pt for pt in ["FCEV"] if pt in self.array.coords["powertrain"].values],
            "fuel mass",
            :,
            :,
        ] = np.reshape(
            (1.078e58 * np.exp(-6.32e-2 * self.array.year.values) + 3.43e2),
            (1, 1, n_year, n_iterations),
        )

        # Correction of CNG mass

    def drop_hybrid(self):
        """
        This method drops the powertrains `PHEV-c-p`, `PHEV-c-d` and `PHEV-e` as they were only used to create the
        `PHEV` powertrain.
        :returns: Does not return anything. Modifies ``self.array`` in place.
        """
        self.array = self.array.sel(
            powertrain=[
                pt
                for pt in [
                    "ICEV-p",
                    "ICEV-d",
                    "ICEV-g",
                    "PHEV-p",
                    "PHEV-d",
                    "FCEV",
                    "BEV",
                    "HEV-p",
                    "HEV-d",
                ]
                if pt in self.array.coords["powertrain"].values
            ]
        )

    def set_electricity_consumption(self):
        """
        This method calculates the total electricity consumption for BEV and plugin-hybrid vehicles
        :returns: Does not return anything. Modifies ``self.array`` in place.
        """

        for pt in self.electric:
            if pt in self.array.coords["powertrain"].values:
                with self(pt) as cpm:
                    cpm["electricity consumption"] = (
                        cpm["TtW energy"] / cpm["battery charge efficiency"]
                    ) / 3600

    def calculate_ttw_energy(self):
        """
        This method calculates the energy required to operate auxiliary services as well
        as to move the car. The sum is stored under the parameter label "TtW energy" in :attr:`self.array`.

        """

        self.energy = xr.DataArray(
            np.zeros(
                (
                    len(self.array.coords["size"]),
                    len(self.array.coords["powertrain"]),
                    3,
                    len(self.array.coords["year"]),
                    len(self.array.coords["value"]),
                    self.ecm.cycle.shape[0],
                )
            ).astype("float32"),
            coords=[
                self.array.coords["size"],
                self.array.coords["powertrain"],
                ["auxiliary energy", "motive energy", "recuperated energy"],
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
            ttw_efficiency=self["TtW efficiency"],
            recuperation_efficiency=self["recuperation efficiency"],
            motor_power=self["electric power"],
        )

        self.energy.loc[dict(parameter="motive energy")] = np.clip(
            motive_energy / 1000, 0, None
        )

        self.energy.loc[dict(parameter="motive energy")] /= self["TtW efficiency"]

        self.energy.loc[dict(parameter="recuperated energy")] = np.clip(
            recuperated_energy / 1000, self["power"].values[..., None] * -1, 0
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

        self["TtW energy, combustion mode"] = self["TtW energy"] * (
            self["combustion power share"] > 0
        )
        self["TtW energy, electric mode"] = self["TtW energy"] * (
            self["combustion power share"] == 0
        )

    def set_fuel_cell_parameters(self):
        """
        Specific setup for fuel cells, which are mild hybrids.
        Must be called after :meth:`.set_power_parameters`.
        """
        for pt in self.fuel_cell:
            if pt in self.array.coords["powertrain"].values:
                with self(pt):
                    self["fuel cell system efficiency"] = (
                        self["fuel cell stack efficiency"]
                        / self["fuel cell own consumption"]
                    )
                    self["fuel cell power share"] = self["fuel cell power share"].clip(
                        min=0, max=1
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
                        * (800 / self["fuel cell power area density"])
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

                    self["oxidation energy stored"] = (
                        self["fuel mass"] * 120 / 3.6
                    )  # kWh
                    self["fuel tank mass"] = (
                        self["oxidation energy stored"]
                        * self["H2 tank mass per energy"]
                    )

    def set_auxiliaries(self):
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

    def set_battery_fuel_cell_replacements(self):
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

        if "FCEV" in self.array.coords["powertrain"].values:
            with self("FCEV") as pt:
                pt["fuel cell lifetime replacements"] = np.ceil(
                    np.clip(
                        (
                            pt["lifetime kilometers"]
                            / (
                                pt.ecm.cycle.sum(axis=0)
                                / pt.ecm.cycle.shape[0]
                                * pt["fuel cell lifetime hours"].T
                            )
                        )
                        - 1,
                        0,
                        5,
                    )
                )

    def set_car_masses(self):
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

        self["total cargo mass"] = (
            self["average passengers"] * self["average passenger mass"]
            + self["cargo mass"]
        )
        self["driving mass"] = self["curb mass"] + self["total cargo mass"]

    def set_power_parameters(self):
        """Set electric and combustion motor powers based on input parameter ``power to mass ratio``."""
        # Convert from W/kg to kW
        self["power"] = self["power to mass ratio"] * self["curb mass"] / 1000
        self["combustion power share"] = self["combustion power share"].clip(
            min=0, max=1
        )
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
        ) * (self["electric power"] > 0)
        self["powertrain mass"] = (
            self["power"] * self["powertrain mass per power"]
            + self["powertrain fixed mass"]
        )

    def set_electric_utility_factor(self, uf=None):
        """Set the electric utility factor according to a sampled values in Germany (ICTT 2020)
        https://theicct.org/sites/default/files/publications/PHEV-white%20paper-sept2020-0.pdf

        Real-world range in simulation 20 km 30 km 40 km 50 km 60 km 70 km 80 km
        Observed UF for Germany (Sample-size weighted regression ± 2 standard errors)
        Observed UF private (in %) 30±2 41±2 50±3 58±3 65±3 71±3 75±3

        which correlated the share of km driven in electric-mode to the capacity of the battery
        (the range that can be driven in battery-depleting mode).

        The argument `uf` is used to override this relation, if needed.
        `uf` must be a ratio between 0 and .75, for each."""
        if "PHEV-e" in self.array.coords["powertrain"].values:
            with self("PHEV-e") as cpm:
                if uf is None:
                    cpm["electric utility factor"] = np.clip(
                        np.interp(
                            cpm["range"],
                            [0, 20, 30, 40, 50, 60, 70, 80],
                            [0, 0.3, 0.41, 0.5, 0.58, 0.65, 0.71, 0.75],
                        ),
                        0,
                        0.75,
                    )
                else:
                    cpm["electric utility factor"] = np.array(uf).reshape([1, -1, 1])

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

    def set_battery_properties(self):
        pt_list = [
            pt
            for pt in ["ICEV-p", "HEV-p", "HEV-d", "ICEV-g", "ICEV-d"]
            if pt in self.array.coords["powertrain"].values
        ]
        self.array.loc[:, pt_list, "battery power"] = self.array.loc[
            :, pt_list, "electric power"
        ]

        self.array.loc[:, pt_list, "battery cell mass"] = (
            self.array.loc[:, pt_list, "battery power"]
            / self.array.loc[:, pt_list, "battery cell power density"]
        )

        self["battery cell mass share"] = self["battery cell mass share"].clip(
            min=0, max=1
        )
        self.array.loc[:, pt_list, "battery BoP mass", :, :] = (
            self.array.loc[
                :,
                pt_list,
                "battery cell mass",
            ]
            * (1 - self.array.loc[:, pt_list, "battery cell mass share", :, :])
        )

        list_pt_el = [
            pt
            for pt in ["BEV", "PHEV-c-p", "PHEV-c-d", "PHEV-e"]
            if pt in self.array.coords["powertrain"].values
        ]

        self.array.loc[:, list_pt_el, "battery cell mass"] = (
            self.array.loc[:, list_pt_el, "energy battery mass"]
            * self.array.loc[:, list_pt_el, "battery cell mass share"]
        )

        self.array.loc[:, list_pt_el, "battery BoP mass"] = self.array.loc[
            :, list_pt_el, "energy battery mass"
        ] * (1 - self.array.loc[:, list_pt_el, "battery cell mass share"])

    def set_range(self):

        list_pt = [
            pt
            for pt in [
                "ICEV-p",
                "HEV-p",
                "HEV-d",
                "PHEV-c-p",
                "PHEV-c-d",
                "ICEV-d",
                "ICEV-g",
                "FCEV",
            ]
            if pt in self.array.coords["powertrain"].values
        ]

        list_pt_el = [
            pt
            for pt in ["BEV", "PHEV-e"]
            if pt in self.array.coords["powertrain"].values
        ]

        fuel_mass = self.array.loc[:, list_pt, "fuel mass"]
        lhv = self.array.loc[:, list_pt, "LHV fuel MJ per kg"]

        energy_stored = self.array.loc[:, list_pt_el, "electric energy stored"]
        battery_DoD = self.array.loc[:, list_pt_el, "battery DoD"]

        TtW_el = self.array.loc[:, list_pt_el, "TtW energy"]
        TtW = self.array.loc[:, list_pt, "TtW energy"]

        self.array.loc[:, list_pt, "range"] = ne.evaluate(
            "(fuel_mass * lhv * 1000) / TtW"
        )
        self.array.loc[:, list_pt_el, "range"] = ne.evaluate(
            "(energy_stored * battery_DoD * 3.6 * 1000) / TtW_el"
        )

    def set_energy_stored_properties(self):

        list_combustion = [
            pt
            for pt in ["ICEV-p", "HEV-p", "HEV-d", "PHEV-c-p", "PHEV-c-d", "ICEV-d"]
            if pt in self.array.coords["powertrain"].values
        ]

        self.array.loc[:, list_combustion, "oxidation energy stored"] = (
            self.array.loc[:, list_combustion, "fuel mass"]
            * self.array.loc[:, list_combustion, "LHV fuel MJ per kg"]
            / 3.6
        )
        self.array.loc[:, list_combustion, "fuel tank mass"] = (
            self.array.loc[:, list_combustion, "oxidation energy stored"]
            * self.array.loc[:, list_combustion, "fuel tank mass per energy"]
        )

        if "ICEV-g" in self.array.coords["powertrain"].values:
            self.array.loc[:, "ICEV-g", "oxidation energy stored"] = (
                self.array.loc[:, "ICEV-g", "fuel mass"]
                * self.array.loc[:, "ICEV-g", "LHV fuel MJ per kg"]
                / 3.6
            )

            self.array.loc[:, "ICEV-g", "fuel tank mass"] = (
                self.array.loc[:, "ICEV-g", "oxidation energy stored"]
                * self.array.loc[:, "ICEV-g", "CNG tank mass slope"]
                + self.array.loc[:, "ICEV-g", "CNG tank mass intercept"]
            )

        for pt in self.battery:
            if pt in self.array.coords["powertrain"].values:
                with self(pt) as cpm:
                    battery_tech_label = (
                        "battery cell energy density, "
                        + self.energy_storage["electric"][pt]
                    )
                    cpm["electric energy stored"] = (
                        cpm["battery cell mass"] * cpm[battery_tech_label]
                    )

        for pt in self.electric_hybrid:
            if pt in self.array.coords["powertrain"].values:
                with self(pt) as cpm:
                    battery_tech_label = (
                        "battery cell energy density, "
                        + self.energy_storage["electric"][pt]
                    )
                    cpm["electric energy stored"] = (
                        cpm["battery cell mass"] * cpm[battery_tech_label]
                    )
                    cpm["fuel tank mass"] = (
                        cpm["fuel mass"]
                        * cpm["LHV fuel MJ per kg"]
                        / 3.6
                        * cpm["fuel tank mass per energy"]
                    )

        # kWh electricity/kg battery cell
        self["battery cell production energy electricity share"] = self[
            "battery cell production energy electricity share"
        ].clip(min=0, max=1)
        self["battery cell production electricity"] = (
            self["battery cell production energy"]
            * self["battery cell production energy electricity share"]
        )
        # MJ heat/kg battery cell
        self["battery cell production heat"] = (
            self["battery cell production energy"]
            - self["battery cell production electricity"]
        ) * 3.6

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
            self["energy battery cost per kWh"] * self["electric energy stored"]
        )
        self["fuel tank cost"] = self["fuel tank cost per kg"] * self["fuel mass"]
        # Per km
        self["energy cost"] = self["energy cost per kWh"] * self["TtW energy"] / 3600

        # For battery, need to divide cost of electricity in battery by efficiency of charging
        for pt in self.battery:
            if pt in self.array.coords["powertrain"].values:
                with self(pt):
                    self["energy cost"] /= self["battery charge efficiency"]

        self["component replacement cost"] = (
            self["energy battery cost"] * self["battery lifetime replacements"]
            + self["fuel cell cost"] * self["fuel cell lifetime replacements"]
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

        self[to_markup] *= self["markup factor"]

        # calculate costs per km:
        self["lifetime"] = self["lifetime kilometers"] / self["kilometers per year"]
        i = self["interest rate"]
        lifetime = self["lifetime"]
        amortisation_factor = ne.evaluate("i + (i / ((1 + i) ** lifetime - 1))")

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

        self["purchase cost"] = self[purchase_cost_list].sum(axis=2)

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

    def set_ttw_efficiency(self):
        _ = lambda array: np.where(array == 0, 1, array)
        # TODO> check if battery charge efficiency should be added
        self["TtW efficiency"] = (
            _(self["battery discharge efficiency"])
            * _(self["fuel cell system efficiency"])
            * self["drivetrain efficiency"]
            * self["engine efficiency"]
        )

    def set_hot_emissions(self):
        """
        Calculate hot pollutant emissions based on ``driving cycle``.
        The driving cycle is passed to the :class:`HotEmissionsModel` class and :meth:`get_emissions_per_powertrain`
        return emissions per substance per second of driving cycle.
        :return: Does not return anything. Modifies ``self.array`` in place.
        """
        hem = HotEmissionsModel(self.ecm.cycle, self.ecm.cycle_name)

        list_direct_emissions = [
            "Hydrocarbons direct emissions, urban",
            "Carbon monoxide direct emissions, urban",
            "Nitrogen oxides direct emissions, urban",
            "Particulate matters direct emissions, urban",
            "Methane direct emissions, urban",
            "NMVOC direct emissions, urban",
            "Dinitrogen oxide direct emissions, urban",
            "Ammonia direct emissions, urban",
            "Lead direct emissions, urban",
            "Benzene direct emissions, urban",
            "Ethane direct emissions, urban",
            "Propane direct emissions, urban",
            "Butane direct emissions, urban",
            "Pentane direct emissions, urban",
            "Hexane direct emissions, urban",
            "Cyclohexane direct emissions, urban",
            "Heptane direct emissions, urban",
            "Ethene direct emissions, urban",
            "Propene direct emissions, urban",
            "1-Pentene direct emissions, urban",
            "Toluene direct emissions, urban",
            "m-Xylene direct emissions, urban",
            "o-Xylene direct emissions, urban",
            "Formaldehyde direct emissions, urban",
            "Acetaldehyde direct emissions, urban",
            "Benzaldehyde direct emissions, urban",
            "Acetone direct emissions, urban",
            "Methyl ethyl ketone direct emissions, urban",
            "Acrolein direct emissions, urban",
            "Styrene direct emissions, urban",
            "PAH, polycyclic aromatic hydrocarbons direct emissions, urban",
            "Arsenic direct emissions, urban",
            "Selenium direct emissions, urban",
            "Zinc direct emissions, urban",
            "Copper direct emissions, urban",
            "Nickel direct emissions, urban",
            "Chromium direct emissions, urban",
            "Chromium VI direct emissions, urban",
            "Mercury direct emissions, urban",
            "Cadmium direct emissions, urban",
            "Hydrocarbons direct emissions, suburban",
            "Carbon monoxide direct emissions, suburban",
            "Nitrogen oxides direct emissions, suburban",
            "Particulate matters direct emissions, suburban",
            "Methane direct emissions, suburban",
            "NMVOC direct emissions, suburban",
            "Dinitrogen oxide direct emissions, suburban",
            "Ammonia direct emissions, suburban",
            "Lead direct emissions, suburban",
            "Benzene direct emissions, suburban",
            "Ethane direct emissions, suburban",
            "Propane direct emissions, suburban",
            "Butane direct emissions, suburban",
            "Pentane direct emissions, suburban",
            "Hexane direct emissions, suburban",
            "Cyclohexane direct emissions, suburban",
            "Heptane direct emissions, suburban",
            "Ethene direct emissions, suburban",
            "Propene direct emissions, suburban",
            "1-Pentene direct emissions, suburban",
            "Toluene direct emissions, suburban",
            "m-Xylene direct emissions, suburban",
            "o-Xylene direct emissions, suburban",
            "Formaldehyde direct emissions, suburban",
            "Acetaldehyde direct emissions, suburban",
            "Benzaldehyde direct emissions, suburban",
            "Acetone direct emissions, suburban",
            "Methyl ethyl ketone direct emissions, suburban",
            "Acrolein direct emissions, suburban",
            "Styrene direct emissions, suburban",
            "PAH, polycyclic aromatic hydrocarbons direct emissions, suburban",
            "Arsenic direct emissions, suburban",
            "Selenium direct emissions, suburban",
            "Zinc direct emissions, suburban",
            "Copper direct emissions, suburban",
            "Nickel direct emissions, suburban",
            "Chromium direct emissions, suburban",
            "Chromium VI direct emissions, suburban",
            "Mercury direct emissions, suburban",
            "Cadmium direct emissions, suburban",
            "Hydrocarbons direct emissions, rural",
            "Carbon monoxide direct emissions, rural",
            "Nitrogen oxides direct emissions, rural",
            "Particulate matters direct emissions, rural",
            "Methane direct emissions, rural",
            "NMVOC direct emissions, rural",
            "Dinitrogen oxide direct emissions, rural",
            "Ammonia direct emissions, rural",
            "Lead direct emissions, rural",
            "Benzene direct emissions, rural",
            "Ethane direct emissions, rural",
            "Propane direct emissions, rural",
            "Butane direct emissions, rural",
            "Pentane direct emissions, rural",
            "Hexane direct emissions, rural",
            "Cyclohexane direct emissions, rural",
            "Heptane direct emissions, rural",
            "Ethene direct emissions, rural",
            "Propene direct emissions, rural",
            "1-Pentene direct emissions, rural",
            "Toluene direct emissions, rural",
            "m-Xylene direct emissions, rural",
            "o-Xylene direct emissions, rural",
            "Formaldehyde direct emissions, rural",
            "Acetaldehyde direct emissions, rural",
            "Benzaldehyde direct emissions, rural",
            "Acetone direct emissions, rural",
            "Methyl ethyl ketone direct emissions, rural",
            "Acrolein direct emissions, rural",
            "Styrene direct emissions, rural",
            "PAH, polycyclic aromatic hydrocarbons direct emissions, rural",
            "Arsenic direct emissions, rural",
            "Selenium direct emissions, rural",
            "Zinc direct emissions, rural",
            "Copper direct emissions, rural",
            "Nickel direct emissions, rural",
            "Chromium direct emissions, rural",
            "Chromium VI direct emissions, rural",
            "Mercury direct emissions, rural",
            "Cadmium direct emissions, rural",
        ]

        l_y = []
        for y in self.array.year.values:
            # European emission standards funciton of registration year
            if y < 1993:
                l_y.append(1)
            if 1993 <= y < 1997:
                l_y.append(1)
            if 1997 <= y < 2001:
                l_y.append(2)
            if 2001 <= y < 2006:
                l_y.append(3)
            if 2006 <= y < 2011:
                l_y.append(4)
            if 2011 <= y < 2015:
                l_y.append(5)
            if 2015 <= y < 2017:
                l_y.append(6.0)
            if 2017 <= y < 2019:
                l_y.append(6.1)
            if 2019 <= y < 2021:
                l_y.append(6.2)
            if y >= 2021:
                l_y.append(6.3)

        l_pt = []
        for pt in self.array.powertrain.values:
            if pt in ["ICEV-d", "PHEV-c-d", "HEV-d", "PHEV-d"]:
                l_pt.append("ICEV-d")
            if pt in ["ICEV-p", "PHEV-c-p", "HEV-p", "PHEV-p"]:
                l_pt.append("ICEV-p")
            if pt == "ICEV-g":
                l_pt.append("ICEV-g")

        self.array.loc[
            dict(
                powertrain=[
                    pt
                    for pt in self.array.powertrain.values
                    if pt not in ["FCEV", "BEV", "PHEV-e"]
                ],
                parameter=list_direct_emissions,
            )
        ] = hem.get_hot_emissions(
            powertrain_type=l_pt,
            euro_class=l_y,
            lifetime_km=self.array.loc[
                dict(
                    powertrain=[
                        pt
                        for pt in self.array.powertrain.values
                        if pt not in ["FCEV", "BEV", "PHEV-e"]
                    ],
                    parameter="lifetime kilometers",
                )
            ],
            energy_consumption=self.energy.sel(
                powertrain=[
                    pt
                    for pt in self.array.powertrain.values
                    if pt not in ["FCEV", "BEV", "PHEV-e"]
                ],
                parameter=["motive energy", "auxiliary energy", "recuperated energy"],
            ).sum(dim="parameter"),
            yearly_km=self.array.sel(
                parameter="kilometers per year",
                powertrain=[
                    pt
                    for pt in self.array.powertrain.values
                    if pt not in ["FCEV", "BEV", "PHEV-e"]
                ],
            ),
        )

    def set_noise_emissions(self):
        """
        Calculate noise emissions based on ``driving cycle``.
        The driving cycle is passed to the :class:`NoiseEmissionsModel` class and :meth:`get_sound_power_per_compartment`
        returns emissions per compartment type ("rural", "non-urban" and "urban") per second of driving cycle.

        Noise emissions are not differentiated by size classes at the moment, but only by powertrain "type"
        (e.g., combustion, hybrid and electric)

        :return: Does not return anything. Modifies ``self.array`` in place.
        """
        nem = NoiseEmissionsModel(self.ecm.cycle, self.ecm.cycle_name)

        list_noise_emissions = [
            "noise, octave 1, day time, urban",
            "noise, octave 2, day time, urban",
            "noise, octave 3, day time, urban",
            "noise, octave 4, day time, urban",
            "noise, octave 5, day time, urban",
            "noise, octave 6, day time, urban",
            "noise, octave 7, day time, urban",
            "noise, octave 8, day time, urban",
            "noise, octave 1, day time, suburban",
            "noise, octave 2, day time, suburban",
            "noise, octave 3, day time, suburban",
            "noise, octave 4, day time, suburban",
            "noise, octave 5, day time, suburban",
            "noise, octave 6, day time, suburban",
            "noise, octave 7, day time, suburban",
            "noise, octave 8, day time, suburban",
            "noise, octave 1, day time, rural",
            "noise, octave 2, day time, rural",
            "noise, octave 3, day time, rural",
            "noise, octave 4, day time, rural",
            "noise, octave 5, day time, rural",
            "noise, octave 6, day time, rural",
            "noise, octave 7, day time, rural",
            "noise, octave 8, day time, rural",
        ]

        self.array.loc[
            :,
            [
                pt
                for pt in list(self.combustion)
                if pt in self.array.coords["powertrain"].values
            ],
            list_noise_emissions,
            :,
            :,
        ] = nem.get_sound_power_per_compartment("combustion").reshape((24, 1, 1))

        self.array.loc[
            :,
            [
                pt
                for pt in list(self.electric)
                if pt in self.array.coords["powertrain"].values
            ],
            list_noise_emissions,
            :,
            :,
        ] = nem.get_sound_power_per_compartment("electric").reshape((24, 1, 1))

        self.array.loc[
            :,
            [
                pt
                for pt in list(self.fuel_cell)
                if pt in self.array.coords["powertrain"].values
            ],
            list_noise_emissions,
            :,
            :,
        ] = nem.get_sound_power_per_compartment("electric").reshape((24, 1, 1))

        self.array.loc[
            :,
            [
                pt
                for pt in ["HEV-p", "HEV-d"]
                if pt in self.array.coords["powertrain"].values
            ],
            list_noise_emissions,
            :,
            :,
        ] = nem.get_sound_power_per_compartment("hybrid").reshape((24, 1, 1))

    def calculate_cost_impacts(self, sensitivity=False, scope=None):
        """
        This method returns an array with cost values per vehicle-km, sub-divided into the following groups:

            * Purchase
            * Maintentance
            * Component replacement
            * Energy
            * Total cost of ownership

        :return: A xarray array with cost information per vehicle-km
        :rtype: xarray.core.dataarray.DataArray
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

        response = xr.DataArray(
            np.zeros(
                (
                    len(scope["size"]),
                    len(scope["powertrain"]),
                    len(list_cost_cat),
                    len(scope["year"]),
                    len(self.array.coords["value"].values),
                )
            ),
            coords=[
                scope["size"],
                scope["powertrain"],
                ["purchase", "maintenance", "component replacement", "energy", "total"],
                scope["year"],
                self.array.coords["value"].values.tolist(),
            ],
            dims=["size", "powertrain", "cost_type", "year", "value"],
        )

        response.loc[
            :,
            :,
            ["purchase", "maintenance", "component replacement", "energy", "total"],
            :,
            :,
        ] = self.array.sel(
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
        ).values

        if not sensitivity:
            return response
        else:
            return response / response.sel(value="reference")
