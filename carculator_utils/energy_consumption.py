"""
energy_consumption.py contains the class EnergyConsumption Model
which exposes two methods:
one for calculating the auxiliary energy needs,
and another one for calculating the motive
energy needs.
"""

import csv
from typing import Any, List, Tuple, Union

import numexpr as ne
import numpy as np
import pandas as pd
import xarray as xr
from xarray import DataArray

from . import DATA_DIR
from .driving_cycles import (
    get_driving_cycle_specs,
    get_standard_driving_cycle_and_gradient,
)

MONTHLY_AVG_TEMP = "monthly_avg_temp.csv"


def _(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, -1)

    return obj


def __(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a heading dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, 0)

    return obj


def get_default_driving_cycle_name(vehicle_type) -> str:
    """Get the default driving cycle name"""
    return list(get_driving_cycle_specs()["columns"][vehicle_type].keys())[0]


class EnergyConsumptionModel:
    """
    Calculate energy consumption of a vehicle for a
    given driving cycle and vehicle parameters.

    Based on a selected driving cycle, this class calculates
    the acceleration needed and provides
    two methods:

        - :func:`~energy_consumption.EnergyConsumptionModel.aux_energy_per_km` calculates the energy needed to power auxiliary services
        - :func:`~energy_consumption.EnergyConsumptionModel.motive_energy_per_km` calculates the energy needed to move the vehicle over 1 km

    Acceleration is calculated as the difference between velocity
    at t_2 and velocity at t_0, divided by 2.
    See for example: http://www.unece.org/fileadmin/DAM/trans/doc/2012/wp29grpe/WLTP-DHC-12-07e.xls

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
    :type cycle: np.ndarray
    :param rho_air: Mass per unit volume of air. Set to (1.225 kg/m3) by default.
    :type rho_air: float
    :param gradient: Road gradient per second of driving, in degrees.
    None by default. Should be passed as an array of length equal
    to the length of the driving cycle.
    :type gradient: numpy.ndarray

    :ivar rho_air: Mass per unit volume of air. Value of 1.204 at 23C (test temperature for WLTC).
    :vartype rho_air: float
    :ivar velocity: Time series of speed values, in meters per second.
    :vartype velocity: numpy.ndarray
    :ivar acceleration: Time series of acceleration, calculated as
    increment in velocity per interval of 1 second, in meter per second^2.
    :vartype acceleration: numpy.ndarray

    """

    def __init__(
        self,
        vehicle_type: str,
        vehicle_size: List[str],
        cycle: Union[str, np.ndarray],
        gradient: Union[str, np.ndarray],
        rho_air: float = 1.204,
        country: str = "CH",
    ) -> None:

        self.rho_air = rho_air
        self.cycle_name = (
            cycle
            if isinstance(cycle, str)
            else get_default_driving_cycle_name(vehicle_type)
        )
        self.country = country
        self.vehicle_type = vehicle_type
        self.vehicle_size = vehicle_size

        # If a string is passed, the corresponding driving cycle is retrieved

        if self.cycle_name:
            try:
                cycle, gradient = get_standard_driving_cycle_and_gradient(
                    vehicle_type, vehicle_size, self.cycle_name
                )

            except KeyError as err:
                raise KeyError(
                    "The driving cycle specified could not be found."
                ) from err
        elif isinstance(cycle, np.ndarray):
            self.cycle_name = "custom"

        else:
            raise TypeError("The format of the driving cycle is not valid.")

        self.cycle = cycle
        self.gradient = gradient

        if "Micro" in vehicle_size:
            idx = vehicle_size.index("Micro")
            self.cycle[:, idx] = np.clip(self.cycle[:, idx], 0, 90)

        assert len(self.cycle) == len(
            self.gradient
        ), "The length of the driving cycle and the gradient must be the same."

        # Unit conversion km/h to m/s
        self.velocity = np.where(np.isnan(self.cycle), 0, (self.cycle * 1000) / 3600)
        self.velocity = self.velocity[:, None, None, None, :]

        # Model acceleration as difference in velocity between
        # time steps (1 second)
        # Zero at first value
        self.acceleration = np.zeros_like(self.velocity)
        self.acceleration[1:-1] = (self.velocity[2:, ...] - self.velocity[:-2, ...]) / 2

    def get_country_temperature(self, country):
        """
        Retrieves mothly average temperature
        :type country: country for which to retrieve temperature values
        :return:
        """

        with open(DATA_DIR / MONTHLY_AVG_TEMP) as f:
            reader = csv.reader(f, delimiter=";")
            for row in reader:
                if row[2] == country:
                    return np.array([float(i) for i in row[3:]])

        print(
            f"Could not find monthly average temperature series for {country}. "
            f"Uses those for CH instead."
        )

        with open(DATA_DIR / MONTHLY_AVG_TEMP) as f:
            reader = csv.reader(f, delimiter=";")
            for row in reader:
                if row[2] == "CH":
                    return np.array([int(i) for i in row[3:]])

    def calculate_hvac_energy(
        self,
        hvac_power,
        battery_cooling_unit,
        battery_heating_unit,
        ambient_temp,
        indoor_temp,
    ) -> tuple[Any, Any, Any, Any]:

        if ambient_temp is not None:
            ambient_temp = np.resize(ambient_temp, (12,))
        else:
            ambient_temp = self.get_country_temperature(self.country)

        # use ambient temperature if provided, otherwise
        # monthly temperature average (12 values)
        # relation between ambient temperature
        # and HVAC power required
        # from https://doi.org/10.1016/j.energy.2018.12.064
        amb_temp_data_points = np.array([-30, -20, -10, 0, 10, 20, 30, 40])
        pct_power_HVAC = np.array([0.95, 0.54, 0.29, 0.13, 0.04, 0.08, 0.45, 0.7])

        # Heating power as long as ambient temperature, in W
        # is below the comfort indoor temperature
        p_heating = (
            np.where(
                ambient_temp < indoor_temp,
                np.interp(ambient_temp, amb_temp_data_points, pct_power_HVAC),
                0,
            ).mean()
            * hvac_power
        ).values

        # Cooling power as long as ambient temperature, in W
        # is above the comfort indoor temperature
        p_cooling = (
            np.where(
                ambient_temp >= indoor_temp,
                np.interp(ambient_temp, amb_temp_data_points, pct_power_HVAC),
                0,
            ).mean()
            * hvac_power
        ).values

        # We want to add power draw for battery cooling
        # and battery heating

        # battery cooling occurring above 20C, in W
        p_battery_cooling = np.where(ambient_temp > 20, _(battery_cooling_unit), 0)
        p_battery_cooling = p_battery_cooling.mean(-1)

        # battery heating occurring below 5C, in W
        p_battery_heating = np.where(ambient_temp < 5, _(battery_heating_unit), 0)
        p_battery_heating = p_battery_heating.mean(-1)

        return p_cooling, p_heating, p_battery_cooling, p_battery_heating

    def aux_energy_per_km(
        self,
        aux_power: Union[xr.DataArray, np.array],
        efficiency: Union[xr.DataArray, np.array],
        hvac_power: Union[xr.DataArray, np.array] = None,
        battery_cooling_unit: Union[xr.DataArray, np.array] = None,
        battery_heating_unit: Union[xr.DataArray, np.array] = None,
        heat_pump_cop_cooling: Union[xr.DataArray, np.array] = None,
        heat_pump_cop_heating: Union[xr.DataArray, np.array] = None,
        cooling_consumption: Union[xr.DataArray, np.array] = None,
        heating_consumption: Union[xr.DataArray, np.array] = None,
        ambient_temp: float = None,
        indoor_temp: float = 20.0,
    ) -> Union[float, np.ndarray]:
        """
        Calculate energy used other than motive energy per km driven.

        :param aux_power: Total power needed for auxiliaries, heating, and cooling (W)
        :param efficiency: Efficiency of electricity generation (dimensionless, between 0.0 and 1.0).
                Battery electric vehicles should have efficiencies of one here, as we account for
                battery efficiencies elsewhere.
        :returns: total auxiliary energy in kJ/km

        """

        _o = lambda x: np.where(x == 0, 1, x)

        if hvac_power is not None:

            (
                p_cooling,
                p_heating,
                p_battery_cooling,
                p_battery_heating,
            ) = self.calculate_hvac_energy(
                hvac_power=hvac_power,
                battery_cooling_unit=battery_cooling_unit,
                battery_heating_unit=battery_heating_unit,
                ambient_temp=ambient_temp,
                indoor_temp=indoor_temp,
            )

            aux_energy = (
                aux_power
                + (p_cooling / _o(heat_pump_cop_cooling) * cooling_consumption)
                + (p_heating / _o(heat_pump_cop_heating) * heating_consumption)
                + p_battery_cooling
                + p_battery_heating
            ).T.values * np.ones_like(self.velocity)

            return aux_energy

        _c = lambda x: x.values if isinstance(x, xr.DataArray) else x

        # Provide energy in kJ / km (1 J = 1 Ws)
        auxiliary_energy = (
            _c(aux_power).T[None, ...]  # Watts
            * np.ones_like(self.velocity)
            * 1000  # m / km
            / 1000  # 1 / (J / kJ)
        )

        efficiency = _c(efficiency)[..., None]
        efficiency[efficiency == 0] = 1

        return auxiliary_energy / _o(efficiency).T

    def convert_to_xr(self, data):

        return xr.DataArray(
            data,
            dims=["second", "value", "year", "powertrain", "size", "parameter"],
            coords={
                "second": range(0, data.shape[0]),
                "value": range(0, data.shape[1]),
                "year": range(0, data.shape[2]),
                "powertrain": range(0, data.shape[3]),
                "size": range(0, data.shape[4]),
                "parameter": [
                    "rolling resistance",
                    "air resistance",
                    "gradient resistance",
                    "kinetic energy",
                    "motive energy at wheels",
                    "motive energy",
                    "negative motive energy",
                    "recuperated energy",
                    "power load",
                    "auxiliary energy",
                    "transmission efficiency",
                    "engine efficiency",
                    "velocity",
                ],
            },
        )

    def motive_energy_per_km(
        self,
        driving_mass: Union[xr.DataArray, np.array],
        rr_coef: Union[xr.DataArray, np.array],
        drag_coef: Union[xr.DataArray, np.array],
        frontal_area: Union[xr.DataArray, np.array],
        electric_motor_power: Union[xr.DataArray, np.array],
        engine_power: Union[xr.DataArray, np.array],
        recuperation_efficiency: Union[xr.DataArray, np.array],
        aux_power: Union[xr.DataArray, np.array],
        engine_efficiency: Union[xr.DataArray, np.array],
        transmission_efficiency: Union[xr.DataArray, np.array],
        battery_charge_eff: Union[xr.DataArray, np.array],
        battery_discharge_eff: Union[xr.DataArray, np.array],
        fuel_cell_system_efficiency: Union[xr.DataArray, np.array] = None,
        hvac_power: Union[xr.DataArray, np.array] = None,
        battery_cooling_unit: Union[xr.DataArray, np.array] = None,
        battery_heating_unit: Union[xr.DataArray, np.array] = None,
        heat_pump_cop_cooling: Union[xr.DataArray, np.array] = None,
        heat_pump_cop_heating: Union[xr.DataArray, np.array] = None,
        cooling_consumption: Union[xr.DataArray, np.array] = None,
        heating_consumption: Union[xr.DataArray, np.array] = None,
        ambient_temp: float = None,
        indoor_temp: float = 20.0,
    ) -> DataArray:
        """
        Calculate energy used and recuperated for a given vehicle per km driven.

        :param driving_mass: Mass of vehicle (kg)
        :param rr_coef: Rolling resistance coefficient (dimensionless, between 0.0 and 1.0)
        :param drag_coef: Aerodynamic drag coefficient (dimensionless, between 0.0 and 1.0)
        :param frontal_area: Frontal area of vehicle (m2)
        :param sizes: size classes of the vehicles
        :param electric_motor_power: Electric motor power (watts). Optional.
        :returns: net motive energy (in kJ/km)

        Power to overcome rolling resistance is calculated by:

        .. math::

            g v M C_{r}

        where :math:`g` is 9.81 (m/s2), :math:`v` is velocity (m/s), :math:`M` is mass (kg),
        and :math:`C_{r}` is the rolling resistance coefficient (dimensionless).

        Power to overcome air resistance is calculated by:

        .. math::

            \frac{1}{2} \rho_{air} v^{3} A C_{d}


        where :math:`\rho_{air}` is 1.225 (kg/m3), :math:`v` is velocity (m/s), :math:`A` is frontal area (m2), and :math:`C_{d}`
        is the aerodynamic drag coefficient (dimensionless).

        """

        _c = lambda x: x.values if isinstance(x, xr.DataArray) else x
        _o = lambda x: np.where((x == 0) | (x == np.nan), 1, x)

        # Calculate the energy used for each second of the drive cycle
        ones = np.ones_like(self.velocity)

        # Resistance from the tire rolling: rolling resistance coefficient * driving mass * 9.81
        rolling_resistance = _c((driving_mass * rr_coef * 9.81).T) * ones

        # Resistance from the drag: frontal area * drag coefficient * air density * 1/2 * velocity^2
        air_resistance = _c((frontal_area * drag_coef * self.rho_air / 2).T) * np.power(
            self.velocity, 2
        )

        # Resistance from road gradient: driving mass * 9.81 * sin(gradient)
        gradient_resistance = _c((driving_mass * 9.81).T) * np.sin(
            np.nan_to_num(self.gradient)[:, None, None, None, :]
        )

        # Inertia: driving mass * acceleration
        inertia = self.acceleration * _c(driving_mass).T

        total_resistance = (
            rolling_resistance + air_resistance + gradient_resistance + inertia
        )

        engine_power = xr.where(engine_power == 0, 1, engine_power)
        engine_efficiency = xr.where(engine_efficiency == 0, 1, engine_efficiency)
        transmission_efficiency = xr.where(
            transmission_efficiency == 0, 1, transmission_efficiency
        )
        recuperation_efficiency = xr.where(
            recuperation_efficiency == 0, 1, recuperation_efficiency
        )

        if fuel_cell_system_efficiency is None:
            fuel_cell_system_efficiency = np.array([1.0])

        fuel_cell_system_efficiency = xr.where(
            fuel_cell_system_efficiency == 0, 1, fuel_cell_system_efficiency
        )

        motive_energy_at_wheels = xr.where(total_resistance < 0, 0, total_resistance)

        if fuel_cell_system_efficiency is None:
            fuel_cell_system_efficiency = np.array([1.0])

        motive_energy = (
            motive_energy_at_wheels
            / _o(_c(engine_efficiency.T))[None, ...]
            / _o(_c(transmission_efficiency.T))[None, ...]
            / _o(_c(fuel_cell_system_efficiency.T))[None, ...]
        )

        negative_motive_energy = xr.where(total_resistance > 0, 0, total_resistance)
        recuperated_energy = (
            negative_motive_energy
            * _c(recuperation_efficiency).T[None, ...]
            * _c(battery_charge_eff).T[None, ...]
            * _c(battery_discharge_eff).T[None, ...]
            * (_c(electric_motor_power).T[None, ...] > 0)
        )

        auxiliary_energy = self.aux_energy_per_km(
            aux_power,
            engine_efficiency,
            hvac_power,
            battery_cooling_unit,
            battery_heating_unit,
            heat_pump_cop_cooling,
            heat_pump_cop_heating,
            cooling_consumption,
            heating_consumption,
            ambient_temp,
            indoor_temp,
        )
        auxiliary_energy *= self.velocity > 0

        all_arrays = np.concatenate(
            [
                _(rolling_resistance),
                _(air_resistance),
                _(gradient_resistance),
                _(inertia),
                _(motive_energy_at_wheels),
                _(motive_energy),
                _(negative_motive_energy),
                _(recuperated_energy),
                _((motive_energy + recuperated_energy) / (_c(engine_power).T * 1000)),
                _(auxiliary_energy),
                _(__(_c(transmission_efficiency).T) * np.ones_like(motive_energy)),
                _(__(_c(engine_efficiency).T) * np.ones_like(motive_energy)),
                _(self.velocity * np.ones_like(motive_energy)),
            ],
            axis=-1,
        )

        all_arrays[..., :-5] /= 1000
        all_arrays[..., -4] /= 1000
        all_arrays[..., :-5] *= _(self.velocity)
        all_arrays[..., -5] = np.clip(all_arrays[..., -5], 0, 1)
        all_arrays[..., 4] = np.where(
            all_arrays[..., 4] > _(engine_power).T,
            _(engine_power).T,
            all_arrays[..., 4],
        )
        all_arrays[..., 5] = np.where(
            all_arrays[..., 5] > _(engine_power).T,
            _(engine_power).T,
            all_arrays[..., 5],
        )
        all_arrays[..., 7] = np.where(
            all_arrays[..., 7] < _(electric_motor_power).T * -1,
            _(electric_motor_power).T * -1,
            all_arrays[..., 7],
        )

        return self.convert_to_xr(all_arrays).fillna(0)
