"""
energy_consumption.py contains the class EnergyConsumption Model which exposes two methods:
one for calculating the auxiliary energy needs, and another one for calculating the motive
energy needs.
"""

from typing import Any, Union

import numexpr as ne
import numpy as np
import xarray as xr

from .driving_cycles import get_standard_driving_cycle


def _(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, -1)

    return obj


class EnergyConsumptionModel:
    """
    Calculate energy consumption of a vehicle for a given driving cycle and vehicle parameters.

    Based on a selected driving cycle, this class calculates the acceleration needed and provides
    two methods:

        - :func:`~energy_consumption.EnergyConsumptionModel.aux_energy_per_km` calculates the energy needed to power auxiliary services
        - :func:`~energy_consumption.EnergyConsumptionModel.motive_energy_per_km` calculates the energy needed to move the vehicle over 1 km

    Acceleration is calculated as the difference between velocity at t_2 and velocity at t_0, divided by 2.
    See for example: http://www.unece.org/fileadmin/DAM/trans/doc/2012/wp29grpe/WLTP-DHC-12-07e.xls

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: np.ndarray
    :param rho_air: Mass per unit volume of air. Set to (1.225 kg/m3) by default.
    :type rho_air: float
    :param gradient: Road gradient per second of driving, in degrees. None by default. Should be passed as an array of
                    length equal to the length of the driving cycle.
    :type gradient: numpy.ndarray

    :ivar rho_air: Mass per unit volume of air. Value of 1.204 at 23C (test temperature for WLTC).
    :vartype rho_air: float
    :ivar velocity: Time series of speed values, in meters per second.
    :vartype velocity: numpy.ndarray
    :ivar acceleration: Time series of acceleration, calculated as increment in velocity per interval of 1 second,
        in meter per second^2.
    :vartype acceleration: numpy.ndarray


    """

    def __init__(
        self,
        cycle: Union[str, np.ndarray],
        rho_air: float = 1.204,
        gradient: Union[str, np.ndarray] = None,
    ) -> None:
        # If a string is passed, the corresponding driving cycle is retrieved
        if isinstance(cycle, str):
            try:
                self.cycle_name = cycle
                cycle = get_standard_driving_cycle(cycle)

            except KeyError as err:
                raise KeyError(
                    "The driving cycle specified could not be found."
                ) from err
        elif isinstance(cycle, np.ndarray):
            self.cycle_name = "custom"

        else:
            raise TypeError("The format of the driving cycle is not valid.")

        self.cycle = cycle
        self.rho_air = rho_air

        if gradient is not None:
            try:
                assert isinstance(gradient, np.ndarray)
            except AssertionError as err:
                raise AssertionError(
                    "The type of the gradient array is not valid. Required: numpy.ndarray."
                ) from err
            try:
                assert len(gradient) == len(self.cycle)
            except AssertionError as err:
                raise AssertionError(
                    "The length of the gradient array does not equal the length of the driving cycle."
                ) from err
            self.gradient = gradient
        else:
            self.gradient = np.zeros_like(cycle)

    def aux_energy_per_km(
        self, aux_power: Union[float, np.ndarray], efficiency: float = 1.0
    ) -> Union[float, np.ndarray]:
        """
        Calculate energy used other than motive energy per km driven.

        :param aux_power: Total power needed for auxiliaries, heating, and cooling (W)
        :param efficiency: Efficiency of electricity generation (dimensionless, between 0.0 and 1.0).
                Battery electric vehicles should have efficiencies of one here, as we account for
                battery efficiencies elsewhere.
        :returns: total auxiliary energy in kJ/km

        """
        # Unit conversion km/h to m/s
        velocity = (self.cycle * 1000) / 3600

        # Provide energy in kJ / km (1 J = 1 Ws)
        auxiliary_energy = (
            aux_power  # Watt
            * velocity.size  # Number of seconds -> Ws -> J
            / velocity.sum()  # m/s * 1s = m -> J/m
            * 1000  # m / km
            / 1000  # 1 / (J / kJ)
        )

        return auxiliary_energy / efficiency

    def motive_energy_per_km(
        self,
        driving_mass: Union[float, np.ndarray, xr.DataArray],
        rr_coef: Union[float, np.ndarray, xr.DataArray],
        drag_coef: Union[float, np.ndarray, xr.DataArray],
        frontal_area: Union[float, np.ndarray, xr.DataArray],
        sizes: Union[str, np.ndarray, xr.DataArray],
        motor_power: Union[float, np.ndarray, xr.DataArray] = 0,
    ) -> tuple[Union[float, Any], Any, Union[float, Any]]:
        """
        Calculate energy used and recuperated for a given vehicle per km driven.

        :param driving_mass: Mass of vehicle (kg)
        :param rr_coef: Rolling resistance coefficient (dimensionless, between 0.0 and 1.0)
        :param drag_coef: Aerodynamic drag coefficient (dimensionless, between 0.0 and 1.0)
        :param frontal_area: Frontal area of vehicle (m2)
        :param sizes: size classes of the vehicles
        :param motor_power: Electric motor power (watts). Optional.
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

        # correct the driving cycle speed profile
        # if a micro car is present
        # as those can only drive at 90 km/h max
        self.cycle = np.repeat(self.cycle.reshape(1, -1), len(sizes), axis=0).T

        if "Micro" in sizes:
            idx = sizes.tolist().index("Micro")
            self.cycle[idx] = np.clip(self.cycle[idx], 0, 90)

        # Unit conversion km/h to m/s
        velocity = (self.cycle * 1000) / 3600

        # Model acceleration as difference in velocity between time steps (1 second)
        # Zero at first value
        acceleration = np.zeros_like(velocity)
        acceleration[1:-1] = (velocity[2:] - velocity[:-2]) / 2

        # Convert to km; velocity is m/s, times 1 second
        # Distance WLTC 3.2 = 4.75 km
        distance = velocity.sum(axis=0) / 1000

        # Total power required at the wheel to meet acceleration requirement,
        # and overcome air and rolling resistance.
        # This number is generally positive (power is needed), but can be negative
        # if the vehicle is decelerating.
        # Power is in watts (kg m2 / s3)

        ones = np.ones_like(velocity).T[:, None, None, None]
        driving_mass = _(driving_mass)
        rr_coeff = _(rr_coef)
        frontal_area = _(frontal_area)
        drag_coef = _(drag_coef)
        motor_power = _(motor_power)
        velocity = velocity.T[:, None, None, None]
        acceleration = acceleration.T[:, None, None, None]
        gradient = self.gradient
        self.rho_air = self.rho_air

        # rolling resistance + air resistance + kinetic energy + gradient resistance
        rolling_resistance = ones * driving_mass * rr_coeff * 9.81
        air_resistance = (
            np.power(velocity, 2) * frontal_area * drag_coef * self.rho_air / 2
        )

        kinetic_energy = (acceleration * driving_mass) + (
            driving_mass * 9.81 * np.sin(gradient)
        )
        total_force = rolling_resistance + air_resistance + kinetic_energy

        total_force_velocity = total_force * velocity

        # Can only recuperate when power is less than zero, limited by recuperation efficiency
        # Motor power in kW, other power in watts

        recuperated_power = ne.evaluate(
            "where(total_force_velocity < (-1000 * motor_power), (-1000 * motor_power), where(total_force_velocity>0, 0, total_force_velocity))"
        )

        return total_force_velocity, recuperated_power, distance
