"""
.. module: energy_consumption.py

"""

from .driving_cycles import get_standard_driving_cycle
import numpy as np
import pandas as pd
import xarray


def _(o):
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(o, (np.ndarray, xarray.DataArray)):
        return np.expand_dims(o, -1)
    else:
        return o


class EnergyConsumptionModel:
    """
    Calculate energy consumption of a vehicle for a given driving cycle and vehicle parameters.

    Based on a selected driving cycle, this class calculates the acceleration needed and provides
    two methods:
    - aux_energy_per_km() calculates the energy needed to power auxiliary services
    - motive_energy_per_km() calculates the energy needed to move the vehicle over 1 km

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series
    :param rho_air: Mass per unit volume of air. Set to (1.225 kg/m3) by default.
    :type rho_air: float

    :ivar rho_air: Mass per unit volume of air.
    :vartype rho_air: float
    :ivar velocity: Time series of speed values, in meters per second.
    :vartype velocity: numpy.ndarray
    :ivar acceleration: Time series of acceleration, calculated as increment in velocity per interval of 1 second,
        in meter per second^2.
    :vartype acceleration: numpy.ndarray


    """

    def __init__(self, cycle, rho_air=1.225):
        # If a string is passed, the corresponding driving cycle is retrieved
        if isinstance(cycle, str):
            cycle = get_standard_driving_cycle(cycle).values
        # If the parameter is not a string nor a Pandas Series, this is an issue.
        elif not isinstance(cycle, pd.Series):
            try:
                raise TypeError('Invalid format for "cycle" parameter. str or panda.Series expected.')
            except TypeError as error:
                print('Caught this error: ' + repr(error))


        self.rho_air = rho_air
        # Unit conversion km/h to m/s
        self.velocity = cycle * 1000 / 3600

        # Model acceleration as difference in velocity between time steps (1 second)
        # Zero at first value
        self.acceleration = np.zeros_like(self.velocity)
        self.acceleration[1:] = self.velocity[1:] - self.velocity[:-1]

    def aux_energy_per_km(self, aux_power, efficiency=1):
        """Calculate energy used other than motive energy per km driven.

        :param aux_power: Total power needed for auxiliaries, heating, and cooling (W)
        :type aux_power: int
        :param efficiency: Efficiency of electricity generation (dimensionless, between 0.0 and 1.0).
                Battery electric vehicles should have efficiencies of one here, as we account for
                battery efficiencies elsewhere.
        :type: efficiency: float

        :returns: total auxiliary energy in kJ/km
        :rtype: float

        """
        # Provide energy in kJ / km (1 J = 1 Ws)
        auxiliary_energy = (
            aux_power               # Watt
            * self.velocity.size    # Number of seconds -> Ws -> J
            / self.velocity.sum()   # m/s * 1s = m -> J/m
            * 1000                  # m / km
            / 1000                  # 1 / (J / kJ)
        )

        return auxiliary_energy / efficiency

    def motive_energy_per_km(
        self,
        driving_mass,
        rr_coef,
        drag_coef,
        frontal_area,
        ttw_efficiency,
        recuperation_efficiency=0,
        motor_power=0,
    ):
        """Calculate energy used and recuperated for a given vehicle per km driven.

        :param driving_mass: Mass of vehicle (kg)
        :type driving_mass: int
        :param rr_coef: Rolling resistance coefficient (dimensionless, between 0.0 and 1.0)
        :type rr_coef: float
        :param drag_coef: Aerodynamic drag coefficient (dimensionless, between 0.0 and 1.0)
        :type drag_coef: float
        :param frontal_area: Frontal area of vehicle (m2)
        :type frontal_area: float
        :param ttw_efficiency: Efficiency of translating potential energy into motion (dimensionless,
                between 0.0 and 1.0)
        :type ttw_efficiency: float
        :param recuperation_efficiency: Fraction of energy that can be recuperated (dimensionless, between 0.0 and 1.0).
                Optional.
        :type recuperation_efficiency: float
        :param motor_power: Electric motor power (watts). Optional.
        :type motor_power: int


        Power to overcome rolling resistance is calculated by:

        .. math::

            g v M C_{r}

        where :math:`g` is 9.81 (m/s2), *v* is velocity (m/s), *M* is mass (kg), and :math:`C_{r}` is rolling
        resistance coefficient (dimensionless).

        Power to overcome air resistance is calculated by:

        .. math::

            \frac{1}{2} \rho_{air} v^{3} A C_{d}

        where :math:`\rho_{air}` is 1.225 (kg/m3), *v* is velocity (m/s), *A* is frontal area (m2), and :math:`C_{d}` is
        aerodynamic drag coefficient (dimensionless).

        Returns net motive energy in kJ/km

        """

        # Total power required at the wheel to meet acceleration requirement,
        # and overcome air and rolling resistance.
        # This number is generally positive (power is needed), but can be negative
        # if the vehicle is decelerating.
        # Power is in watts (kg m2 / s3)
        power_rolling_resistance = np.ones_like(self.velocity) * _(driving_mass) * _(rr_coef) * 9.81
        power_aerodynamic = (self.velocity ** 2 * _(frontal_area) * _(drag_coef)
               * self.rho_air / 2)
        power_kinetic = self.acceleration * self.velocity * _(driving_mass)

        net_power = (
            np.where(power_kinetic > 0, power_kinetic, 0)
            + power_rolling_resistance
            + power_aerodynamic
        )
        total_power = power_kinetic + power_rolling_resistance + power_aerodynamic

        # Can only recuperate when power is less than zero, limited by recuperation efficiency
        # Motor power in kW, other power in watts
        self.recuperated_power = (
            np.clip(total_power, -1000 * _(motor_power), 0) * _(recuperation_efficiency)
        )

        return (
            (net_power + self.recuperated_power) # watt
            / self.velocity.sum()                # m/s -> Ws/m -> J/m
            * 1000                               # m / km -> J/km
            / 1000                               # 1 / (J / kJ) -> kJ/km
            / _(ttw_efficiency)
        )
