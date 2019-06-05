from .driving_cycles import get_standard_driving_cycle
import numpy as np
import xarray


def _(o):
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(o, (np.ndarray, xarray.DataArray)):
        return np.expand_dims(o, -1)
    else:
        return o


class EnergyConsumptionModel:
    def __init__(self, cycle, rho_air=1.225):
        """Calculate energy consumption of a vehicle for a given driving cycle and vehicle parameters.

        Input parameters:

        * ``cycle``: 1-d Numpy array of second-by-second speeds (km/h)
        * ``rho_air``: Density of air (kg/m3). Optional.

        Creates ``self.velocity`` (m/s) and ``self.acceleration`` (m/s2)

        """
        if isinstance(cycle, str):
            cycle = get_standard_driving_cycle(cycle).values

        self.rho_air = rho_air
        # Unit conversion km/h to m/s
        self.velocity = cycle * 1000 / 3600

        # Model acceleration as difference in velocity between time steps (1 second)
        # Zero at first value
        self.acceleration = np.zeros_like(self.velocity)
        self.acceleration[1:] = self.velocity[1:] - self.velocity[:-1]

    def aux_energy_per_km(self, aux_power, efficiency=1):
        """Calculate energy used other than motive energy.

        Input parameters:

        * ``aux_power``: Total power needed for auxiliaries, heating, and cooling (W)
        * ``efficiency``: Efficiency of electricity generation (dimensionless)

        Battery electric vehicles should have efficiencies of one here, as we
        account for battery efficiencies elsewhere.

        Returns total auxiliary energy in kJ/km

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
        """Calculate energy used and recuperated for a given vehicle.

        Input parameters:

        * ``driving_mass``: Mass of vehicle (kg)
        * ``rr_coef``: Rolling resistance coefficient (dimensionless)
        * ``drag_coef``: Aerodynamic drag coefficient (dimensionless)
        * ``frontal_area``: Frontal area of vehicle (m2)
        * ``ttw_efficiency``: Efficiency of translating potential energy into motion (dimensionless)
        * ``recuperation_efficiency``: Fraction of energy that can be recuperated (dimensionless). Optional.
        * ``motor_power``: Electric motor power (watts). Optional.

        Power to overcome rolling resistance is calculated by:

        .. math::

            g v M C_{r}

        where :math:`g` is 9.81 (m/s2), *v* is velocity (m/s), *M* is mass (kg), and :math:`C_{r}` is rolling resistance coefficient (dimensionless).

        Power to overcome air resistance is calculated by:

        .. math::

            \frac{1}{2} \rho_{air} v^{3} A C_{d}

        where :math:`\rho_{air}` is 1.225 (kg/m3), *v* is velocity (m/s), *A* is frontal area (m2), and :math:`C_{d}` is aerodynamic drag coefficient (dimensionless).

        Returns net motive energy in kJ/km

        """

        # Total power required at the wheel to meet acceleration requirement,
        # and overcome air and rolling resistance.
        # This number is generally positive (power is needed), but can be negative
        # if the vehicle is decelerating.
        # Power is in Watts (kg m2 / s3)
        positive_acceleration = np.where(self.acceleration > 0, self.acceleration, 0)
        power = (
            positive_acceleration * self.velocity * _(driving_mass) # Kinetic
            + np.ones_like(self.velocity) * _(driving_mass) * _(rr_coef) * 9.81  # Rolling resistance
            + (self.velocity ** 2 * _(frontal_area) * _(drag_coef) # Air resistance
               * self.rho_air / 2)
        )

        # Can only recuperate when power is less than zero, limited by recuperation efficiency
        self.recuperated_power = (
            np.clip(power, -1000 * _(motor_power), 0) * _(recuperation_efficiency)
        )

        return (
            (power + self.recuperated_power) # Watt
            / self.velocity.sum()            # m/s -> Ws/m -> J/m
            * 1000                           # m / km -> J/km
            / 1000                           # 1 / (J / kJ) -> kJ/km
            / _(ttw_efficiency)
        )
