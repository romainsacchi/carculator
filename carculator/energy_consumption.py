import numexpr as ne
import numpy as np
import xarray

from .driving_cycles import get_standard_driving_cycle


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

        - :func:`~energy_consumption.EnergyConsumptionModel.aux_energy_per_km` calculates the energy needed to power auxiliary services
        - :func:`~energy_consumption.EnergyConsumptionModel.motive_energy_per_km` calculates the energy needed to move the vehicle over 1 km

    Acceleration is calculated as the difference between velocity at t_2 and velocity at t_0, divided by 2.
    See for example: http://www.unece.org/fileadmin/DAM/trans/doc/2012/wp29grpe/WLTP-DHC-12-07e.xls

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series
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

    def __init__(self, cycle, rho_air=1.204, gradient=None):
        # If a string is passed, the corresponding driving cycle is retrieved
        if isinstance(cycle, str):
            try:
                self.cycle_name = cycle
                cycle = get_standard_driving_cycle(cycle)

            except KeyError:
                raise "The driving cycle specified could not be found."
        elif isinstance(cycle, np.ndarray):
            self.cycle_name = "custom"
            pass
        else:
            raise "The format of the driving cycle is not valid."

        self.cycle = cycle
        self.rho_air = rho_air

        if gradient is not None:
            try:
                assert isinstance(gradient, np.ndarray)
            except AssertionError:
                raise AssertionError(
                    "The type of the gradient array is not valid. Required: numpy.ndarray."
                )
            try:
                assert len(gradient) == len(self.cycle)
            except AssertionError:
                raise AssertionError(
                    "The length of the gradient array does not equal the length of the driving cycle."
                )
            self.gradient = gradient
        else:
            self.gradient = np.zeros_like(cycle)

        # Unit conversion km/h to m/s
        self.velocity = (cycle * 1000) / 3600

        # Model acceleration as difference in velocity between time steps (1 second)
        # Zero at first value
        self.acceleration = np.zeros_like(self.velocity)
        self.acceleration[1:-1] = (self.velocity[2:] - self.velocity[:-2]) / 2

    def aux_energy_per_km(self, aux_power, efficiency=1):
        """
        Calculate energy used other than motive energy per km driven.

        :param aux_power: Total power needed for auxiliaries, heating, and cooling (W)
        :type aux_power: int
        :param efficiency: Efficiency of electricity generation (dimensionless, between 0.0 and 1.0).
                Battery electric vehicles should have efficiencies of one here, as we account for
                battery efficiencies elsewhere.
        :type efficiency: float

        :returns: total auxiliary energy in kJ/km
        :rtype: float

        """
        # Provide energy in kJ / km (1 J = 1 Ws)
        auxiliary_energy = (
            aux_power  # Watt
            * self.velocity.size  # Number of seconds -> Ws -> J
            / self.velocity.sum()  # m/s * 1s = m -> J/m
            * 1000  # m / km
            / 1000  # 1 / (J / kJ)
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
        """
        Calculate energy used and recuperated for a given vehicle per km driven.

        :param driving_mass: Mass of vehicle (kg)
        :type driving_mass: int
        :param rr_coef: Rolling resistance coefficient (dimensionless, between 0.0 and 1.0)
        :type rr_coef: float
        :param drag_coef: Aerodynamic drag coefficient (dimensionless, between 0.0 and 1.0)
        :type drag_coef: float
        :param frontal_area: Frontal area of vehicle (m2)
        :type frontal_area: float
        :param ttw_efficiency: Efficiency of translating potential energy into motion (dimensionless, between 0.0 and 1.0)
        :type ttw_efficiency: float
        :param recuperation_efficiency: Fraction of energy that can be recuperated (dimensionless, between 0.0 and 1.0). Optional.
        :type recuperation_efficiency: float
        :param motor_power: Electric motor power (watts). Optional.
        :type motor_power: int

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

        :returns: net motive energy (in kJ/km)
        :rtype: float

        """

        # Convert to km; velocity is m/s, times 1 second
        # Distance WLTC 3.2 = 4.75 km
        distance = self.velocity.sum() / 1000

        # Total power required at the wheel to meet acceleration requirement,
        # and overcome air and rolling resistance.
        # This number is generally positive (power is needed), but can be negative
        # if the vehicle is decelerating.
        # Power is in watts (kg m2 / s3)

        # We opt for simpler variable names to be accepted by `numexpr`
        ones = np.ones_like(self.velocity)
        dm = _(driving_mass)
        rr = _(rr_coef)
        fa = _(frontal_area)
        dc = _(drag_coef)
        v = self.velocity
        a = self.acceleration
        g = self.gradient
        rho_air = self.rho_air
        ttw_eff = _(ttw_efficiency)
        mp = _(motor_power)
        re = _(recuperation_efficiency)

        # rolling resistance + air resistance + kinetic energy + gradient resistance
        total_force = np.float16(
            ne.evaluate(
                "(ones * dm * rr * 9.81) + (v ** 2 * fa * dc * rho_air / 2) + (a * dm) + (dm * 9.81 * sin(g))"
            )
        )

        tv = ne.evaluate("total_force * v")

        # Can only recuperate when power is less than zero, limited by recuperation efficiency
        # Motor power in kW, other power in watts

        recuperated_power = ne.evaluate(
            "where(tv < (-1000 * mp), (-1000 * mp), where(tv>0, 0, tv))"
        )
        # braking_power = pd.w - recuperated_power

        # self.recuperated_power = recuperated_power/distance/1000
        # self.braking_power = braking_power/distance/1000
        # self.power_rolling_resistance = pa.r / distance / 1000
        # self.power_aerodynamic = pa.a / distance / 1000
        # self.power_kinetic = pa.k / distance / 1000
        # self.total_power = pa.w / distance / 1000

        # t_e = ne.evaluate("where(total_force<0, 0, tv)") #
        # t_e = np.where(total_force<0, 0, tv)

        # results = ne.evaluate(
        #    "((where(total_force<0, 0, tv) / (distance * 1000)) + (recuperated_power / distance / 1000))/ ttw_eff"
        # )

        return tv, recuperated_power, distance
