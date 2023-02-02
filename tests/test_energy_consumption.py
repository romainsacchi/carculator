import numpy as np
from carculator_utils.energy_consumption import EnergyConsumptionModel


def test_acceleration():
    ecm = EnergyConsumptionModel(
        vehicle_type="car", vehicle_size=["Medium"], cycle="WLTC", gradient="WLTC"
    )
    assert ecm.cycle.shape == ecm.gradient.shape


def test_aux_power():
    ecm = EnergyConsumptionModel(
        vehicle_type="car", vehicle_size=["Medium"], cycle="WLTC", gradient="WLTC"
    )
    aux = ecm.aux_energy_per_km(aux_power=np.array([100]), efficiency=np.array([0.5]))
    assert isinstance(aux, np.ndarray)
    assert aux.reshape(ecm.cycle.shape).shape == ecm.cycle.shape
    assert np.allclose(aux, np.repeat(100 / 0.5, 1801))


def test_motive_energy():
    ecm = EnergyConsumptionModel(
        vehicle_type="car", vehicle_size=["Medium"], cycle="WLTC", gradient="WLTC"
    )

    array = ecm.motive_energy_per_km(
        driving_mass=np.array([1879]),
        rr_coef=np.array([0.0093]),
        drag_coef=np.array([0.3]),
        frontal_area=np.array([2.4]),
        electric_motor_power=np.array([0.0]),
        engine_power=np.array([100.0]),
        recuperation_efficiency=np.array([0.0]),
        aux_power=np.array([100]),
        engine_efficiency=np.array([0.29]),
        transmission_efficiency=np.array([0.9]),
        battery_charge_eff=np.array([0.9]),
        battery_discharge_eff=np.array([0.9]),
    )

    dist = array.sel(parameter="velocity").sum(dim="second") / 1000

    motive = (
        array.sel(
            parameter=["motive energy", "auxiliary energy", "recuperated energy"]
        ).sum(dim=["second", "parameter"])
        / dist
    ).T

    print(motive)

    assert motive / 36000 * 100 < 10.0
    assert motive / 36000 * 100 > 5.0
