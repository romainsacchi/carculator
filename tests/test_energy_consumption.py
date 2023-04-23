import numpy as np
from carculator_utils.energy_consumption import EnergyConsumptionModel


def test_acceleration():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size=["Medium"],
        cycle="WLTC",
        gradient="WLTC",
        powertrains=["ICEV-d"],
    )
    assert ecm.cycle.shape == ecm.gradient.shape


def test_aux_power():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size=["Medium"],
        cycle="WLTC",
        gradient="WLTC",
        powertrains=["ICEV-d"],
    )
    aux = ecm.aux_energy_per_km(aux_power=np.array([100]), efficiency=np.array([0.5]))
    assert isinstance(aux, np.ndarray)
    assert aux.reshape(ecm.cycle.shape).shape == ecm.cycle.shape
    assert np.allclose(aux, np.repeat(100 / 0.5, 1801))
