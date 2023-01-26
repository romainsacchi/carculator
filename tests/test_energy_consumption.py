import numpy as np

from carculator_utils.energy_consumption import EnergyConsumptionModel


def test_acceleration():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size="Medium",
        cycle="WLTC",
        gradient="WLTC",
    )
    assert ecm.cycle.shape == ecm.gradient.shape


def test_aux_power():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size="Medium",
        cycle="WLTC",
        gradient="WLTC",
    )
    aux = ecm.auxiliary_power_per_km(
        aux_power=100,
        efficiency=0.5
    )
    assert isinstance(aux, np.ndarray)
    assert aux.shape == ecm.cycle.shape
    assert np.allclose(aux, np.repeat(100/0.5, 1801))


def test_motive_energy():
    ecm = EnergyConsumptionModel("WLTC")
    aux = np.full(1801, 94) / 0.29 / 1000
    motive, recup, dist = ecm.motive_energy_per_km(
        driving_mass=1879,
        rr_coef=0.01,
        drag_coef=0.3,
        frontal_area=2.4,
        sizes=["Medium"],
    )

    motive = np.clip(motive / 1000, 0, None)

    assert int(motive.sum() - recup.sum() + aux.sum()) / 0.24 / dist > 2500
    assert int(motive.sum() - recup.sum() + aux.sum()) / 0.24 / dist < 2700
