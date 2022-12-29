import numpy as np

from carculator_utils.energy_consumption import EnergyConsumptionModel


def test_acceleration():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size=["Medium"],
        cycle="WLTC",
        gradient="WLTC",
    )
    assert ecm.cycle.shape == ecm.gradient.shape


def test_aux_power():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size=["Medium"],
        cycle="WLTC",
        gradient="WLTC",
    )
    aux = ecm.aux_energy_per_km(
        aux_power=np.array(100),
        efficiency=np.array(0.5)
    )
    assert isinstance(aux, np.ndarray)
    assert np.squeeze(aux).shape == np.squeeze(ecm.cycle).shape
    assert np.allclose(aux, np.repeat(100/0.5, 1801))


def test_motive_energy():
    ecm = EnergyConsumptionModel(
        vehicle_type="car",
        vehicle_size=["Medium"],
        cycle="WLTC",
        gradient="WLTC",
    )
    aux = np.full(1801, 94) / 0.29 / 1000
    energy = ecm.motive_energy_per_km(
        driving_mass=np.array(1736),
        rr_coef=np.array(0.0093),
        drag_coef=np.array(0.296),
        frontal_area=np.array(2.247826),
        electric_motor_power=np.array(0),
        engine_power=np.array(172 * 0.74),
        recuperation_efficiency=np.array(0),
        aux_power=np.array(100),
        engine_efficiency=np.array(0.3),
        transmission_efficiency=np.array(0.9),
        battery_charge_eff=np.array(1),
        battery_discharge_eff=np.array(1),
    )

    distance = energy.sel(parameter="velocity").sum(dim="second") / 1000
    motive = energy.sel(parameter="motive energy")

    print(distance)
    print(motive.sum(dim="second") / distance / 36000 * 100)

    assert motive.sum(dim="second") / distance / 36000 * 100 < 10.0
    assert motive.sum(dim="second") / distance / 36000 * 100 > 5.0
