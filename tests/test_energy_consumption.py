from carculator.energy_consumption import EnergyConsumptionModel


def test_acceleration():
    ecm = EnergyConsumptionModel('WLTC')
    assert ecm.acceleration.sum() == 0
    assert int(ecm.velocity.sum()) == 23262

def test_aux_power():
    ecm = EnergyConsumptionModel('WLTC')
    assert int(ecm.aux_energy_per_km(300)) == 23

def test_motive_energy():
    ecm = EnergyConsumptionModel('WLTC')
    motive = ecm.motive_energy_per_km(
        driving_mass=1800,
        rr_coef=0.05,
        drag_coef=0.01,
        frontal_area=2,
        ttw_efficiency=0.1).sum()
    assert int(motive) == 9615



