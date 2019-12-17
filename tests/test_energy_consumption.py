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
    aux = ecm.aux_energy_per_km(93.75)
    aux /= 0.29
    motive = ecm.motive_energy_per_km(
        driving_mass=1879,
        rr_coef=0.01,
        drag_coef=0.3,
        frontal_area=2.4,
        ttw_efficiency=0.232).sum()
    assert int(motive + aux) > 2500
    assert int(motive + aux) < 2700