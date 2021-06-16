from carculator import *

def test_monte_carlo():
    cip = CarInputParameters()
    cip.stochastic(5)
    dcts, array = fill_xarray_from_input_parameters(cip,
                                                    scope={"powertrain": ["ICEV-d", "BEV"],
                                                           "size": ["Lower medium"],
                                                           "year": [2020]})
    cm = CarModel(array, cycle='WLTC 3.4')
    cm.set_all()
    ic = InventoryCalculation(cm.array)
    results = ic.calculate_impacts()

    assert len(cm.array.value.values) == 5
    assert len(results.value.values) == 5