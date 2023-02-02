import numpy as np
import pytest

from carculator import (
    CarInputParameters,
    CarModel,
    InventoryCar,
    fill_xarray_from_input_parameters,
)

# generate vehicle parameters
cip = CarInputParameters()
cip.static()

# fill in array with vehicle parameters
scope = {"powertrain": ["ICEV-d", "ICEV-p", "BEV"], "size": ["Medium"]}
_, array = fill_xarray_from_input_parameters(cip, scope=scope)

# build CarModel object
cm = CarModel(array, cycle="WLTC")
# build vehicles
cm.set_all()


def test_scope():
    """Test if scope works as expected"""
    ic = InventoryCar(
        cm,
        method="recipe",
        indicator="midpoint",
    )
    results = ic.calculate_impacts()

    assert "Large" not in results.coords["size"].values
    assert "FCEV" not in results.coords["powertrain"].values


def test_plausibility_of_GWP():
    """Test if GWP scores make sense"""

    for method in ["recipe", "ilcd"]:
        ic = InventoryCar(cm, method=method, indicator="midpoint")
        results = ic.calculate_impacts()

        if method == "recipe":
            m = "climate change"
        else:
            m = "climate change - climate change total"

        gwp_icev = results.sel(
            impact_category=m,
            powertrain=["ICEV-d", "ICEV-p"],
            value=0,
            year=2020,
            size="Medium",
        )

        # Are the medium ICEVs between 0.28 and 0.35 kg CO2-eq./vkm?

        if method == "recipe":
            assert (gwp_icev.sum(dim="impact") > 0.24).all() and (
                gwp_icev.sum(dim="impact") < 0.31
            ).all()

            # Are the medium ICEVs direct emissions between 0.13 and  0.18 kg CO2-eq./vkm?
            assert (gwp_icev.sel(impact="direct - exhaust") > 0.13).all() and (
                gwp_icev.sel(impact="direct - exhaust") < 0.19
            ).all()

        # Are the ICEVs glider emissions between 0.04 and 0.05 kg CO2-eq./vkm?
        # print(m)
        # print(gwp_icev.sel(impact="glider"))
        # assert (gwp_icev.sel(impact="glider") > 0.04).all() and (
        #    gwp_icev.sel(impact="glider") < 0.05
        # ).all()

        # Is the GWP score for batteries of Medium BEVs between 0.025 and 0.035 kg Co2-eq./vkm?
        gwp_bev = results.sel(
            impact_category=m, powertrain="BEV", value=0, year=2020, size="Medium"
        )

        assert (gwp_bev.sel(impact="energy storage") > 0.02).all() and (
            gwp_bev.sel(impact="energy storage") < 0.04
        ).all()

        assert gwp_bev.sel(impact="direct - exhaust") == 0


def test_fuel_blend():
    """Test if fuel blends defined by the user are considered"""

    bc = {
        "petrol": {
            "primary": {
                "type": "petrol",
                "share": [0.9, 0.9, 0.9, 0.9, 0.9, 0.9],
            },
        },
        "diesel": {
            "primary": {
                "type": "diesel",
                "share": [0.93, 0.93, 0.93, 0.93, 0.93, 0.93],
            },
        },
        "hydrogen": {
            "primary": {
                "type": "electrolysis",
                "share": [0.9, 0.9, 0.9, 0.9, 0.9, 0.9],
            },
        },
        "cng": {
            "primary": {
                "type": "biogas - sewage sludge",
                "share": [
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                ],
            }
        },
    }

    cm = CarModel(array, cycle="WLTC", fuel_blend=bc)
    cm.set_all()

    assert np.array_equal(
        cm.fuel_blend["petrol"]["primary"]["share"],
        np.array(
            [
                0.9,
                0.9,
                0.9,
                0.9,
                0.9,
                0.9,
            ]
        ),
    )

    assert np.array_equal(
        cm.fuel_blend["diesel"]["primary"]["share"],
        np.array(
            [
                0.93,
                0.93,
                0.93,
                0.93,
                0.93,
                0.93,
            ]
        ),
    )
    assert np.array_equal(
        cm.fuel_blend["cng"]["primary"]["share"], np.array([1, 1, 1, 1, 1, 1])
    )
    assert np.allclose(np.sum(cm.fuel_blend["cng"]["secondary"]["share"]), np.zeros(6))

    for fuels in [
        ("petrol", "diesel", "electrolysis", "cng"),
        (
            "bioethanol - wheat straw",
            "biodiesel - palm oil",
            "smr - natural gas",
            "biogas - sewage sludge",
        ),
        (
            "bioethanol - forest residues",
            "biodiesel - rapeseed oil",
            "smr - natural gas with CCS",
            "biogas - biowaste",
        ),
        (
            "bioethanol - maize starch",
            "biodiesel - cooking oil",
            "wood gasification with EF with CCS",
            "biogas - biowaste",
        ),
        (
            "bioethanol - sugarbeet",
            "biodiesel - algae",
            "atr - biogas",
            "biogas - biowaste",
        ),
        (
            "synthetic gasoline - energy allocation",
            "synthetic diesel - energy allocation",
            "wood gasification with EF with CCS",
            "syngas",
        ),
    ]:
        bc["petrol"]["primary"]["type"] = fuels[0]
        bc["diesel"]["primary"]["type"] = fuels[1]
        bc["hydrogen"]["primary"]["type"] = fuels[2]
        bc["cng"]["primary"]["type"] = fuels[3]

        cm = CarModel(array, cycle="WLTC", fuel_blend=bc)
        cm.set_all()
        ic = InventoryCar(cm)
        ic.calculate_impacts()


def test_countries():
    """Test that calculation works with all countries"""
    for c in [
        "AO",
        "AT",
        "AU",
        "BE",
    ]:
        ic = InventoryCar(
            cm,
            method="recipe",
            indicator="midpoint",
            background_configuration={
                "country": c,
                "energy storage": {"electric": {"type": "NMC-622"}, "origin": c},
            },
        )
        ic.calculate_impacts()


def test_endpoint():
    """Test if the correct impact categories are considered"""
    ic = InventoryCar(cm, method="recipe", indicator="endpoint")
    results = ic.calculate_impacts()
    assert "human health" in [i.lower() for i in results.impact_category.values]
    assert len(results.impact_category.values) == 4
    #
    #     """Test if it errors properly if an incorrect method type is give"""
    with pytest.raises(ValueError) as wrapped_error:
        ic = InventoryCar(cm, method="recipe", indicator="endpint")
        ic.calculate_impacts()
    assert wrapped_error.type == ValueError


def test_sulfur_concentration():
    ic = InventoryCar(cm, method="recipe", indicator="endpoint")
    ic.get_sulfur_content("RER", "diesel", 2000)
    ic.get_sulfur_content("foo", "diesel", 2000)

    with pytest.raises(ValueError) as wrapped_error:
        ic.get_sulfur_content("FR", "diesel", "jku")
    assert wrapped_error.type == ValueError


def test_custom_electricity_mix():
    """Test if a wrong number of electricity mixes throws an error"""

    # Passing four mixes instead of 6
    mix_1 = np.zeros((5, 15))
    mix_1[:, 0] = 1

    mixes = [mix_1]

    for mix in mixes:
        with pytest.raises(ValueError) as wrapped_error:
            InventoryCar(
                cm,
                method="recipe",
                indicator="endpoint",
                background_configuration={"custom electricity mix": mix},
            )
        assert wrapped_error.type == ValueError


def test_export_to_bw():
    """Test that inventories export successfully"""
    ic = InventoryCar(cm, method="recipe", indicator="endpoint")
    #

    for b in (
        "3.5",
        "3.6",
        "3.7",
        "3.8",
    ):
        for c in (True, False):
            ic.export_lci(
                ecoinvent_version=b,
            )


def test_export_to_excel():
    """Test that inventories export successfully to Excel/CSV"""
    ic = InventoryCar(cm, method="recipe", indicator="endpoint")

    for b in ("3.5", "3.6", "3.7", "3.7.1", "3.8"):
        for s in ("brightway2", "simapro"):
            for d in ("file", "bw2io"):
                #
                ic.export_lci(
                    ecoinvent_version=b,
                    format=d,
                    software=s,
                    directory="directory",
                )
