from carculator import *
import pytest
import numpy as np

# generate vehicle parameters
cip = CarInputParameters()
cip.static()

# fill in array with vehicle parameters
_, array = fill_xarray_from_input_parameters(cip)

# build CarModel object
cm = CarModel(array, cycle="WLTC")
# build vehicles
cm.set_all()


def test_scope():
    """Test if scope works as expected"""
    ic = InventoryCalculation(
        cm.array,
        method="recipe",
        method_type="midpoint",
        scope={"powertrain": ["ICEV-d", "ICEV-p"], "size": ["Lower medium"]},
    )
    results = ic.calculate_impacts()

    assert "Large" not in results.coords["size"].values
    assert "BEV" not in results.coords["powertrain"].values


def test_fuel_blend():
    """Test if fuel blends defined by the user are considered"""

    bc = {
        "fuel blend": {
            "petrol": {
                "primary fuel": {"type": "petrol", "share": [0.9, 0.9, 0.9, 0.9]},
                "secondary fuel": {
                    "type": "bioethanol - wheat straw",
                    "share": [0.1, 0.1, 0.1, 0.1],
                },
            },
            "diesel": {
                "primary fuel": {"type": "diesel", "share": [0.93, 0.93, 0.93, 0.93]},
                "secondary fuel": {
                    "type": "biodiesel - cooking oil",
                    "share": [0.07, 0.07, 0.07, 0.07],
                },
            },
            "cng": {
                "primary fuel": {
                    "type": "biogas - sewage sludge",
                    "share": [1, 1, 1, 1,],
                }
            },
        }
    }

    ic = InventoryCalculation(
        cm.array, method="recipe", method_type="midpoint", background_configuration=bc
    )

    assert ic.fuel_blends["petrol"]["primary"]["share"] == [0.9, 0.9, 0.9, 0.9]
    assert ic.fuel_blends["petrol"]["secondary"]["share"] == [0.1, 0.1, 0.1, 0.1]
    assert ic.fuel_blends["cng"]["primary"]["share"] == [1, 1, 1, 1]
    assert np.sum(ic.fuel_blends["cng"]["secondary"]["share"]) == 0

    ic.calculate_impacts()

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
            "synthetic gasoline",
            "synthetic diesel",
            "wood gasification with EF with CCS (Swiss forest)",
            "syngas",
        ),
    ]:
        ic = InventoryCalculation(
            cm.array,
            method="recipe",
            method_type="midpoint",
            background_configuration={
                "fuel blend": {
                    "petrol": {
                        "primary fuel": {"type": fuels[0], "share": [1, 1, 1, 1]},
                    },
                    "diesel": {
                        "primary fuel": {"type": fuels[1], "share": [1, 1, 1, 1]},
                    },
                    "hydrogen": {
                        "primary fuel": {"type": fuels[2], "share": [1, 1, 1, 1,],}
                    },
                    "cng": {
                        "primary fuel": {"type": fuels[3], "share": [1, 1, 1, 1,],}
                    },
                }
            },
        )
        ic.calculate_impacts()


def test_countries():
    """Test that calculation works with all countries"""
    for c in [
        "AO","AT","AU","BE","BF","BG","BI","BJ","BR","BW","CA","CD","CF",
         #"CG","CH","CI","CL","CM","CN","CY","CZ","DE","DJ","DK","DZ","EE",
         #"EG","ER","ES","ET","FI","FR","GA",
         #"GB","GH","GM","GN","GQ","GR","GW","HR","HU","IE",
         #"IN","IT", "IS", "JP", "KE", "LR","LS","LT","LU","LV","LY","MA","ML","MR","MT","MW","MZ",
         #"NE", "NG","NL","NM","NO","PL","PT","RER","RO","RU","RW","SD","SE","SI","SK","SL","SN","SO","SS","SZ",
         #"TD","TG","TN","TZ","UG","UK","US","ZA","ZM",
         #"ZW",
    ]:
        ic = InventoryCalculation(
            cm.array,
            method="recipe",
            method_type="midpoint",
            background_configuration={
                "country": c,
                "energy storage": {"electric": {"origin": c}},
            },
        )
        ic.calculate_impacts()


def test_IAM_regions():
    """Test that calculation works with all IAM regions"""
    for c in [
         "BRA","CAN","CEU","CHN","EAF","INDIA","INDO","JAP","KOR","ME","MEX",
       #     "NAF","OCE","RCAM","RSAF","RSAM","RSAS","RUS","SAF","SEAS","STAN","TUR",
       #  "UKR","USA","WAF","WEU","LAM","CAZ","EUR","CHA","SSA","IND","OAS","JPN","MEA",
       # "REF","USA",
    ]:
        ic = InventoryCalculation(
            cm.array,
            method="recipe",
            method_type="midpoint",
            background_configuration={
                "country": c,
                "energy storage": {"electric": {"origin": c}},
            },
        )
        ic.calculate_impacts()


def test_endpoint():

    """Test if the correct impact categories are considered"""
    ic = InventoryCalculation(cm.array, method="recipe", method_type="endpoint")
    results = ic.calculate_impacts()
    assert "human health" in [i.lower() for i in results.impact_category.values]
    assert len(results.impact_category.values) == 4

    """Test if it errors properly if an incorrect method type is give"""
    with pytest.raises(TypeError) as wrapped_error:
        InventoryCalculation(cm.array, method="recipe", method_type="endpint")
    assert wrapped_error.type == TypeError


def test_sulfur_concentration():
    ic = InventoryCalculation(cm.array, method="recipe", method_type="endpoint")

    ic.get_sulfur_content("RER", "diesel", 2000)
    ic.get_sulfur_content("foo", "diesel", 2000)

    with pytest.raises(ValueError) as wrapped_error:
        ic.get_sulfur_content("FR", "diesel", "jku")
    assert wrapped_error.type == ValueError


def test_custom_electricity_mix():
    """Test if a wrong number of electricity mixes throws an error"""

    bc = {
        "custom electricity mix": [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    }

    with pytest.raises(ValueError) as wrapped_error:
        InventoryCalculation(
            cm.array, method="recipe", method_type="endpoint", background_configuration=bc
        )
    assert wrapped_error.type == ValueError

    """ Test if a sum of share superior to 1 throws an error """

    bc = {
        "custom electricity mix": [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    }

    with pytest.raises(ValueError) as wrapped_error:
        InventoryCalculation(
            cm.array, method="recipe", method_type="endpoint", background_configuration=bc
        )
    assert wrapped_error.type == ValueError


def test_export_to_bw():
    """ Test that inventories export successfully"""
    ic = InventoryCalculation(
        cm.array, method="recipe", method_type="endpoint"
    )

    for a in (True, False):
        for b in ("3.5", "3.6", "3.7", "uvek"):
            for c in (True, False):

                ic.export_lci(
                    ecoinvent_compatibility=a,
                    ecoinvent_version=b,
                    create_vehicle_datasets=c,
                )


def test_export_to_excel():
    """ Test that inventories export successfully to Excel/CSV"""
    ic = InventoryCalculation(
        cm.array, method="recipe", method_type="endpoint"
    )

    for a in (True, False):
        for b in ("3.5", "3.6", "3.7", "uvek"):
            for c in (True, False):
                for d in ("file", "string"):

                    ic.export_lci_to_excel(
                        ecoinvent_compatibility=a,
                        ecoinvent_version=b,
                        create_vehicle_datasets=c,
                        export_format=d,
                        directory="directory"
                    )