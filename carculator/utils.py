import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from . import DATA_DIR
import itertools
import csv
import numpy as np
from pathlib import Path

pd.options.mode.chained_assignment = None

REGION_MAPPING_FILEPATH = DATA_DIR / "regionmappingH12.csv"
IAM_ELEC_MARKETS = DATA_DIR / "electricity_markets.csv"
IEA_DIESEL_SHARE = DATA_DIR / "diesel_share_oecd.csv"


def extract_biofuel_shares_from_IAM(
    model, fp, IAM_region, years, allocate_all_synfuel=False
):
    """
    This function extracts biofuel shares from a IAM file provided.

    :param model:
    :param fp: file path to IAM file
    :type fp: str
    :param IAM_region: IAM region for which to extract the biofuel shares
    :type IAM_region: str
    :param years: the list of years for which to extract biofuel shares
    :param allocate_all_synfuel: Temporary workaround. If True, then all synfuel in the transport sector
    is allocated to passenger cars.
    :return: a dictionary that contains fuel types as keys and lists of fuel shares as values
    """

    if model == "remind":
        df = pd.read_csv(fp, delimiter=";", index_col=["Region", "Variable", "Unit"]).drop(
            columns=["Model", "Scenario"]
        )

    if model == "image":
        df = pd.read_excel(fp, index_col=[2, 3, 4]).drop(columns=["Model", "Scenario"])

    df = df.reset_index()
    df = df.loc[df["Region"] == IAM_region]
    df = df.loc[:, : str(2050)]
    df["Variable"] = df["Variable"].str.replace("|", "-")

    if model == "remind":
        if allocate_all_synfuel:

            # get shares of synthetic fuel
            df_total = df.loc[df["Variable"] == "FE-Transport-Pass-Road-LDV-Liquids"]
            df_total.index = df.loc[df["Variable"] == "SE-Liquids-Hydrogen"].index
            share = np.clip(
                (
                    df.loc[df["Variable"] == "SE-Liquids-Hydrogen", "2005":].divide(
                        df_total.loc[:, "2005":], axis=0
                    )
                ).values,
                0,
                1,
            )

            var = ["FE-Transport-Liquids-Oil", "FE-Transport-Liquids-Biomass"]

            df_liquids = df.loc[df["Variable"].isin(var), :]
            df_liquids.iloc[:, 3:] /= df_liquids.iloc[:, 3:].sum(axis=0)

            df_liquids.loc[:, "2005":] *= 1 - share
            to_append = [IAM_region, "FE-Transport-Liquids-Hydrogen", "EJ/yr"] + share[
                0
            ].tolist()
            a_series = pd.Series(to_append, index=df_liquids.columns)
            df_liquids = df_liquids.append(a_series, ignore_index=True)

        else:

            var = [
                "FE-Transport-Liquids-Oil",
                "FE-Transport-Liquids-Biomass",
                "FE-Transport-Liquids-Hydrogen",
            ]

            df_liquids = df.loc[df["Variable"].isin(var), :]


            if len(df_liquids) == 0:
                df_liquids = pd.DataFrame(0,
                                          columns=["Region", "Variable", "Unit"] + [str(i) for i in list(range(2005, 2050, 5))],
                                          index=range(len(var))
                                          )
                df_liquids["Region"] = IAM_region
                df_liquids["Variable"] = var
                df_liquids["Unit"] = "EJ/yr"

                # Set 100% fossil if missing data
                df_liquids.iloc[0, 3:] = 1

            df_liquids.iloc[:, 3:] /= df_liquids.iloc[:, 3:].sum(axis=0)



        var = [
            "FE-Transport-Gases-Non-Biomass",
            "FE-Transport-Gases-Biomass"
            ]

        df_gas = df.loc[df["Variable"].isin(var), :]
        df_gas.iloc[:, 3:] /= df_gas.iloc[:, 3:].sum(axis=0)

        if len(df_gas) == 0:
            df_gas = pd.DataFrame(0,
                                      columns=["Region", "Variable", "Unit"] + [str(i) for i in list(range(2005, 2055, 5))],
                                      index=range(0, len(var))
                                      )
            df_gas["Region"] = IAM_region
            df_gas["Variable"] = var
            df_gas["Unit"] = "EJ/yr"

            # Set 100% fossil if missing data
            df_gas.iloc[0, 3:] = 1

        d_map_fuels = {
            "FE-Transport-Liquids-Oil": "liquid - fossil",
            "FE-Transport-Liquids-Biomass": "liquid - biomass",
            "FE-Transport-Liquids-Hydrogen": "liquid - synfuel",
            "FE-Transport-Gases-Non-Biomass": "gas - fossil",
            "FE-Transport-Gases-Biomass": "gas - biomass",
        }

    if model == "image":
        var = [
            "Final Energy-Transportation-Freight-Liquids-Oil",
            "Final Energy-Transportation-Freight-Liquids-Biomass",
        ]

        df_liquids = df.loc[df["Variable"].isin(var), :]

        if len(df_liquids) == 0:
            df_liquids = pd.DataFrame(0,
                                      columns=["Region", "Variable", "Unit"] + [str(i) for i in list(range(2005, 2050, 5))],
                                      index=range(len(var))
                                      )
            df_liquids["Region"] = IAM_region
            df_liquids["Variable"] = var
            df_liquids["Unit"] = "EJ/yr"

            # Set 100% fossil if missing data
            df_liquids.iloc[0, 3:] = 1

        df_liquids.iloc[:, 3:] /= df_liquids.iloc[:, 3:].sum(axis=0)

        var = ["Final Energy-Transportation-Freight-Gases"]
        df_gas = df.loc[df["Variable"].isin(var), :]

        if len(df_gas) == 0:
            df_gas = pd.DataFrame(0,
                                      columns=["Region", "Variable", "Unit"] + [str(i) for i in list(range(2005, 2050, 5))],
                                      index=range(len(var))
                                      )
            df_gas["Region"] = IAM_region
            df_gas["Variable"] = var
            df_gas["Unit"] = "EJ/yr"

            # Set 100% fossil if missing data
            df_gas.iloc[0, 3:] = 1

        df_gas.iloc[:, 3:] /= df_gas.iloc[:, 3:].sum(axis=0)

        d_map_fuels = {
            "Final Energy-Transportation-Freight-Liquids-Oil": "liquid - fossil",
            "Final Energy-Transportation-Freight-Liquids-Biomass": "liquid - biomass",
            "Final Energy-Transportation-Freight-Gases": "gas - fossil",
        }

    new_df = pd.concat([df_liquids, df_gas])
    new_df["Variable"] = new_df["Variable"].map(d_map_fuels)
    new_df.columns = (
        new_df.columns[:3].tolist() + new_df.columns[3:].astype(int).tolist()
    )

    print(new_df)

    new_df = new_df.rename(columns={"Variable": "fuel_type"})
    new_df = new_df.groupby("fuel_type").sum()

    print(new_df)

    arr = (
        new_df.to_xarray()
        .to_array()
        .interp(variable=years, kwargs={"fill_value": "extrapolate"})
    )

    arr = np.clip(arr, 0, 1)
    return arr

def get_iam_electricity_market_labels(model):
    """
    Loads a csv file into a dictionary. This dictionary contains labels of electricity markets
    in the IAM.

    :return: dictionary that contains market names equivalence
    :rtype: dict
    """

    d = dict()
    with open(IAM_ELEC_MARKETS) as f:
        reader = csv.reader(f, delimiter=";")
        for row in reader:
            if row[0] == model:
                d[row[1]] = row[2]
    return d

def extract_electricity_mix_from_IAM_file(model, fp, IAM_region, years):
    """
    This function extracts electricity mixes from a IAM file provided.

    :param model:
    :param fp: file path to IAM file
    :type fp: str
    :param IAM_region: IAM region for which to extract the electricity mix
    :type IAM_region: str
    :param years: the list of years for which to extract electricity mixes
    :type years: list
    :return: list of lists of electricity mixes, that can be consumed by `InventnoryCalculation`
    :rtype: list
    """

    electricity_markets = get_iam_electricity_market_labels(model)
    rev_tech = {v: k for k, v in electricity_markets.items()}

    if model == "remind":

        df = pd.read_csv(
            fp, delimiter=";", index_col=["Region", "Variable", "Unit"]
        ).drop(columns=["Model", "Scenario"])

        d_var = {
            "Biomass IGCC CCS": "Biomass CCS",
            "Biomass IGCC": "Biomass",
            "Biomass CHP": "Biomass",
            "Coal IGCC": "Coal",
            "Coal IGCC CCS": "Coal CCS",
            "Coal PC": "Coal",
            "Coal PC CCS": "Coal CCS",
            "Coal CHP": "Coal",
            "Gas CCS": "Gas CCS",
            "Gas CC": "Gas",
            "Gas OC": "Gas",
            "Gas CHP": "Gas",
            "Hydrogen": "Hydrogen",
            "Oil ST": "Oil",
            "Nuclear": "Nuclear",
            "Geothermal": "Geothermal",
            "Hydro": "Hydro",
            "Solar CSP": "Solar",
            "Solar PV": "Solar",
            "Wind": "Wind",
        }

    if model == "image":
        df = pd.read_excel(fp, index_col=[2, 3, 4]).drop(columns=["Model", "Scenario"])

        d_var = {
            "Biomass CHP": "Biomass",
            "Biomass CHP CCS": "Biomass CCS",
            "Biomass IGCC CCS": "Biomass CCS",
            "Biomass IGCC": "Biomass",
            "Biomass ST": "Biomass",
            "Coal PC": "Coal",
            "Coal IGCC": "Coal",
            "Coal IGCC CCS": "Coal CCS",
            "Coal CHP": "Coal",
            "Coal CHP CCS": "Coal",
            "Gas OC": "Gas",
            "Gas CC": "Gas",
            "Gas CHP": "Gas",
            "Gas CC CCS": "Gas CCS",
            "Gas CHP CCS": "Gas CCS",
            "Geothermal": "Geothermal",
            "Hydro": "Hydro",
            "Nuclear": "Nuclear",
            "Oil CC CCS": "Coal CCS",
            "Oil CHP CCS": "Coal CCS",
            "Oil ST": "Oil",
            "Oil CC": "Oil",
            "Oil CHP": "Oil",
            "Solar CSP": "Solar",
            "Solar PV Centralized": "Solar",
            "Solar PV Residential": "Solar",
            "Wind Onshore": "",
            "Wind Offshore": "Wind",
        }

    df = df.reset_index()
    df = df.loc[df["Region"] == IAM_region]
    df = df.loc[:, : str(2050)]

    print(electricity_markets)
    print(len(electricity_markets))
    print(df.head())
    print(df["Variable"].unique())
    df = df.loc[df["Variable"].isin(electricity_markets.values())]
    df["Variable"] = df["Variable"].map(rev_tech)
    df.columns = df.columns[:3].tolist() + df.columns[3:].astype(int).tolist()
    df.iloc[:, 3:] /= df.iloc[:, 3:].sum(axis=0)

    df["Variable"] = df["Variable"].map(d_var)

    list_tech = [
        "Hydro",
        "Nuclear",
        "Gas",
        "Solar",
        "Wind",
        "Biomass",
        "Coal",
        "Oil",
        "Geothermal",
        "Waste",
        "Biogas CCS",
        "Biomass CCS",
        "Coal CCS",
        "Gas CCS",
        "Wood CCS",
    ]

    for row in [i for i in list_tech if i not in df["Variable"].unique()]:
        df.loc[df.index[-1] + 1] = [IAM_region, row, "EJ/yr"] + [0] * (df.shape[1] - 3)

    df = df.groupby("Variable").sum().loc[list_tech]

    arr = (
        df.to_xarray()
        .to_array()
        .interp(variable=years, kwargs={"fill_value": "extrapolate"})
        .values
    )
    arr = np.clip(arr, 0, 1)
    arr /= np.sum(arr, axis=1)[:, None]

    return arr

def create_fleet_composition_from_IAM_file(
    fp
):
    """
    This function creates a consumable fleet composition array from a CSV file.
    The array returned is consumed by `InventoryCalculation`.

    :param fp: Path file path
    :type fp: Path
    :return: fleet composition array
    :rtype: xarray.DataArray
    """

    if isinstance(fp, str):
        fp = Path(fp)

    if not fp.is_file():
        raise FileNotFoundError("Could not locate {}".format(fp))

    # Read the fleet composition CSV file
    df = pd.read_csv(fp, delimiter=",")
    df = df.fillna(0)

    # Filter out unecessary columns
    df = df[["year", "IAM_region", "powertrain", "size", "vintage_year", "vintage_demand_vkm"]]

    #df_gr = df.groupby(["IAM_region", "powertrain", "size", "year", "vintage_year"]).sum()
    #df_gr = df_gr.groupby(level=[0, 1, 3]).apply(lambda x: x / float(x.sum()))

    #df = df_gr.reset_index()

    # Turn the dataframe into a pivot table
    df = df.pivot_table(
        index=["IAM_region", "powertrain", "size", "vintage_year"], columns=["year"], aggfunc=np.sum
    )["vintage_demand_vkm"]

    # xarray.DataArray is returned
    return df.to_xarray().fillna(0).to_array().round(3)


def build_fleet_array(fp, scope):
    """
    Receives a file path that points to a CSV file that contains the fleet composition
    Checks that the fleet composition array is valid.

    Specifically:

    * the years specified in the fleet must be present in scope["year"]
    * the powertrains specified in the fleet must be present in scope["powertrain"]
    * the sizes specified in the fleet must be present in scope["size"]
    * the sum for each year-powertrain-size set must equal 1

    :param scope:
    :param fp: filepath to an array that contains fleet composition
    :type fp: str
    :return array: fleet composition array
    :rtype array: xarray.DataArray
    """
    arr = pd.read_csv(fp, delimiter=";", header=0, index_col=[0, 1, 2])
    arr = arr.fillna(0)
    arr.columns = [int(c) for c in arr.columns]

    new_cols = [c for c in scope["year"] if c not in arr.columns]
    arr[new_cols] = pd.DataFrame([[0] * len(new_cols)], index=arr.index)

    a = [scope["powertrain"]] + [scope["size"]] + [scope["year"]]

    for row in [i for i in list(itertools.product(*a)) if i not in arr.index]:
        arr.loc[row] = 0

    array = arr.to_xarray()
    array = array.rename(
        {"level_0": "powertrain", "level_1": "size", "level_2": "vintage_year"}
    )

    if not set(list(array.data_vars)).issubset(scope["year"]):
        raise ValueError("The fleet years differ from {}".format(scope["year"]))

    if set(scope["year"]) != set(array.coords["vintage_year"].values.tolist()):
        raise ValueError(
            "The list of vintage years differ from {}.".format(self.scope["year"])
        )

    if not set(array.coords["powertrain"].values.tolist()).issubset(
        scope["powertrain"]
    ):
        raise ValueError(
            "The fleet powertrain list differs from {}".format(scope["powertrain"])
        )

    if not set(array.coords["size"].values.tolist()).issubset(scope["size"]):
        raise ValueError("The fleet size list differs from {}".format(scope["size"]))

    return array.to_array().fillna(0)
