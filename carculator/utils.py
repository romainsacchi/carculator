import pandas as pd
import pycountry
from . import DATA_DIR
import itertools
import csv
import numpy as np
from pathlib import Path
import bisect

pd.options.mode.chained_assignment = None

REGION_MAPPING_FILEPATH = DATA_DIR / "regionmappingH12.csv"
IAM_ELEC_MARKETS = DATA_DIR / "electricity_markets.csv"
IEA_DIESEL_SHARE = DATA_DIR / "diesel_share_oecd.csv"


def extract_biofuel_shares_from_IAM(
    model, fp, IAM_region, years, allocate_all_synfuel=False
):
    """
    This function extracts biofuel shares from a IAM file provided.

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
                                      columns=["Region", "Variable", "Unit"] + [str(i) for i in list(range(2005, 2050, 5))],
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

    new_df = new_df.rename(columns={"Variable": "fuel_type"})
    new_df = new_df.groupby("fuel_type").sum()

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

def create_fleet_composition_from_REMIND_file(
    fp, remind_region, fleet_year=2020, normalized=True
):
    """
    This function creates a consumable fleet composition array from a CSV file.
    The array returned is consumed by `InventoryCalculation`.

    :param fp: Path file path
    :type fp: Path
    :param remind_region: REMIND region for which to extract fleet composition
    :type remind_region: str
    :param fleet_year: the year for which to extract fleet composition
    :type fleet_year: int
    :return: fleet composition array
    :rtype: xarray.DataArray
    """

    if isinstance(fp, str):
        fp = Path(fp)

    if not fp.is_file():
        raise FileNotFoundError("Could not locate {}".format(fp))

    # Load diesel shares from IEA on OECD countries
    share_diesel = pd.read_csv(
        IEA_DIESEL_SHARE, sep=";", usecols=range(1, 10), index_col=[0]
    )
    share_diesel.columns = [int(c) for c in share_diesel.columns]
    shares = pd.concat(
        [
            share_diesel,
            pd.DataFrame(columns=range(2019, 2051)),
            pd.DataFrame(columns=range(2000, 2011)),
        ]
    )
    shares = shares.T.fillna(method="ffill")
    shares = shares.fillna(method="bfill")

    # Read the fleet composition CSV file
    df = pd.read_csv(fp, delimiter=",")
    df = df.fillna(0)

    def get_iso_2_to_REMIND_map():
        """ Generate a dictionary of shape {`iso alpha-2 country code`:`REMIND region`}"""
        with open(REGION_MAPPING_FILEPATH) as f:
            f.readline()
            csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
            d = {x[1]: x[2] for x in csv_list}
        return d

    # Generate a dictionary of shape {`iso alpha-3 country code`:`iso alpha-2 country code`}
    d_iso3_iso2 = {
        c: pycountry.countries.get(alpha_3=c).alpha_2 for c in df["iso"].unique()
    }

    df["iso_2"] = df["iso"].map(d_iso3_iso2)
    d_iso_2_remind = get_iso_2_to_REMIND_map()
    df["REMIND_region"] = df["iso_2"].map(d_iso_2_remind)

    # Filter out rows for which no REMIND region equivalence can be found.
    if len(df.loc[df["REMIND_region"].isnull(), "iso_2"].unique()) > 0:
        print(
            "No REMIND region could be mapped to the following ISO alpha2 country codes: {}".format(
                df.loc[df["REMIND_region"].isnull(), "iso_2"].unique()
            )
        )
        print("They will be skipped.")
        df = df.loc[~df["REMIND_region"].isnull()]

    # Dictionary map of shape {`REMIND car size`:`carculator car size`}
    d_map_sizes = {
        "Compact Car": "Lower medium",
        "Large Car and SUV": "Large",
        "Mini Car": "Mini",
        "Subcompact Car": "Small",
        "Van": "Van",
        "Multipurpose Vehicle": "Van",
        "Midsize Car": "Medium",
        "Large Car": "Large",
        "Light Truck and SUV": "SUV",
    }

    # Dictionary map of shape {`REMIND car powertrain`:`carculator car powertrain`}
    d_map_tech = {
        "BEV": "BEV",
        "FCEV": "FCEV",
        "Hybrid Electric": "PHEV",
        "Hybrid Liquids": "HEV",
        "Liquids": "ICEV",
        "NG": "ICEV-g",
    }

    df["powertrain"] = df["technology"].map(d_map_tech)
    df["size"] = df["vehicle_type"].map(d_map_sizes)

    # Clean up
    df["variable"] = df["variable"].str.replace(r"\D", "").astype(int)

    # Rename `variable`, since it is a built-in name used by `xarray`
    df = df.rename(columns={"variable": "vintage_year"})

    # Check if we have a fleet composition for the fleet year requested
    if len(df.loc[df["year"] == fleet_year]) == 0:
        # Find the nearest smallest year instead
        list_years = df["year"].unique()
        index = bisect.bisect(list_years, fleet_year)
        substitute_year = list_years[index-1]
        print("Fleet information for {} is not available. We'll use {} instead.".format(fleet_year, substitute_year))
        df = df.loc[df["year"] == substitute_year]
    else:
        df = df.loc[df["year"] == fleet_year]

    df = df.loc[df["REMIND_region"] == remind_region]

    # Associate diesel shares
    def get_diesel_factor(row):

        if row["iso_2"] in shares.columns:
            return float(shares.loc[row["vintage_year"], row["iso_2"]])
        else:
            return 0


    df["diesel_factor"] = df.apply(get_diesel_factor, axis=1)

    new_df_p = df.loc[df["powertrain"].isin(("HEV", "PHEV", "ICEV"))]
    new_df_p.loc[:, "powertrain"] += "-p"
    new_df_p.loc[:, "vintage_demand_vkm"] *= 1 - new_df_p.loc[:, "diesel_factor"]
    new_df_d = df.loc[df["powertrain"].isin(("HEV", "PHEV", "ICEV"))]
    new_df_d.loc[:, "powertrain"] += "-d"
    new_df_d.loc[:, "vintage_demand_vkm"] *= new_df_p.loc[:, "diesel_factor"]

    df = pd.concat([df, new_df_p, new_df_d])
    df = df.loc[~df["powertrain"].isin(("HEV", "PHEV", "ICEV"))]

    # Filter out unecessary columns
    df = df[["powertrain", "size", "year", "vintage_year", "vintage_demand_vkm",]]

    if len(df) == 0:
        raise ValueError("This fleet year is not available.")

    # Distribute the transport demand of 2010 to anterior years

    distr_km = {
        # EURO 4 - 60%
        2010: 0.15,
        2009: 0.15,
        2008: 0.10,
        2007: 0.10,
        2006: 0.10,
        # EURO 3 - 25%
        2005: 0.05,
        2004: 0.05,
        2003: 0.05,
        2002: 0.05,
        2001: 0.05,
        # EURO 2 - 10%
        2000: 0.025,
        1999: 0.025,
        1998: 0.025,
        1997: 0.025,
        # EURO 1 - 5%
        1996: 0.05,
    }

    new_df = pd.DataFrame()
    for y in range(2010, 1995, -1):
        temp_df = df.loc[df["vintage_year"] == 2010]
        temp_df["vintage_year"] = y
        temp_df["vintage_demand_vkm"] *= distr_km[y]
        new_df = pd.concat([new_df, temp_df])

    df = df.loc[df["vintage_year"] != 2010]

    df = pd.concat([df, new_df])

    # Turn the dataframe into a pivot table
    df = df.pivot_table(
        index=["powertrain", "size", "vintage_year"], columns=["year"], aggfunc=np.sum
    )["vintage_demand_vkm"]

    # The table needs to be square: so we add powertrain-size-vintage_year combinations and set their value to 0
    new_cols = [
        c
        for c in df.index.get_level_values("vintage_year").unique()
        if c not in df.columns
    ]
    df[new_cols] = pd.DataFrame([[0] * len(new_cols)], index=df.index)
    pt = df.index.get_level_values("powertrain").unique().tolist()
    y = df.index.get_level_values("vintage_year").unique().tolist()
    y.append(y[-1] + 1)
    s = df.index.get_level_values("size").unique().tolist()
    a = [pt] + [s] + [y]

    for row in [i for i in list(itertools.product(*a)) if i not in df.index]:
        df.loc[row] = 0

    if normalized:
        # The vkm values are normalized to 1, year-wise
        df /= df.sum(axis=0)

    # xarray.DataArray is returned
    return df.to_xarray().fillna(0).to_array()


def build_fleet_array(fp, scope):
    """
    Receives a file path that points to a CSV file that contains the fleet composition
    Checks that the fleet composition array is valid.

    Specifically:

    * the years specified in the fleet must be present in scope["year"]
    * the powertrains specified in the fleet must be present in scope["powertrain"]
    * the sizes specified in the fleet must be present in scope["size"]
    * the sum for each year-powertrain-size set must equal 1

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
