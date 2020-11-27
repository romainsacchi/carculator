import pandas as pd
import pycountry
from . import DATA_DIR
import itertools
import csv
import numpy as np
import xarray

pd.options.mode.chained_assignment = None

REGION_MAPPING_FILEPATH = DATA_DIR / "regionmappingH12.csv"
REMIND_ELEC_MARKETS = DATA_DIR / "remind_electricity_markets.csv"
IEA_DIESEL_SHARE = DATA_DIR / "diesel_share_oecd.csv"


def extract_biofuel_shares_from_REMIND(fp, remind_region, years, allocate_all_synfuel=False):
    """
    This function extracts biofuel shares from a REMIND file provided.

    :param fp: file path to REMIND file
    :type fp: str
    :param remind_region: REMIND region for which to extract the biofuel shares
    :type remind_region: str
    :param years: the list of years for which to extract biofuel shares
    :param allocate_all_synfuel: Temporary workaround. If True, then all synfuel in the transport sector
    is allocated to passenger cars.
    :return: a dictionary that contains fuel types as keys and lists of fuel shares as values
    """

    df = pd.read_csv(fp, delimiter=";", index_col=["Region", "Variable", "Unit"]).drop(
        columns=["Model", "Scenario"]
    )
    df = df.reset_index()
    df = df.loc[df["Region"] == remind_region]
    df = df.loc[:, : str(2050)]
    df["Variable"] = df["Variable"].str.replace("|", "-")

    if allocate_all_synfuel:

        # get shares of synthetic fuel
        df_total = df.loc[df["Variable"] == "FE-Transport-Pass-Road-LDV-Liquids"]
        df_total.index = df.loc[df["Variable"] == "SE-Liquids-Hydrogen"].index
        share = np.clip(
            (df.loc[df["Variable"] == "SE-Liquids-Hydrogen", "2005":].divide(
                df_total.loc[:, "2005":], axis=0)
            ).values,
            0, 1)

        var = ["FE-Transport-Liquids-Oil",
               "FE-Transport-Liquids-Biomass"]

        df_liquids = df.loc[df["Variable"].isin(var), :]
        df_liquids.iloc[:, 3:] /= df_liquids.iloc[:, 3:].sum(axis=0)

        df_liquids.loc[:, "2005":] *= (1 - share)
        to_append = [remind_region, "FE-Transport-Liquids-Hydrogen", "EJ/yr"] + share[0].tolist()
        a_series = pd.Series(to_append, index=df_liquids.columns)
        df_liquids = df_liquids.append(a_series, ignore_index=True)

    else:

        var = ["FE-Transport-Liquids-Oil",
               "FE-Transport-Liquids-Biomass",
               "FE-Transport-Liquids-Hydrogen"]

        df_liquids = df.loc[df["Variable"].isin(var), :]
        df_liquids.iloc[:, 3:] /= df_liquids.iloc[:, 3:].sum(axis=0)

    var = ["FE-Transport-Gases-Non-Biomass", "FE-Transport-Gases-Biomass"]

    df_gas = df.loc[df["Variable"].isin(var), :]
    df_gas.iloc[:, 3:] /= df_gas.iloc[:, 3:].sum(axis=0)

    d_map_fuels = {
        "FE-Transport-Liquids-Oil": "liquid - fossil",
        "FE-Transport-Liquids-Biomass": "liquid - biomass",
        "FE-Transport-Liquids-Hydrogen": "liquid - synfuel",
        "FE-Transport-Gases-Non-Biomass": "gas - fossil",
        "FE-Transport-Gases-Biomass": "gas - biomass",
    }

    new_df = pd.concat([df_liquids, df_gas])
    new_df["Variable"] = new_df["Variable"].map(d_map_fuels)
    new_df.columns = (
        new_df.columns[:3].tolist() + new_df.columns[3:].astype(int).tolist()
    )

    new_df = new_df.rename(columns={"Variable": "fuel_type"})
    new_df = new_df.groupby("fuel_type").sum()

    return new_df.to_xarray().to_array().interp(variable=years)

def extract_electricity_mix_from_REMIND_file(fp, remind_region, years):
    """
    This function extracts electricity mixes from a REMIND file provided.

    :param fp: file path to REMIND file
    :type fp: str
    :param remind_region: REMIND region for which to extract the electricity mix
    :type remind_region: str
    :param years: the list of years for which to extract electricity mixes
    :type years: list
    :return: list of lists of electricity mixes, that can be consumed by `InventnoryCalculation`
    :rtype: list
    """
    with open(REMIND_ELEC_MARKETS) as f:
        electricity_markets = dict(filter(None, csv.reader(f, delimiter=";")))

    rev_tech = {v: k for k, v in electricity_markets.items()}

    df = pd.read_csv(fp, delimiter=";", index_col=["Region", "Variable", "Unit"]).drop(
        columns=["Model", "Scenario"]
    )

    df = df.reset_index()
    df = df.loc[df["Region"] == remind_region]
    df = df.loc[:, : str(2050)]
    df = df.loc[df["Variable"].isin(electricity_markets.values())]
    df["Variable"] = df["Variable"].map(rev_tech)
    df.columns = df.columns[:3].tolist() + df.columns[3:].astype(int).tolist()
    df.iloc[:, 3:] /= df.iloc[:, 3:].sum(axis=0)

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
        "Oil": "Oil",
        "Nuclear": "Nuclear",
        "Geothermal": "Geothermal",
        "Hydro": "Hydro",
        "Solar CSP": "Solar",
        "Solar PV": "Solar",
        "Wind": "Wind",
    }

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
        df.loc[df.index[-1] + 1] = [remind_region, row, "EJ/yr"] + [0] * (
            df.shape[1] - 3
        )

    df = df.groupby("Variable").sum().loc[list_tech]

    return df.to_xarray().to_array().interp(variable=years).values

def create_fleet_composition_from_REMIND_file(fp, remind_region, fleet_year=2020, normalized=True):
    """
    This function creates a consumable fleet composition array from a CSV file.
    The array returned is consumed by `InventoryCalculation`.

    :param fp: filepath to CSV file
    :type fp: str
    :param remind_region: REMIND region for which to extract fleet composition
    :type remind_region: str
    :param fleet_year: the year for which to extract fleet composition
    :type fleet_year: int
    :return: fleet composition array
    :rtype: xarray.DataArray
    """

    # Load diesel shares from IEA on OECD countries
    share_diesel = pd.read_csv(IEA_DIESEL_SHARE,
                               sep=";", usecols=range(1, 10), index_col=[0])
    share_diesel.columns = [int(c) for c in share_diesel.columns]
    shares = pd.concat([
        share_diesel,
        pd.DataFrame(columns=range(2019, 2051)),
        pd.DataFrame(columns=range(2000, 2011))])
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
    new_df_p.loc[:, "vintage_demand_vkm"] *= (1 - new_df_p.loc[:, "diesel_factor"])
    new_df_d = df.loc[df["powertrain"].isin(("HEV", "PHEV", "ICEV"))]
    new_df_d.loc[:, "powertrain"] += "-d"
    new_df_d.loc[:, "vintage_demand_vkm"] *= new_df_p.loc[:, "diesel_factor"]

    df = pd.concat([df, new_df_p, new_df_d])
    df = df.loc[~df["powertrain"].isin(("HEV", "PHEV", "ICEV"))]



    # Filter out unecessary columns
    df = df[
        [
            "powertrain",
            "size",
            "year",
            "vintage_year",
            "vintage_demand_vkm",
        ]
    ]

    if len(df)==0:
        raise ValueError("This fleet year is not available.")


    # Distribute the transport demand of 2010 to anterior years

    distr_km = {
        2010: .25,
        2009: .15,
        2008: .15,
        2007: .15,
        2006: .15,
        2005: .15
    }

    new_df = pd.DataFrame()
    for y in range(2010, 2004, -1):
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
    y.append(y[-1]+1)
    s = df.index.get_level_values("size").unique().tolist()
    a = [pt] + [s] + [y]

    for row in [i for i in list(itertools.product(*a)) if i not in df.index]:
        df.loc[row] = 0

    if normalized:
        # The vkm values are normalized to 1, year-wise
        df /= df.sum(axis=0)

    # xarray.DataArray is returned
    return df.to_xarray().fillna(0).to_array()
