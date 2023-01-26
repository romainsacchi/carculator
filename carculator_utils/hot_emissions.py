"""
hot_emissions.py contains HotEmissionModel
which calculates fuel-related exhaust emissions.
"""

from pathlib import Path
from typing import Any, List, Union

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from xarray import DataArray

from . import DATA_DIR

FILEPATH_DC_SPECS = DATA_DIR / "driving cycle" / "dc_specs.yaml"

MAP_PWT = {
    "Human": "BEV",
    "BEV": "BEV",
    "BEV-depot": "BEV",
    "BEV-opp": "BEV",
    "BEV-motion": "BEV",
    "PHEV-e": "BEV",
    "PHEV-c-d": "ICEV-d",
    "PHEV-d": "ICEV-d",
    "PHEV-c-p": "ICEV-p",
    "PHEV-p": "ICEV-p",
    "ICEV-d": "ICEV-d",
    "ICEV-g": "ICEV-g",
    "ICEV-p": "ICEV-p",
    "HEV-p": "ICEV-p",
    "HEV-d": "ICEV-d",
    "FCEV": "BEV",
}


def _(obj: Union[np.ndarray, xr.DataArray]) -> Union[np.ndarray, xr.DataArray]:
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(obj, (np.ndarray, xr.DataArray)):
        return np.expand_dims(obj, -1)
    return obj


def get_emission_factors(filepath) -> [Any, None]:
    """Hot emissions factors extracted for passenger cars from HBEFA 4.1
    detailed by size, powertrain and EURO class for each substance.
    """

    try:
        df = pd.read_csv(filepath, sep=",")
        cols = ["powertrain", "component"]

        if "euro_class" in df.columns:
            cols.append("euro_class")

        if "size" in df.columns:
            cols.insert(1, "size")
        if "type" in df.columns:
            cols.insert(-2, "type")

        ef = df.groupby(cols)["value"].mean().to_xarray()

    except FileNotFoundError:
        return None

    return ef.fillna(0.0)


def get_mileage_degradation_factor(
    lifetime_km: xr.DataArray,
    euro_class: List[int],
    powertrains: List[str],
    vehicle_type: str,
) -> DataArray | None:
    """
    Catalyst degrade overtime, leading to increased emissions
    of CO, HC and NOX. We apply a correction factor from HBEFA 4.1
    to reflect this.
    :return:
    """

    corr = get_emission_factors(
        filepath=DATA_DIR / "emission_factors" / vehicle_type / f"degradation_EF.csv",
    )

    if corr is None:
        return None

    corr = corr.fillna(1.0)
    corr = corr.expand_dims({"km": 2}).copy()

    if vehicle_type == "car":
        max_km = 200000
    elif vehicle_type in ["bus", "truck"]:
        max_km = 890000
    else:
        max_km = 200000

    corr = corr.assign_coords({"km": np.array([0, max_km])})

    corr = corr.sel(
        powertrain=[
            p for p in lifetime_km.powertrain.values if p in corr.powertrain.values
        ],
        euro_class=euro_class,
    )

    corr = corr.expand_dims({"size": len(lifetime_km.coords["size"])})
    corr = corr.assign_coords({"size": lifetime_km.coords["size"].values})

    corr = corr.transpose("euro_class", "powertrain", "size", "component", "km")

    corr = corr.interp(
        km=lifetime_km.max(), method="linear", kwargs={"fill_value": "extrapolate"}
    )

    corr.values = np.where(corr < 1, 1.0, corr)

    return corr


def get_driving_cycle_compartments(cycle_name, vehicle_type) -> dict:

    with open(FILEPATH_DC_SPECS, "r") as f:
        return yaml.safe_load(f)["environments"][vehicle_type][cycle_name]


class HotEmissionsModel:
    """
    Calculate hot pollutants emissions based on HBEFA 4.1 data, function of fuel consumption
    for vehicles with a combustion engine.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(
        self,
        powertrains: np.array,
        sizes: np.array,
        velocity: np.ndarray,
        cycle_name: str,
        vehicle_type: str,
    ) -> None:

        self.powertrains = powertrains
        self.sizes = sizes
        self.velocity = velocity / 1000 * 3600  # m/s to km/h
        self.cycle_name = cycle_name
        self.vehicle_type = vehicle_type
        self.exhaust = get_emission_factors(
            filepath=DATA_DIR
            / "emission_factors"
            / vehicle_type
            / f"EF_HBEFA42_exhaust.csv"
        )
        self.non_exhaust = get_emission_factors(
            filepath=DATA_DIR
            / "emission_factors"
            / vehicle_type
            / f"EF_HBEFA42_non_exhaust.csv"
        )
        self.nmhc_species = get_emission_factors(
            filepath=DATA_DIR / "emission_factors" / vehicle_type / f"NMHC_species.csv"
        )
        self.engine_wear = get_emission_factors(
            filepath=DATA_DIR / "emission_factors" / vehicle_type / f"engine_wear.csv"
        )

    def get_hot_emissions(
        self,
        euro_class: List[int],
        lifetime_km: xr.DataArray,
        energy_consumption: xr.DataArray,
        yearly_km: xr.DataArray,
    ) -> np.ndarray:
        """
        Calculate hot pollutants emissions given a powertrain type (i.e., diesel, petrol, CNG)
        and a EURO pollution class, per air sub-compartment (i.e., urban, suburban and rural).
        Note that Nh3 and N2O emissions do not depend on the speed level. For those, average
        values observed across different traffic situations are used instead.
        Also includes cold start emissions. Cold start emissions are also given by HBEFA 4.1
        and expressed in given in g/start. Cold start emissions are divided by the average trip
        length in Europe (20 km), to normalize them per km.

        The emission sums are further divided into `air compartments`: urban, suburban and rural.

            * *urban*: from 0 to 50 km/k
            * *suburban*: from 51 km/h to 80 km/h
            * *rural*: above 80 km/h

        :param powertrain_type: "diesel", "petrol" or "CNG"
        :param euro_class: integer, corresponding to the EURO pollution class
        :param energy_consumption: tank-to-wheel energy consumption for each second of the driving cycle
        :param yearly_km: annual mileage, to calculate cold start emissions
        :return: Pollutants emission per km driven, per air compartment.
        """

        hot_emissions = self.exhaust.sel(
            powertrain=[
                MAP_PWT[pt] if MAP_PWT[pt] in self.exhaust.powertrain.values else "BEV"
                for pt in lifetime_km.coords["powertrain"].values
            ],
            euro_class=euro_class,
        )

        if "size" not in hot_emissions.dims:
            hot_emissions = hot_emissions.expand_dims(
                {"size": len(lifetime_km.coords["size"])}
            )
            hot_emissions.coords["size"] = lifetime_km.coords["size"]

        hot_emissions = hot_emissions.transpose(
            "euro_class", "powertrain", "size", "component"
        )

        hot_emissions = hot_emissions.sel(size=lifetime_km.coords["size"].values)

        energy_consumption = energy_consumption.sel(
            size=lifetime_km.coords["size"].values,
            powertrain=lifetime_km.coords["powertrain"].values,
            year=lifetime_km.coords["year"].values,
        )

        distance = self.velocity.sum(dim="second") / 1000

        # Emissions for each second of the driving cycle equal:
        # a * energy consumption
        # with a being a coefficient given by fitting HBEFA 4.1 data
        # the fitting of emissions function of energy consumption
        # is described in the notebook
        # `HBEFA trucks.ipynb` in the folder `dev`.

        # energy consumption is given in kj for each second
        # emissions are in grams per MJ

        emissions = xr.DataArray(
            hot_emissions.values * _(energy_consumption / 1000),
            dims=["second", "value", "year", "powertrain", "size", "component"],
            coords={
                "second": energy_consumption.coords["second"],
                "value": energy_consumption.coords["value"],
                "year": energy_consumption.coords["year"],
                "powertrain": energy_consumption.coords["powertrain"],
                "size": energy_consumption.coords["size"],
                "component": hot_emissions.coords["component"],
            },
        )

        # bit of a manual calibration for N2O and NH3
        # as they do not correlate with fuel consumption

        if self.vehicle_type == "car":
            emissions.loc[dict(component="Dinitrogen oxide")] *= 0.5
            emissions.loc[dict(component="Ammonia")] *= 0.5
        elif self.vehicle_type == "truck":
            emissions.loc[dict(component="Dinitrogen oxide")] *= 10
            emissions.loc[dict(component="Ammonia")] *= 10
        elif self.vehicle_type == "bus":
            emissions.loc[dict(component="Dinitrogen oxide")] /= 15
            emissions.loc[dict(component="Ammonia")] /= 12
        else:
            pass

        # apply a mileage degradation factor for CO, HC and NOx
        degradation_correction = get_mileage_degradation_factor(
            lifetime_km=lifetime_km,
            euro_class=euro_class,
            powertrains=emissions.powertrain.values,
            vehicle_type=self.vehicle_type,
        )

        if degradation_correction is not None:
            emissions.loc[
                dict(
                    powertrain=degradation_correction.powertrain.values,
                    component=degradation_correction.component.values,
                )
            ] *= degradation_correction.values

        if self.non_exhaust is not None:
            non_exhaust = self.non_exhaust.sel(
                powertrain=[
                    MAP_PWT[pt]
                    if MAP_PWT[pt] in self.non_exhaust.powertrain.values
                    else "BEV"
                    for pt in emissions.powertrain.values
                ],
                euro_class=euro_class,
                component=[
                    c
                    for c in emissions.component.values
                    if c in self.non_exhaust.component.values
                ],
            )

            if "size" not in non_exhaust.dims:
                non_exhaust = non_exhaust.expand_dims(
                    {"size": len(energy_consumption.coords["size"])}
                )
                non_exhaust.coords["size"] = energy_consumption.coords["size"]

            non_exhaust = non_exhaust.transpose(
                "euro_class", "powertrain", "size", "component", "type"
            )

            non_exhaust = non_exhaust.sel(size=lifetime_km.coords["size"].values)

            start_per_day = 2.3  # source for

            # Cold start and soak emissions are defined per start and stop
            # Therefore, we need to normalize per km
            # And add cold start emissions to the first second of the driving cycle

            yearly_km = yearly_km.transpose("value", "year", "powertrain", "size")

            emissions.loc[
                dict(
                    second=emissions.second.values[0],
                    component=non_exhaust.component.values,
                )
            ] += (
                _(distance / yearly_km * start_per_day * 365)
                * non_exhaust.sel(type="cold start").values
            )

            # And add soak emissions to the last second of the driving cycle
            emissions.loc[
                dict(
                    second=emissions.second.values[-1],
                    component=non_exhaust.component.values,
                )
            ] += (
                _(distance / yearly_km * start_per_day * 365)
                * non_exhaust.loc[dict(type="soak")].values
            )

            # Diurnal emissions are defined in g/day
            # And need to be evenly distributed
            # throughout the driving cycle

            daily_km_to_year = distance / (yearly_km / 365)

            emissions.loc[dict(component=non_exhaust.component.values)] += (
                _(daily_km_to_year)
                * non_exhaust.sel(type="diurnal").values
                / emissions.shape[0]
            )

            # Running losses are in g/km (no conversion needed)
            # And need to be evenly distributed throughout the driving cycle

            emissions.loc[dict(component=non_exhaust.component.values)] += (
                _(distance)
                * non_exhaust.sel(type="running losses").values
                / emissions.shape[0]
            )

        # Add additional species derived from NMHC emissions
        # Toluene, Xylene, Formaldehyde, Acetaldehyde, etc.
        # Also heavy metals
        # And add them to the emissions array

        nmhc = self.nmhc_species.sel(
            powertrain=[
                MAP_PWT[pt] if pt in self.nmhc_species.powertrain.values else "BEV"
                for pt in emissions.powertrain.values
            ]
        )
        nmhc = nmhc.assign_coords({"powertrain": emissions.powertrain.values})

        if "size" not in nmhc.dims:
            nmhc = nmhc.expand_dims({"size": len(energy_consumption.coords["size"])})
            nmhc.coords["size"] = energy_consumption.coords["size"]

        if "year" not in nmhc.dims:
            nmhc = nmhc.expand_dims({"year": len(energy_consumption.coords["year"])})
            nmhc.coords["year"] = energy_consumption.coords["year"]

        nmhc = nmhc.transpose(
            "year",
            "powertrain",
            "size",
            "component",
        )
        nmhc = nmhc.sel(size=energy_consumption.coords["size"].values)

        emissions.loc[dict(component=nmhc.coords["component"].values)] = (
            nmhc.values * emissions.sel(component=["Non-methane hydrocarbon"]).values
        )

        emissions.loc[dict(component="Non-methane hydrocarbon")] *= nmhc.sum(
            dim="component"
        ).values

        # Heavy metals emissions are dependent of fuel consumption
        # given in grams of emission per kj

        engine_wear = self.engine_wear.sel(
            powertrain=[
                MAP_PWT[pt] if pt in self.engine_wear.powertrain.values else "BEV"
                for pt in emissions.powertrain.values
            ]
        )
        engine_wear = engine_wear.assign_coords(
            {"powertrain": emissions.powertrain.values}
        )

        if "size" not in engine_wear.dims:
            engine_wear = engine_wear.expand_dims(
                {"size": len(energy_consumption.coords["size"])}
            )
            engine_wear.coords["size"] = energy_consumption.coords["size"]

        if "year" not in engine_wear.dims:
            engine_wear = engine_wear.expand_dims(
                {"year": len(energy_consumption.coords["year"])}
            )
            engine_wear.coords["year"] = energy_consumption.coords["year"]

        engine_wear = engine_wear.transpose(
            "year",
            "powertrain",
            "size",
            "component",
        )

        emissions.loc[dict(component=engine_wear.coords["component"].values)] = (
            energy_consumption * engine_wear
        ).values

        # reorder the component coordinates by alphabetical order
        emissions = emissions.reindex(component=sorted(emissions.component.values))

        # urban emissions are the sum of emissions
        # along the ``second`` dimension
        # where velocity is below 50 km/h
        # converted to kg/km
        urban_emissions = (
            np.where(_(self.velocity) <= 50, emissions, 0).sum(axis=0)
            / _(distance)
            / 1000
        )

        # rural emissions are the sum of emissions
        # along the ``second`` dimension
        # where velocity is between 50 km/h
        # and 80 km/h
        rural_emissions = (
            np.where(
                (_(self.velocity) > 50) & (_(self.velocity) <= 80), emissions, 0
            ).sum(axis=0)
            / _(distance)
            / 1000
        )

        # highway emissions are the sum of emissions
        # along the ``second`` dimension
        # where velocity is above 80 km/h
        highway_emissions = (
            np.where(_(self.velocity) > 80, emissions, 0).sum(axis=0)
            / _(distance)
            / 1000
        )

        res = np.concatenate(
            (urban_emissions, rural_emissions, highway_emissions), axis=-1
        )

        return res.transpose(3, 2, 4, 1, 0)
