import pickle

import numpy as np
import xarray as xr

from . import DATA_DIR


def _(o):
    """Add a trailing dimension to make input arrays broadcast correctly"""
    if isinstance(o, (np.ndarray, xr.DataArray)):
        return np.expand_dims(o, -1)
    else:
        return o


def get_hot_emission_factors():
    """Hot emissions factors extracted for trucks from HBEFA 4.1
    detailed by size, powertrain and EURO class for each substance.
    """
    fp = DATA_DIR / "hot.pickle"

    with open(fp, "rb") as f:
        hot = pickle.load(f)

    return hot


def get_non_hot_emission_factors():
    """Non hot emissions factors (cold start, evaporation, soak emissions) extracted for trucks from HBEFA 4.1
    detailed by size, powertrain and EURO class for each substance.
    """
    fp = DATA_DIR / "non_hot.pickle"

    with open(fp, "rb") as f:
        non_hot = pickle.load(f)

    return non_hot


def get_mileage_degradation_factor(powertrain_type, euro_class, lifetime_km):
    """
    Catalyst degrade overtime, leading to increased emissions
    of CO, HC and NOX. We apply a correction factor
    to reflect this.
    :return:
    """

    d_corr = {
        "ICEV-p": {
            "CO": {
                1: 1.9,
                2: 1.6,
                3: 1.75,
                4: 1.9,
                5: 2,
                6: 1.3,
                6.1: 1.3,
                6.2: 1.3,
                6.3: 1.3,
            },
            "HC": {1: 1.59, 2: 1.59, 3: 1.02, 4: 1.02},
            "NOx": {
                1: 2.5,
                2: 2.3,
                3: 2.9,
                4: 2,
                5: 2.5,
                6: 1.3,
                6.1: 1.3,
                6.2: 1.3,
                6.3: 1.3,
            },
        },
        "ICEV-d": {
            "CO": {
                4: 1.3,
                5: 1.3,
                6: 1.4,
                6.1: 1.4,
                6.2: 1.4,
                6.3: 1.4,
            },
            "NOx": {
                2: 1.25,
                3: 1.2,
                4: 1.06,
                5: 1.03,
                6: 1.15,
                6.1: 1.15,
                6.2: 1.15,
                6.3: 1.15,
            },
        },
    }

    shape = list(lifetime_km.shape)
    shape[-1] = 3
    corr = np.ones(tuple(shape))

    for p, pt in enumerate(powertrain_type):
        for e, ec in enumerate(euro_class):
            for c, co in enumerate(["CO", "HC", "NOx"]):
                try:
                    val = d_corr[pt][co][ec]
                    y_max = 120000 if co == "HC" else 200000
                    corr[:, p, e, c] = np.clip(
                        np.interp(lifetime_km[:, p, e, 0] / 2, [0, y_max], [1, val]),
                        1,
                        None,
                    )

                except KeyError:
                    pass

    return corr.transpose(3, 0, 1, 2)


class HotEmissionsModel:
    """
    Calculate hot pollutants emissions based on HBEFA 4.1 data, function of fuel consumption
    for vehicles with a combustion engine.

    :param cycle: Driving cycle. Pandas Series of second-by-second speeds (km/h) or name (str)
        of cycle e.g., "WLTC","WLTC 3.1","WLTC 3.2","WLTC 3.3","WLTC 3.4","CADC Urban","CADC Road",
        "CADC Motorway","CADC Motorway 130","CADC","NEDC".
    :type cycle: pandas.Series

    """

    def __init__(self, cycle, cycle_name):

        self.cycle = cycle
        self.cycle_name = cycle_name

        self.cycle_environment = {
            "WLTC": {
                "urban start": 0,
                "urban stop": 590,
                "suburban start": 591,
                "suburban stop": 1023,
                "rural start": 1024,
                "rural stop": 1801,
            },
            "WLTC 3.1": {"urban start": 0, "urban stop": 590},
            "WLTC 3.2": {"suburban start": 0, "suburban stop": 433},
            "WLTC 3.3": {"rural start": 0, "rural stop": 455},
            "WLTC 3.4": {"rural start": 0, "rural stop": 323},
            "CADC Urban": {"urban start": 0, "urban stop": 994},
            "CADC Road": {"suburban start": 0, "suburban stop": 1082},
            "CADC Motorway": {"rural start": 0, "rural stop": 1068},
            "CADC Motorway 130": {"rural start": 0, "rural stop": 1068},
            "CADC": {
                "urban start": 0,
                "urban stop": 994,
                "suburban start": 995,
                "suburban stop": 2077,
                "rural start": 2078,
                "rural stop": 3144,
            },
            "NEDC": {
                "urban start": 0,
                "urban stop": 780,
                "rural start": 781,
                "rural stop": 1180,
            },
        }

        self.hot = get_hot_emission_factors()
        self.non_hot = get_non_hot_emission_factors()

    def get_hot_emissions(
        self, powertrain_type, euro_class, lifetime_km, energy_consumption, yearly_km
    ):
        """
        Calculate hot pollutants emissions given a powertrain type (i.e., diesel, petrol, CNG) and a EURO pollution class, per air sub-compartment
        (i.e., urban, suburban and rural).
        Note that Nh3 and N2O emissions do not depend on the speed level. FOr those, average values observed across
        different traffic situations are used instead.
        Also includes cold start emissions. Cold start emissions are also given by HBEFA 4.1 and expressed in given in g/start.
        Cold start emissions are divided by the average trip length in Europe (20 km), to normalize them per km.

        The emission sums are further divided into `air compartments`: urban, suburban and rural.

            * *urban*: from 0 to 50 km/k
            * *suburban*: from 51 km/h to 80 km/h
            * *rural*: above 80 km/h

        :param powertrain_type: "diesel", "petrol" or "CNG"
        :type powertrain_type: str
        :param euro_class: integer, corresponding to the EURO pollution class
        :type euro_class: float
        :param energy_consumption: tank-to-wheel energy consumption for each second of the driving cycle
        :type energy_consumption: xarray
        :param yearly_km: annual mileage, to calculate cold start emissions
        :return: Pollutants emission per km driven, per air compartment.
        :rtype: numpy.array
        """

        # Check if the powertrains passed are valid
        if set(powertrain_type).intersection({"FCEV", "BEV", "PHEV-e"}):
            raise TypeError("Wrong powertrain!")

        hot_emissions = self.hot.sel(
            powertrain=powertrain_type,
            euro_class=euro_class,
            component=[
                "HC",
                "CO",
                "NOx",
                "PM2.5",
                "CH4",
                "NMHC",
                "N2O",
                "NH3",
                "Pb",
                "Benzene",
            ],
        ).transpose("component", "powertrain", "euro_class", "variable")

        distance = self.cycle.sum() / 3600

        # Emissions for each second of the driving cycle equal:
        # a * energy consumption
        # with a being a coefficient given by fitting HBEFA 4.1 data
        # the fitting of emissions function of energy consumption is described in the notebook
        # `HBEFA trucks.ipynb` in the folder `dev`.

        # energy conusmption is given in kj for each second
        # emissions are in grams per MJ
        hot = (
            hot_emissions.sel(variable="a").values[:, None, :, :, None, None]
            * energy_consumption.values
        )
        # bit of a manual calibration for N2O and NH3
        # as they do not correlate with fuel consumption
        hot[6] *= 0.5
        hot[7] *= 0.5

        # apply a mileage degradation factor for CO, HC and NOx
        corr = get_mileage_degradation_factor(powertrain_type, euro_class, lifetime_km)
        hot[:3] *= corr[:, :, :, :, None, None]

        non_hot_emissions = self.non_hot.sel(
            powertrain=powertrain_type,
            euro_class=euro_class,
            Component=[
                "HC",
                "CO",
                "NOx",
                "PM2.5",
                "CH4",
                "NMHC",
                "N2O",
                "NH3",
                "Pb",
                "Benzene",
            ],
        ).transpose("Component", "powertrain", "euro_class", "type")

        start_per_day = 2.3  # source for

        # Cold start and soak emissions are defined per start and stop
        # Therefore, we need to normalize per km
        # And add cold start emissions to the first second of the driving cycle
        hot[..., 0] += (start_per_day * 365 / yearly_km).values * non_hot_emissions.loc[
            dict(type="cold start")
        ].values[:, None, ..., None]
        # And add soak emissions to the last second of the driving cycle
        hot[..., -1] += (
            start_per_day * 365 / yearly_km
        ).values * non_hot_emissions.loc[dict(type="soak")].values[:, None, ..., None]

        # Diurnal emissions are defined in g/day
        # And need to be evenly distributed throughout the driving cycle
        hot += (
            (365 / yearly_km).values
            * non_hot_emissions.loc[dict(type="diurnal")].values[:, None, ..., None]
        )[..., None] / len(self.cycle)

        # Running losses are in g/km (no conversion needed)
        # And need to be evenly distributed throughout the driving cycle
        hot += non_hot_emissions.loc[dict(type="running losses")].values[
            :, None, ..., None, None
        ] / len(self.cycle)

        # Add additional species derived from NMHC emissions
        # Toluene, Xylene, Formaldehyde, Acetaldehyde, etc.
        # Also heavy metals
        shape = list(hot.shape)
        shape[0] = 40

        arr = np.zeros(shape)

        arr[
            : hot.shape[0],
            : hot.shape[1],
            : hot.shape[2],
            : hot.shape[3],
            : hot.shape[4],
            : hot.shape[5],
        ] = hot

        final_emissions = xr.DataArray(
            arr,
            coords=[
                [
                    "HC",
                    "CO",
                    "NOx",
                    "PM2.5",
                    "CH4",
                    "NMHC",
                    "N2O",
                    "NH3",
                    "Pb",
                    "Benzene",
                    "Ethane",
                    "Propane",
                    "Butane",
                    "Pentane",
                    "Hexane",
                    "Cyclohexane",
                    "Heptane",
                    "Ethene",
                    "Propene",
                    "1-Pentene",
                    "Toluene",
                    "m-Xylene",
                    "o-Xylene",
                    "Formaldehyde",
                    "Acetaldehyde",
                    "Benzaldehyde",
                    "Acetone",
                    "Methyl ethyl ketone",
                    "Acrolein",
                    "Styrene",
                    "PAH, polycyclic aromatic hydrocarbons",
                    "Arsenic",
                    "Selenium",
                    "Zinc",
                    "Copper",
                    "Nickel",
                    "Chromium",
                    "Chromium VI",
                    "Mercury",
                    "Cadmium",
                ],
                energy_consumption.coords["size"].values,
                energy_consumption.powertrain.values,
                energy_consumption.year.values,
                range(0, arr.shape[4]),
                range(0, arr.shape[-1]),
            ],
            dims=["component", "size", "powertrain", "year", "value", "second"],
        )

        NMHC_ratios_petrol = np.array(
            [
                0.032,
                0.007,
                0.052,
                0.022,
                0.016,
                0.011,
                0.007,
                0.073,
                0.038,
                0.001,
                0.110,
                0.054,
                0.023,
                0.017,
                0.008,
                0.002,
                0.006,
                0.001,
                0.002,
                0.010,
            ]
        )

        NMHC_ratios_diesel = np.array(
            [
                0.0033,
                0.0011,
                0.0011,
                0.0004,
                0.0000,
                0.0065,
                0.0020,
                0.1097,
                0.0360,
                0.0000,
                0.0069,
                0.0061,
                0.0027,
                0.1200,
                0.0647,
                0.0086,
                0.0294,
                0.0120,
                0.0358,
                0.0037,
            ]
        )

        nmvoc_sub_species = [
            "Ethane",
            "Propane",
            "Butane",
            "Pentane",
            "Hexane",
            "Cyclohexane",
            "Heptane",
            "Ethene",
            "Propene",
            "1-Pentene",
            "Toluene",
            "m-Xylene",
            "o-Xylene",
            "Formaldehyde",
            "Acetaldehyde",
            "Benzaldehyde",
            "Acetone",
            "Methyl ethyl ketone",
            "Acrolein",
            "Styrene",
        ]

        final_emissions.loc[
            dict(
                component=nmvoc_sub_species,
                powertrain=[p for p in final_emissions.powertrain.values if "d" in p],
            )
        ] = (
            NMHC_ratios_diesel.reshape((-1, 1, 1, 1, 1, 1))
            * final_emissions.loc[
                dict(
                    component="NMHC",
                    powertrain=[
                        p for p in final_emissions.powertrain.values if "-d" in p
                    ],
                )
            ].values
        )

        final_emissions.loc[
            dict(
                component="NMHC",
                powertrain=[p for p in final_emissions.powertrain.values if "-d" in p],
            )
        ] *= (
            1 - 0.45
        )

        final_emissions.loc[
            dict(
                component=nmvoc_sub_species,
                powertrain=[p for p in final_emissions.powertrain.values if "-p" in p],
            )
        ] = (
            NMHC_ratios_petrol.reshape((-1, 1, 1, 1, 1, 1))
            * final_emissions.loc[
                dict(
                    component="NMHC",
                    powertrain=[
                        p for p in final_emissions.powertrain.values if "-p" in p
                    ],
                )
            ].values
        )

        final_emissions.loc[
            dict(
                component="NMHC",
                powertrain=[p for p in final_emissions.powertrain.values if "-p" in p],
            )
        ] *= (
            1 - 0.492
        )

        # Heavy metals emissions are dependent of fuel consumption
        # given in grams of emission per kj
        heavy_metals = [
            "PAH, polycyclic aromatic hydrocarbons",
            "Arsenic",
            "Selenium",
            "Zinc",
            "Copper",
            "Nickel",
            "Chromium",
            "Chromium VI",
            "Mercury",
            "Cadmium",
        ]
        final_emissions.loc[
            dict(
                component=heavy_metals,
                powertrain=[p for p in final_emissions.powertrain.values if "-p" in p],
            )
        ] = (
            np.array(
                [
                    8.19e-10,
                    7.06e-12,
                    4.71e-12,
                    5.08e-08,
                    9.88e-10,
                    3.06e-10,
                    3.76e-10,
                    7.53e-13,
                    2.05e-10,
                    2.54e-10,
                ]
            ).reshape((-1, 1, 1, 1, 1, 1))
            * energy_consumption.sel(
                powertrain=[
                    p for p in energy_consumption.powertrain.values if "-p" in p
                ]
            ).values
        )

        final_emissions.loc[
            dict(
                component=heavy_metals,
                powertrain=[p for p in final_emissions.powertrain.values if "-d" in p],
            )
        ] = (
            np.array(
                [
                    1.32e-09,
                    2.33e-12,
                    2.33e-12,
                    4.05e-08,
                    4.93e-10,
                    2.05e-10,
                    6.98e-10,
                    1.40e-12,
                    1.23e-10,
                    2.02e-10,
                ]
            ).reshape((-1, 1, 1, 1, 1, 1))
            * energy_consumption.sel(
                powertrain=[
                    p for p in energy_consumption.powertrain.values if "-d" in p
                ]
            ).values
        )

        if self.cycle_name in self.cycle_environment:
            if "urban start" in self.cycle_environment[self.cycle_name]:
                urban_start = self.cycle_environment[self.cycle_name]["urban start"]
                urban_stop = self.cycle_environment[self.cycle_name]["urban stop"]
            else:
                urban_start, urban_stop = [0, 0]

            if "suburban start" in self.cycle_environment[self.cycle_name]:
                suburban_start = self.cycle_environment[self.cycle_name][
                    "suburban start"
                ]
                suburban_stop = self.cycle_environment[self.cycle_name]["suburban stop"]
            else:
                suburban_start, suburban_stop = [0, 0]

            if "rural start" in self.cycle_environment[self.cycle_name]:
                rural_start = self.cycle_environment[self.cycle_name]["rural start"]
                rural_stop = self.cycle_environment[self.cycle_name]["rural stop"]
            else:
                rural_start, rural_stop = [0, 0]

            urban = (
                final_emissions.sel(second=range(urban_start, urban_stop)).sum(axis=-1)
                / distance
            )
            suburban = (
                final_emissions.sel(second=range(suburban_start, suburban_stop)).sum(
                    axis=-1
                )
                / distance
            )
            rural = (
                final_emissions.sel(second=range(rural_start, rural_stop)).sum(axis=-1)
                / distance
            )

        else:
            urban = final_emissions.where(self.cycle <= 50).sum(axis=-1) / distance
            suburban = (
                final_emissions.where((self.cycle > 50) & (self.cycle <= 80)).sum(
                    axis=-1
                )
                / distance
            )
            rural = final_emissions.where(self.cycle > 80).sum(axis=-1) / distance

        urban *= np.isfinite(urban)
        suburban *= np.isfinite(suburban)
        rural *= np.isfinite(rural)

        urban = urban.fillna(0)
        suburban = suburban.fillna(0)
        rural = rural.fillna(0)

        return np.vstack((urban, suburban, rural)).transpose(1, 2, 0, 3, 4) / 1000
