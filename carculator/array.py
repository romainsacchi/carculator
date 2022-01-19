"""
Module array.py exposes fill_xarray_from_input_parameters() which formats
parameters contained in the CarInputParameters object into a multi-dimensional array.
This array host input as well as not yet calculated parameters.
Also, modify_xarray_from_custom_parameters() offers a way to modify some of the default values
used for input parameters.
"""

from typing import Tuple, Union

import numpy as np
import pandas as pd
import stats_arrays as sa
import xarray as xr

from .car_input_parameters import CarInputParameters as c_i_p


def fill_xarray_from_input_parameters(
    cip: c_i_p, sensitivity: bool = False, scope: dict = None
) -> Tuple[Tuple[dict, dict, dict, dict], xr.DataArray]:

    """Create an `xarray` labeled array from the sampled input parameters.

    This function extracts the parameters' names and values contained in the
    `parameters` attribute of the :class:`CarInputParameters` class in
    :mod:`car_input_parameters` and insert them into a
    multi-dimensional numpy-like array from the *xarray* package
    (http://xarray.pydata.org/en/stable/).


    :param cip: Instance of the :class:`CarInputParameters` class in :mod:`car_input_parameters`.
    :param sensitivity: boolean. Whether a sensitivity test is carried out.
    :param scope: a dictionary to narrow down the scope of vehicles to consider
    :returns: `tuple`, `xarray.DataArray`
    - tuple (`size_dict`, `powertrain_dict`, `parameter_dict`, `year_dict`)
    - array

    Dimensions of `array`:

        0. Vehicle size, e.g. "small", "medium". str.
        1. Powertrain, e.g. "ICE-p", "BEV". str.
        2. Year. int.
        3. Samples.

    """

    # Check whether the argument passed is a cip object
    if not isinstance(cip, c_i_p):
        raise TypeError(
            "The argument passed is not an object of the CarInputParameter class"
        )

    if scope is None:
        scope = {"size": cip.sizes, "powertrain": cip.powertrains, "year": cip.years}
    else:
        if "size" not in scope:
            scope["size"] = cip.sizes
        if "powertrain" not in scope:
            scope["powertrain"] = cip.powertrains
        if "year" not in scope:
            scope["year"] = cip.years

    # Make sure to include PHEV-e and PHEV-c-p/d if
    # PHEV-d or PHEV-p are listed

    if "PHEV-p" in scope["powertrain"]:
        for pwt in ["PHEV-e", "PHEV-c-p"]:
            if pwt not in scope["powertrain"]:
                scope["powertrain"].append(pwt)

    if "PHEV-d" in scope["powertrain"]:
        for pwt in ["PHEV-e", "PHEV-c-d"]:
            if pwt not in scope["powertrain"]:
                scope["powertrain"].append(pwt)

    if any(size for size in scope["size"] if size not in cip.sizes):
        raise ValueError("One of the size types is not valid.")

    if any(y for y in scope["year"] if y not in cip.years):
        raise ValueError("One of the years defined is not valid.")

    if any(pwt for pwt in scope["powertrain"] if pwt not in cip.powertrains):
        raise ValueError("One of the powertrain types is not valid.")

    if not sensitivity:
        array = xr.DataArray(
            np.zeros(
                (
                    len(scope["size"]),
                    len(scope["powertrain"]),
                    len(cip.parameters),
                    len(scope["year"]),
                    cip.iterations or 1,
                )
            ).astype("float32"),
            coords=[
                scope["size"],
                scope["powertrain"],
                cip.parameters,
                scope["year"],
                np.arange(cip.iterations or 1),
            ],
            dims=["size", "powertrain", "parameter", "year", "value"],
        )
    else:
        params = ["reference"]
        params.extend(cip.input_parameters)
        array = xr.DataArray(
            np.zeros(
                (
                    len(scope["size"]),
                    len(scope["powertrain"]),
                    len(cip.parameters),
                    len(scope["year"]),
                    len(params),
                )
            ).astype("float32"),
            coords=[
                scope["size"],
                scope["powertrain"],
                cip.parameters,
                scope["year"],
                params,
            ],
            dims=["size", "powertrain", "parameter", "year", "value"],
        )

    size_dict = {k: i for i, k in enumerate(scope["size"])}
    powertrain_dict = {k: i for i, k in enumerate(scope["powertrain"])}
    year_dict = {k: i for i, k in enumerate(scope["year"])}
    parameter_dict = {k: i for i, k in enumerate(cip.parameters)}

    if not sensitivity:

        for param in cip:

            pwt = (
                set(cip.metadata[param]["powertrain"])
                if isinstance(cip.metadata[param]["powertrain"], list)
                else {cip.metadata[param]["powertrain"]}
            )

            size = (
                set(cip.metadata[param]["sizes"])
                if isinstance(cip.metadata[param]["sizes"], list)
                else {cip.metadata[param]["sizes"]}
            )

            year = (
                set(cip.metadata[param]["year"])
                if isinstance(cip.metadata[param]["year"], list)
                else {cip.metadata[param]["year"]}
            )

            if (
                pwt.intersection(scope["powertrain"])
                and size.intersection(scope["size"])
                and year.intersection(scope["year"])
            ):

                array.loc[
                    dict(
                        powertrain=[p for p in pwt if p in scope["powertrain"]],
                        size=[s for s in size if s in scope["size"]],
                        year=[y for y in year if y in scope["year"]],
                        parameter=cip.metadata[param]["name"],
                    )
                ] = cip.values[param]
    else:

        for param in cip.input_parameters:

            names = [n for n in cip.metadata if cip.metadata[n]["name"] == param]

            for name in names:

                pwt = (
                    set(cip.metadata[name]["powertrain"])
                    if isinstance(cip.metadata[name]["powertrain"], list)
                    else {cip.metadata[name]["powertrain"]}
                )

                size = (
                    set(cip.metadata[name]["sizes"])
                    if isinstance(cip.metadata[name]["sizes"], list)
                    else {cip.metadata[name]["sizes"]}
                )

                year = (
                    set(cip.metadata[name]["year"])
                    if isinstance(cip.metadata[name]["year"], list)
                    else {cip.metadata[name]["year"]}
                )

                vals = [
                    cip.values[name] for _ in range(0, len(cip.input_parameters) + 1)
                ]
                vals[cip.input_parameters.index(param) + 1] *= 1.1

                array.loc[
                    dict(
                        powertrain=[p for p in pwt if p in scope["powertrain"]],
                        size=[s for s in size if s in scope["size"]],
                        year=[y for y in year if y in scope["year"]],
                        parameter=cip.metadata[name]["name"],
                    )
                ] = vals

    return (size_dict, powertrain_dict, parameter_dict, year_dict), array


def modify_xarray_from_custom_parameters(
    filepath: Union[dict, str], array: xr.DataArray
) -> None:
    """
    Override default parameters values in `xarray` based on values provided by the user.

    This function allows to override one or several default parameter values by providing either:

        * a file path to an Excel workbook that contains the new values
        * or a dictionary

    The dictionary must be of the following format:

    .. code-block:: python

            {
                (parameter category,
                    powertrain,
                    size,
                    parameter name,
                    uncertainty type): {
                                        (year, 'loc'): value,
                                        (year, 'scale'): value,
                                        (year, 'shape'): value,
                                        (year, 'minimum'): value,
                                        (year, 'maximum'): value
                }

            }

    For example:

    .. code-block:: python

            {
                ('Driving',
                'all',
                'all',
                'lifetime kilometers',
                'none'): {
                    (2018, 'loc'): 150000, (2040, 'loc'): 150000
                    }

            }

    :param filepath: File path of workbook with new values or dictionary.
    :type filepath: str or dict
    :param array: array with vehicle parameters

    """

    if isinstance(filepath, str):
        try:
            param_dictionary = pd.read_excel(
                filepath,
                header=[0, 1],
                index_col=[0, 1, 2, 3, 4],
                sheet_name="Custom_parameters",
            ).to_dict(orient="index")
        except FileNotFoundError as err:
            raise FileNotFoundError("Custom parameters file not found.") from err
    elif isinstance(filepath, dict):
        param_dictionary = filepath
    else:
        raise TypeError("The format passed as parameter is not valid.")

    forbidden_keys = ["Driving cycle", "Background", "Functional unit"]

    for key in param_dictionary:
        if key[0] not in forbidden_keys:
            if not isinstance(key[1], str):
                powertrains = [p.strip() for p in key[1] if p]
                powertrains = [p for p in powertrains if p]
                powertrains = list(powertrains)
            elif key[1] == "all":
                powertrains = array.coords["powertrain"].values
            else:
                if key[1] in array.coords["powertrain"].values:
                    powertrains = [key[1]]
                elif all(
                    p
                    for p in key[1].split(", ")
                    if p in array.coords["powertrain"].values
                ):
                    powertrains = key[1].split(", ")
                else:
                    print(
                        f"{key[1]} is not a recognized powertrain. It will be skipped."
                    )
                    continue

            if not isinstance(key[2], str):
                sizes = [s.strip() for s in key[2] if s]
                sizes = [s for s in sizes if s]
                sizes = list(sizes)
            elif key[2] == "all":
                sizes = array.coords["size"].values
            else:
                if key[2] in array.coords["size"].values:
                    sizes = [key[2]]
                elif all(
                    s for s in key[2].split(", ") if s in array.coords["size"].values
                ):
                    sizes = key[2].split(", ")
                else:
                    print(
                        f"{key[2]} is not a recognized size category. It will be skipped."
                    )
                    continue

            param = key[3]

            if not param in array.coords["parameter"].values:
                print(f"{param} is not a recognized parameter. It will be skipped.")
                continue

            val = param_dictionary[key]

            distr_dic = {
                "triangular": 5,
                "lognormal": 2,
                "normal": 3,
                "uniform": 4,
                "none": 1,
            }
            distr = distr_dic[key[4]]

            years = {v[0] for v in val}

            for year in years:
                # No uncertainty parameters given
                if distr == 1:
                    # There should be at least a `loc`
                    if ~np.isnan(val[(year, "loc")]):
                        for size in sizes:
                            for pwt in powertrains:
                                array.loc[
                                    dict(
                                        powertrain=pwt,
                                        size=size,
                                        year=year,
                                        parameter=param,
                                    )
                                ] = val[(year, "loc")]
                    # Otherwise warn
                    else:
                        print(f"`loc`parameter missing for {param} in {year}.")
                        continue

                elif distr in [2, 3, 4, 5]:

                    # Check if the correct parameters are present
                    # Triangular

                    if distr == 5:
                        if (
                            np.isnan(val[(year, "loc")])
                            or np.isnan(val[(year, "minimum")])
                            or np.isnan(val[(year, "maximum")])
                        ):
                            print(
                                "One or more parameters for the triangular distribution "
                                f"is/are missing for {param} in {year}.\n "
                                "The parameter is skipped and default value applies."
                            )
                            continue

                    # Lognormal
                    if distr == 2:
                        if np.isnan(val[(year, "loc")]) or np.isnan(
                            val[(year, "scale")]
                        ):
                            print(
                                "One or more parameters for the lognormal distribution "
                                f"is/are missing for {param} in {year}.\n "
                                "The parameter is skipped and default value applies"
                            )
                            continue

                    # Normal
                    if distr == 3:
                        if np.isnan(val[(year, "loc")]) or np.isnan(
                            val[(year, "scale")]
                        ):
                            print(
                                "One or more parameters for the normal distribution is/are missing"
                                f" for {param} in {year}.\n"
                                f"The parameter is skipped and default value applies"
                            )
                            continue

                    # Uniform
                    if distr == 4:
                        if np.isnan(val[(year, "minimum")]) or np.isnan(
                            val[(year, "maximum")]
                        ):
                            print(
                                "One or more parameters for the uniform distribution "
                                f"is/are missing for {param} in {year}.\n "
                                "The parameter is skipped and default value applies"
                            )
                            continue

                    uncertainty_params = sa.UncertaintyBase.from_dicts(
                        {
                            "loc": val[year, "loc"],
                            "scale": val[year, "scale"],
                            "shape": val[year, "shape"],
                            "minimum": val[year, "minimum"],
                            "maximum": val[year, "maximum"],
                            "uncertainty_type": distr,
                        }
                    )

                    # Stochastic mode
                    if array.sizes["value"] > 1:

                        rng = sa.MCRandomNumberGenerator(uncertainty_params)

                        for size in sizes:
                            for pwt in powertrains:
                                array.loc[
                                    dict(
                                        powertrain=pwt,
                                        size=size,
                                        year=year,
                                        parameter=param,
                                    )
                                ] = rng.generate(array.sizes["value"]).reshape((-1,))
                    else:

                        dist = sa.uncertainty_choices[distr]
                        median = float(dist.ppf(uncertainty_params, np.array((0.5,))))

                        for size in sizes:
                            for pwt in powertrains:
                                array.loc[
                                    dict(
                                        powertrain=pwt,
                                        size=size,
                                        year=year,
                                        parameter=param,
                                    )
                                ] = median

                else:
                    print(
                        f"The uncertainty type is not recognized for {param} in {year}.\n"
                        f"The parameter is skipped and default value applies"
                    )
                    continue
