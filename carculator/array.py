"""
.. module: array.py

"""


import numpy as np
import xarray as xr
import pandas as pd
import stats_arrays as sa
import carculator.car_input_parameters as c_i_p


def fill_xarray_from_input_parameters(cip):

    """Create an `xarray` labeled array from the sampled input parameters.


    This function extracts the parameters' names and values contained in the
    `parameters` attribute of the CarInputParameters class and insert them into a
    multi-dimensional numpy-like array from the *xarray* package
    (http://xarray.pydata.org/en/stable/).


    :param cip: Instance of the CarInputParameters class.
    :returns: `tuple`, `array`
    - tuple (`size_dict`, `powertrain_dict`, `parameter_dict`, `year_dict`)
    - array

    Dimensions of `array`:

    0. Vehicle size, e.g. "small", "medium". str.
    1. Powertrain, e.g. "ICE-p", "BEV". str.
    2. Year. int.
    3. Samples.


    :rtype size_dict: dict
    :rtype powertrain_dict: dict
    :rtype parameter_dict: dict
    :rtype year_dict: dict
    :rtype array: xarray.DataArray

    """

    # Check whether the argument passed is an cip object
    if not isinstance(cip, c_i_p.CarInputParameters):
        raise TypeError(
            "The argument passed is not an object of the CarInputParameter class"
        )

    array = xr.DataArray(
        np.zeros(
            (
                len(cip.sizes),
                len(cip.powertrains),
                len(cip.parameters),
                len(cip.years),
                cip.iterations or 1,
            )
        ),
        coords=[
            cip.sizes,
            cip.powertrains,
            cip.parameters,
            cip.years,
            np.arange(cip.iterations or 1),
        ],
        dims=["size", "powertrain", "parameter", "year", "value"],
    )

    size_dict = {k: i for i, k in enumerate(cip.sizes)}
    powertrain_dict = {k: i for i, k in enumerate(cip.powertrains)}
    year_dict = {k: i for i, k in enumerate(cip.years)}
    parameter_dict = {k: i for i, k in enumerate(cip.parameters)}

    # TODO: this is taking a lot of time... needs to be optimized.
    # Fixed!
    for param in cip:
        array.loc[
            dict(
                powertrain=cip.metadata[param]["powertrain"],
                size=cip.metadata[param]["sizes"],
                year=cip.metadata[param]["year"],
                parameter=cip.metadata[param]["name"],
            )
        ] = cip.values[param]

    return (size_dict, powertrain_dict, parameter_dict, year_dict), array


def modify_xarray_from_custom_parameters(fp, array):
    """Override default parameters values in `xarray` based on values provided by the user.

        This function allows to override one or several default parameter values by providing either:
        * a file path to an Excel workbook that contains the new values
        * or a dictionary

        The dictionary must be of the following format:
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
        :param fp: File path of workbook with new values or dictionary.
        :type fp: str or dict

        """

    if isinstance(fp, str):
        try:
            d = pd.read_excel(
                fp,
                header=[0, 1],
                index_col=[0, 1, 2, 3, 4],
                sheet_name="Custom_parameters",
            ).to_dict(orient="index")
        except:
            raise FileNotFoundError("Custom parameters file not found.")
    elif isinstance(fp, dict):
        d = fp
    else:
        print("The format passed as parameter is not valid.")
        raise

    for k in d:
        if k[1].lower() != "all":
            pt = [p.strip() for p in k[1].split(",")]
            pt = [p for p in pt if p]
        else:
            pt = array.coords["powertrain"].values

        if k[2].lower() != "all":
            sizes = [s.strip() for s in k[2].split(",") if s]
        else:
            sizes = array.coords["size"].values

        param = k[3]

        if not param in array.coords["parameter"].values:
            print("{} is not a recognized parameter. It will be skipped.".format(param))
            continue

        val = d[k]

        distr_dic = {
            "triangular": 5,
            "lognormal": 2,
            "normal": 3,
            "uniform": 4,
            "none": 1,
        }
        distr = distr_dic[k[4]]

        year = set([v[0] for v in val])

        # Stochastic mode
        if array.sizes["value"] > 1:

            for y in year:

                # No uncertainty parameters given
                if distr == 1:
                    # There should be at least a `loc`
                    if ~np.isnan(val[(y, "loc")]):
                        for s in sizes:
                            for p in pt:
                                array.loc[
                                    dict(powertrain=p, size=s, year=y, parameter=param)
                                ] = val[(y, "loc")]
                    # Otherwise warn
                    else:
                        print("`loc`parameter missing for {} in {}.".format(param, y))
                        continue

                elif distr in [2, 3, 4, 5]:

                    # Check if the correct parameters are present
                    # Triangular

                    if distr == 5:
                        if (
                            np.isnan(val[(y, "loc")])
                            or np.isnan(val[(y, "minimum")])
                            or np.isnan(val[(y, "maximum")])
                        ):
                            print(
                                "One or more parameters for the triangular distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Lognormal
                    if distr == 2:
                        if np.isnan(val[(y, "loc")]) or np.isnan(val[(y, "scale")]):
                            print(
                                "One or more parameters for the lognormal distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Normal
                    if distr == 3:
                        if np.isnan(val[(y, "loc")]) or np.isnan(val[(y, "scale")]):
                            print(
                                "One or more parameters for the normal distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Uniform
                    if distr == 4:
                        if np.isnan(val[(y, "minimum")]) or np.isnan(
                            val[(y, "maximum")]
                        ):
                            print(
                                "One or more parameters for the uniform distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    a = sa.UncertaintyBase.from_dicts(
                        {
                            "loc": val[y, "loc"],
                            "scale": val[y, "scale"],
                            "shape": val[y, "shape"],
                            "minimum": val[y, "minimum"],
                            "maximum": val[y, "maximum"],
                            "uncertainty_type": distr,
                        }
                    )

                    rng = sa.MCRandomNumberGenerator(a)

                    for s in sizes:
                        for p in pt:
                            array.loc[
                                dict(powertrain=p, size=s, year=y, parameter=param)
                            ] = rng.generate(array.sizes["value"]).reshape((-1,))

                else:
                    print(
                        "The uncertainty type is not recognized for {} in {}.\n The parameter is skipped and default value applies".format(
                            param, y
                        )
                    )
                    continue

        # Static mode
        else:
            for y in year:
                if distr == 1:
                    # There should be at least a `loc`
                    if ~np.isnan(val[(y, "loc")]):
                        for s in sizes:
                            for p in pt:
                                array.loc[
                                    dict(powertrain=p, size=s, year=y, parameter=param)
                                ] = val[(y, "loc")]
                    # Otherwise warn
                    else:
                        print("`loc`parameter missing for {} in {}.".format(param, y))
                        continue

                elif distr in [2, 3, 4, 5]:

                    # Check if the correct parameters are present
                    # Triangular

                    if distr == 5:
                        if (
                            np.isnan(val[(y, "loc")])
                            or np.isnan(val[(y, "minimum")])
                            or np.isnan(val[(y, "maximum")])
                        ):
                            print(
                                "One or more parameters for the triangular distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Lognormal
                    if distr == 2:
                        if np.isnan(val[(y, "loc")]) or np.isnan(val[(y, "scale")]):
                            print(
                                "One or more parameters for the lognormal distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Normal
                    if distr == 3:
                        if np.isnan(val[(y, "loc")]) or np.isnan(val[(y, "scale")]):
                            print(
                                "One or more parameters for the normal distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    # Uniform
                    if distr == 4:
                        if np.isnan(val[(y, "minimum")]) or np.isnan(
                            val[(y, "maximum")]
                        ):
                            print(
                                "One or more parameters for the uniform distribution is/are missing for {} in {}.\n The parameter is skipped and default value applies".format(
                                    param, y
                                )
                            )
                            continue

                    a = sa.UncertaintyBase.from_dicts(
                        {
                            "loc": val[y, "loc"],
                            "scale": val[y, "scale"],
                            "shape": val[y, "shape"],
                            "minimum": val[y, "minimum"],
                            "maximum": val[y, "maximum"],
                            "uncertainty_type": distr,
                        }
                    )

                    dist = sa.uncertainty_choices[distr]
                    median = float(dist.ppf(a, np.array((0.5,))))

                    for s in sizes:
                        for p in pt:
                            array.loc[
                                dict(powertrain=p, size=s, year=y, parameter=param)
                            ] = median

                else:
                    print(
                        "The uncertainty type is not recognized for {} in {}.\n The parameter is skipped and default value applies".format(
                            param, y
                        )
                    )
                    continue
