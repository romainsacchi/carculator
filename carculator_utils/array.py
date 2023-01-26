import itertools

import numpy as np
import pandas as pd
import stats_arrays as sa
import xarray as xr

from .vehicle_input_parameters import VehicleInputParameters as vip


def fill_xarray_from_input_parameters(input_parameters, sensitivity=False, scope=None):
    """Create an `xarray` labeled array from the sampled input parameters.


    This function extracts the parameters' names and values contained in the
    `parameters` attribute of the :class:`CarInputParameters` class
    in :mod:`car_input_parameters` and insert them into a
    multi-dimensional numpy-like array from the *xarray* package
    (http://xarray.pydata.org/en/stable/).


    :param sensitivity:
    :param input_parameters: Instance of the :class:`TruckInputParameters` class
    in :mod:`truck_input_parameters`.
    :returns: `tuple`, `xarray.DataArray`
    - tuple (`size_dict`, `powertrain_dict`, `parameter_dict`, `year_dict`)
    - array

    Dimensions of `array`:

        0. Vehicle size, e.g. "3.5t", "7.5t", etc. str.
        1. Powertrain, e.g. "ICE-d", "BEV". str.
        2. Year. int.
        3. Samples.

    """

    # Check whether the argument passed is an instance of :class:`TruckInputParameters`
    if not isinstance(input_parameters, vip):
        raise TypeError(
            "The argument passed is not an object of the TruckInputParameter class"
        )

    if scope is None:
        scope = {
            "size": input_parameters.sizes,
            "powertrain": input_parameters.powertrains,
            "year": input_parameters.years,
        }
    else:
        if "size" not in scope:
            scope["size"] = input_parameters.sizes
        if "powertrain" not in scope:
            scope["powertrain"] = input_parameters.powertrains
        if "year" not in scope:
            scope["year"] = input_parameters.years

    # Make sure to include PHEV-e and PHEV-c-d if
    # PHEV-d is listed

    missing_pwts = [
        ("PHEV-d", "PHEV-e", "PHEV-c-d"),
        ("PHEV-p", "PHEV-e", "PHEV-c-p"),
    ]

    for missing_pwt in missing_pwts:
        if missing_pwt[0] in scope["powertrain"]:
            if not any(p in scope["powertrain"] for p in missing_pwt[1:]):
                scope["powertrain"].extend(missing_pwt[1:])

    if any(s for s in scope["size"] if s not in input_parameters.sizes):
        raise ValueError("One of the size types is not valid.")

    if any(y for y in scope["year"] if y not in input_parameters.years):
        raise ValueError("One of the years defined is not valid.")

    if any(pt for pt in scope["powertrain"] if pt not in input_parameters.powertrains):
        raise ValueError("One of the powertrain types is not valid.")

    # if the purpose is not to do a sensitivity analysis
    # the dimension `value` of the array is as large as
    # the number of iterations to perform
    # that is, 1 in `static` mode, or several in `stochastic` mode.

    d = {v: k for k, v in enumerate(scope["size"])}

    if not sensitivity:
        array = xr.DataArray(
            np.zeros(
                (
                    len(scope["size"]),
                    len(scope["powertrain"]),
                    len(input_parameters.parameters),
                    len(scope["year"]),
                    input_parameters.iterations or 1,
                )
            ),
            coords=[
                sorted(scope["size"], key=lambda x: d[x]),
                scope["powertrain"],
                input_parameters.parameters,
                scope["year"],
                np.arange(input_parameters.iterations or 1),
            ],
            dims=["size", "powertrain", "parameter", "year", "value"],
        ).astype("float32")
    # if the purpose is ot do a sensitivity analysis
    # then the length of the dimensions `value` equals the number of parameters
    else:
        params = ["reference"]
        params.extend([a for a in input_parameters.input_parameters])
        array = xr.DataArray(
            np.zeros(
                (
                    len(scope["size"]),
                    len(scope["powertrain"]),
                    len(input_parameters.parameters),
                    len(scope["year"]),
                    len(params),
                )
            ),
            coords=[
                scope["size"],
                scope["powertrain"],
                input_parameters.parameters,
                scope["year"],
                params,
            ],
            dims=["size", "powertrain", "parameter", "year", "value"],
        ).astype("float32")

    size_dict = {k: i for i, k in enumerate(scope["size"])}
    powertrain_dict = {k: i for i, k in enumerate(scope["powertrain"])}
    year_dict = {k: i for i, k in enumerate(scope["year"])}
    parameter_dict = {k: i for i, k in enumerate(input_parameters.parameters)}

    if not sensitivity:
        for param in input_parameters:
            pwt = (
                set(input_parameters.metadata[param]["powertrain"])
                if isinstance(input_parameters.metadata[param]["powertrain"], list)
                else set([input_parameters.metadata[param]["powertrain"]])
            )

            size = (
                set(input_parameters.metadata[param]["sizes"])
                if isinstance(input_parameters.metadata[param]["sizes"], list)
                else set([input_parameters.metadata[param]["sizes"]])
            )

            year = (
                set(input_parameters.metadata[param]["year"])
                if isinstance(input_parameters.metadata[param]["year"], list)
                else set([input_parameters.metadata[param]["year"]])
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
                        parameter=input_parameters.metadata[param]["name"],
                    )
                ] = input_parameters.values[param]

    else:
        # if `sensitivity` == True, the values of each parameter is
        # incremented by 10% when `value` == `parameter`
        for x, param in enumerate(input_parameters.input_parameters):
            names = [
                n
                for n in input_parameters.metadata
                if input_parameters.metadata[n]["name"] == param
            ]

            pwt = list(
                set(
                    itertools.chain.from_iterable(
                        [
                            input_parameters.metadata[name]["powertrain"]
                            for name in names
                        ]
                    )
                )
            )

            size = list(
                set(
                    itertools.chain.from_iterable(
                        [input_parameters.metadata[name]["sizes"] for name in names]
                    )
                )
            )

            year = [str(input_parameters.metadata[name]["year"]) for name in names]
            year = list(set(year))
            year = [int(y) for y in year]

            for name in names:
                vals = [
                    input_parameters.values[name]
                    for _ in range(0, len(input_parameters.input_parameters) + 1)
                ]
                vals[input_parameters.input_parameters.index(param) + 1] *= 1.1

                array.loc[
                    dict(
                        powertrain=[p for p in pwt if p in scope["powertrain"]],
                        size=[s for s in size if s in scope["size"]],
                        year=[y for y in year if y in scope["year"]],
                        parameter=input_parameters.metadata[name]["name"],
                    )
                ] = vals

    return (size_dict, powertrain_dict, parameter_dict, year_dict), array
