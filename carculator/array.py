"""
.. module: array.py

"""


import numpy as np
import xarray as xr


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

    for param in cip:
        for size in cip.metadata[param]["sizes"]:
            for powertrain in cip.metadata[param]["powertrain"]:
                array[
                    size_dict[size],
                    powertrain_dict[powertrain],
                    parameter_dict[cip.metadata[param]["name"]],
                    year_dict[cip.metadata[param]["year"]],
                    :,
                ] = cip.values[param]

    return (size_dict, powertrain_dict, parameter_dict, year_dict), array