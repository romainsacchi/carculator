import pandas as pd
import numpy as np


def summary_dataframe(array):
    """Calculate average and standard deviation of a parameter.

    ``array`` should be an ``xarray.DataArray`` with ``size`` and ``powertrain`` dimensions.
    You can select the relevant array with something like ``array.sel(parameter="driving mass", year=2017)``.

    Assumes that the last dimension is a set of data, e.g. Monte Carlo samples.

    Facets by powertrain and size."""
    data = [
        np.average(array.sel(size=size, powertrain=pt))
        for size in array["size"].data
        for pt in array["powertrain"].data
    ]
    index = pd.MultiIndex.from_tuples(
        [(size, pt) for size in array["size"].data for pt in array["powertrain"].data],
        names=["size", "powertrain"],
    )
    return pd.DataFrame(data, index=index)
