"""
Module storing reusable test fixtures.

Extensive inspiration and lifting of definitions is taken from the test suite for the
xcdat package (see https://github.com/xCDAT/xcdat/blob/main/tests/fixtures.py)
"""

import cftime
import numpy as np
import xarray as xr

from climepi._xcdat import swap_lon_axis

# Time
time_yearly = xr.DataArray(
    data=np.array(
        [
            cftime.DatetimeNoLeap(2000, 7, 2, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2001, 7, 2, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2002, 7, 2, 12, 0, 0, 0, has_year_zero=True),
        ],
        dtype=object,
    ),
    dims=["time"],
    attrs={
        "axis": "T",
        "long_name": "time",
        "standard_name": "time",
    },
)
time_monthly = xr.DataArray(
    data=np.array(
        [
            cftime.DatetimeNoLeap(2000, 1, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 15, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 3, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 4, 16, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 5, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 6, 16, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 7, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 8, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 9, 16, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 10, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 11, 16, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 12, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2001, 1, 16, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2001, 2, 15, 0, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2001, 12, 16, 12, 0, 0, 0, has_year_zero=True),
        ],
        dtype=object,
    ),
    dims=["time"],
    attrs={
        "axis": "T",
        "long_name": "time",
        "standard_name": "time",
    },
)
time_daily = xr.DataArray(
    data=np.array(
        [
            cftime.DatetimeNoLeap(2000, 1, 28, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 1, 29, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 1, 30, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 1, 31, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 1, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 2, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 3, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 4, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 5, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 6, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 7, 12, 0, 0, 0, has_year_zero=True),
            cftime.DatetimeNoLeap(2000, 2, 8, 12, 0, 0, 0, has_year_zero=True),
        ],
        dtype=object,
    ),
    dims=["time"],
    attrs={
        "axis": "T",
        "long_name": "time",
        "standard_name": "time",
    },
)

time_bnds_yearly = xr.DataArray(
    name="time_bnds",
    data=np.array(
        [
            [
                cftime.DatetimeNoLeap(2000, 1, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2001, 1, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2001, 1, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2002, 1, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2002, 1, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2003, 1, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
        ],
        dtype=object,
    ),
    dims=["time", "bnds"],
    attrs={
        "xcdat_bounds": "True",
    },
)
time_bnds_monthly = xr.DataArray(
    name="time_bnds",
    data=np.array(
        [
            [
                cftime.DatetimeNoLeap(2000, 1, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 3, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 3, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 4, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 4, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 5, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 5, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 6, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 6, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 7, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 7, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 8, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 8, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 9, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 9, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 10, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 10, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 11, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 11, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 12, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 12, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2001, 1, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2001, 1, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2001, 2, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2001, 2, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2001, 3, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2001, 12, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2002, 1, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
        ],
        dtype=object,
    ),
    dims=["time", "bnds"],
    attrs={
        "xcdat_bounds": "True",
    },
)

time_bnds_daily = xr.DataArray(
    name="time_bnds",
    data=np.array(
        [
            [
                cftime.DatetimeNoLeap(2000, 1, 28, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 1, 29, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 1, 29, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 1, 30, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 1, 30, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 1, 31, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 1, 31, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 1, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 1, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 2, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 2, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 3, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 3, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 4, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 4, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 5, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 5, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 6, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 6, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 7, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 7, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 8, 0, 0, 0, 0, has_year_zero=True),
            ],
            [
                cftime.DatetimeNoLeap(2000, 2, 8, 0, 0, 0, 0, has_year_zero=True),
                cftime.DatetimeNoLeap(2000, 2, 9, 0, 0, 0, 0, has_year_zero=True),
            ],
        ],
        dtype=object,
    ),
    dims=["time", "bnds"],
    attrs={
        "xcdat_bounds": "True",
    },
)

# Latitude
lat = xr.DataArray(
    data=np.array([-90, -88.75, 88.75, 90]),
    dims=["lat"],
    attrs={"units": "degrees_north", "axis": "Y", "standard_name": "latitude"},
)
lat_bnds = xr.DataArray(
    name="lat_bnds",
    data=np.array([[-90, -89.375], [-89.375, 0.0], [0.0, 89.375], [89.375, 90]]),
    coords={"lat": lat},
    dims=["lat", "bnds"],
    attrs={"xcdat_bounds": "True"},
)

# Longitude
lon = xr.DataArray(
    data=np.array([0, 1.875, 356.25, 358.125]),
    dims=["lon"],
    attrs={"units": "degrees_east", "axis": "X", "standard_name": "longitude"},
)
lon_bnds = xr.DataArray(
    name="lon_bnds",
    data=np.array(
        [
            [0, 0.9375],
            [0.9375, 179.0625],
            [179.0625, 357.1875],
            [357.1875, 359.0625],
        ]
    ),
    coords={"lon": lon},
    dims=["lon", "bnds"],
    attrs={"xcdat_bounds": "True"},
)


# Dataset generation
def generate_dataset(
    data_var="temperature",
    dtype="float64",
    frequency="yearly",
    lon_0_360=True,
    extra_dims=None,
    has_bounds=True,
    random=True,
):
    """
    Generate a test dataset.

    Parameters
    ----------
    data_var : str or list of str, optional
        Name of the data variable(s) to include in the dataset. Default is
        "temperature".
    dtype : type, optional
        Data type of the data variable(s). Default is "float64"
    frequency : str, optional
        Frequency to compute the group average for (options are "yearly", "monthly" or
        "daily"). Default is "yearly".
    lon_0_360 : bool, optional
        Whether to use a longitude grid from 0 to 360 degrees. Default is True.
        If False, the grid ranges from -180 to 180 degrees.
    extra_dims : Hashable, sequence of Hashable, dict, or None, optional
        Extra dimensions to add to the dataset (used as an argument to xarray's
        DataArray.expand_dims method). Default is None.
    has_bounds : bool, optional
        Whether to include bounds in the dataset. Default is True.
    random : bool, optional
        Whether to fill the data variable(s) with random values. Default is True.
        If False, all values are set to 1.

    Returns
    -------
    xr.Dataset
        Test dataset.
    """
    if isinstance(data_var, str):
        data_var = [data_var]
    if frequency == "monthly":
        time = time_monthly
        time_bnds = time_bnds_monthly
    elif frequency == "yearly":
        time = time_yearly
        time_bnds = time_bnds_yearly
    elif frequency == "daily":
        time = time_daily
        time_bnds = time_bnds_daily
    else:
        raise ValueError(f"Invalid frequency: {frequency}")
    # Create the base dataset.
    da = xr.DataArray(
        data=np.ones((len(time), len(lat), len(lon)), dtype=dtype),
        coords={"time": time.copy(), "lat": lat.copy(), "lon": lon.copy()},
        dims=["time", "lat", "lon"],
    ).expand_dims(extra_dims)
    ds = xr.Dataset(
        data_vars={data_var_curr: da.copy() for data_var_curr in data_var},
    )
    # Randomize data if requested
    if random:
        for data_var_curr in data_var:
            vals_new = np.random.rand(*ds[data_var_curr].shape)
            if np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.bool):
                vals_new = np.round(vals_new)
            ds[data_var_curr].values = vals_new.astype(dtype)

    # Set time encoding
    ds["time"].encoding["calendar"] = "noleap"
    ds["time"].encoding["units"] = "days since 2000-01-01"
    ds["time"].encoding["dtype"] = np.dtype("float64")
    # Add bounds
    if has_bounds:
        ds["lat_bnds"] = lat_bnds.copy()
        ds["lon_bnds"] = lon_bnds.copy()
        ds["time_bnds"] = time_bnds.copy()
        ds["lat"].attrs["bounds"] = "lat_bnds"
        ds["lon"].attrs["bounds"] = "lon_bnds"
        ds["time"].attrs["bounds"] = "time_bnds"
    if not lon_0_360:
        ds = swap_lon_axis(ds, to=(-180, 180))
    return ds
