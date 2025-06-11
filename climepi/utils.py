"""Module containing utility functions for working with xarray datasets."""


def add_var_attrs_from_other(ds, ds_from, var=None):
    """
    Copy variable attributes from one xarray dataset to another.

    Copies the attributes for a variable (or variables) from another xarray dataset
    (ds_from) to a copy of another one (ds), whenever the variable(s) exist in both
    datasets.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to add the variable attributes to.
    ds_from : xarray.Dataset
        The dataset to add the variable attributes from.
    var : str or list, optional
        The name(s) of the variable(s) to add the attributes for (provided both datasets
        contain variable(s) with these names). If None, all variables in ds are used.

    Returns
    -------
    xarray.Dataset
        A copy of ds with the relevant variable attributes from ds_from added.
    """
    ds_out = ds.copy()
    if var is None:
        var = list(ds.data_vars) + list(ds.coords)
    elif isinstance(var, str):
        var = [var]
    for var_curr in var:
        try:
            ds_out[var_curr].attrs = ds_from[var_curr].attrs
        except KeyError:
            pass
    return ds_out


def add_bnds_from_other(ds, ds_from):
    """
    Add latitude, longitude, and time bounds from one xarray dataset to another.

    Adds the latitude, longitude, and time bounds from one xarray dataset (ds_from) to a
    copy of another one (ds), whenever the bounds exist in ds_from but not ds, and the
    corresponding co-ordinate dimension is the same for both datasets.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to add the bounds to.
    ds_from : xarray.Dataset
        The dataset to get the bounds from.

    Returns
    -------
    xarray.Dataset
        A copy of ds with the relevant bounds variables from ds_from added.
    """
    ds_out = ds.copy()
    for var in ["lat", "lon", "time"]:
        bnd_var = var + "_bnds"
        if (
            bnd_var in ds_out.data_vars
            or bnd_var not in ds_from
            or var not in ds_out
            or var not in ds_from
            or not ds_out[var].equals(ds_from[var])
        ):
            continue
        ds_out[bnd_var] = ds_from[bnd_var]
        ds_out[bnd_var].attrs = ds_from[bnd_var].attrs
        ds_out[var].attrs.update(bounds=bnd_var)
    return ds_out


def get_data_var_and_bnds(ds, data_var):
    """
    Get a dataset with only the selected data variable(s) and any bounds variables.

    Returns a new dataset containing only the selected data variable(s) and any
    bounds variables.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to select the data variable(s) from.
    data_var : str or list, optional
        Name(s) of the data variable(s) to select.

    Returns
    -------
    xarray.Dataset
        A new dataset containing the selected data variable(s) and any bounds
        variables.
    """
    if isinstance(data_var, str):
        data_var_list = [data_var]
    elif isinstance(data_var, list):
        data_var_list = data_var
    else:
        raise ValueError("data_var must be a string or list")
    for bnd_var in ["lat_bnds", "lon_bnds", "time_bnds"]:
        if bnd_var in ds:
            data_var_list.append(bnd_var)
    ds_out = ds[data_var_list]
    return ds_out


def list_non_bnd_data_vars(ds):
    """
    List the names of the non-bound variables in the dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to list the non-bound variables from.

    Returns
    -------
    list
        Names of the non-bound variables in the dataset.
    """
    data_vars = list(ds.data_vars)
    bnd_vars = ["lat_bnds", "lon_bnds", "time_bnds"]
    non_bnd_data_vars = [
        data_vars[i] for i in range(len(data_vars)) if data_vars[i] not in bnd_vars
    ]
    return non_bnd_data_vars
