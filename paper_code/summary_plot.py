"""Module for making summary plots of uncertainty in different locations and years."""

import holoviews as hv
import xarray as xr
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure


def make_summary_plot(
    ds,
    ds_historical=None,
    data_var=None,
    years_summary=None,
    ylim=None,
    yticks=None,
    ylabel="",
    legend_location="top_left",
    **kwargs_uncertainty_interval_decomposition,
):
    """
    Make a summary plot.

    Makes a grouped bar plot showing contributions to uncertainty from internal climate
    variability, model uncertainty, and scenario uncertainty in different locations
    and years.
    """
    ds_plot = (
        ds.climepi.uncertainty_interval_decomposition(
            data_var, **kwargs_uncertainty_interval_decomposition
        )[data_var]
        .sel(time=ds["time.year"].isin(years_summary))
        .to_dataset(dim="level")
        .compute()
    )
    ds_plot = ds_plot.assign_coords(time=ds_plot["time.year"].astype("str"))

    #

    ds_plot_historical = (
        ds_historical[[data_var]]
        .assign_coords(time=ds_historical.time.dt.year.astype("str"))
        .rename({data_var: "historical"})
    )
    ds_plot = xr.merge([ds_plot, ds_plot_historical])

    #

    source = ColumnDataSource(ds_plot.to_dataframe(dim_order=["location", "time"]))
    colors = hv.Cycle().values
    p = figure(
        x_range=FactorRange(*source.data["location_time"]), width=700, height=300
    )
    p.line(  # Dummy line to get legend
        x=[p.x_range.start],
        y=[p.y_range.start],
        line_color="black",
        line_width=2,
        legend_label="Mean",
    )
    for bottom, top, color, legend_label in (
        ("internal_lower", "internal_upper", colors[0], "Internal variability"),
        ("model_lower", "internal_lower", colors[1], "Model uncertainty"),
        ("internal_upper", "model_upper", colors[1], None),
        ("scenario_lower", "model_lower", colors[2], "Scenario uncertainty"),
        ("model_upper", "scenario_upper", colors[2], None),
    ):
        kwargs_vbar = {
            "x": "location_time",
            "source": source,
            "bottom": bottom,
            "top": top,
            "fill_color": color,
            "line_width": 0.5,
            "width": 0.8,
        }
        if legend_label is not None:
            kwargs_vbar["legend_label"] = legend_label
        p.vbar(**kwargs_vbar)
    p.vbar(
        x="location_time",
        source=source,
        bottom="baseline",
        top="baseline",
        line_color="black",
        line_width=2,
        width=0.8,
    )
    p.scatter(
        x="location_time",
        y="historical",
        source=source,
        color="black",
        marker="x",
        size=15,
        line_width=2,
    )
    # Formatting
    if ylim is not None:
        p.y_range.start, p.y_range.end = ylim
    if yticks is not None:
        p.yaxis.ticker = yticks
    p.yaxis.axis_label = ylabel
    p.legend.location = legend_location
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    return p
