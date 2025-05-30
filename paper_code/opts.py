"""Module defining data, model and plotting options."""

import pathlib

import numpy as np


def get_opts():
    """Return a dictionary of data, model and plotting options."""
    year_range = [2030, 2100]  # inclusive
    location_default = "London"
    locations_additional = ["Paris", "Istanbul", "Cape Town", "Los Angeles"]
    epi_model_name_default = "mordecai_ae_albopictus_niche"
    epi_model_species_default = "Ae. albopictus"
    epi_model_name_additional = "mordecai_ae_aegypti_niche"
    epi_model_species_additional = "Ae. aegypti"
    ensemble_stats_kwargs = {
        "internal_variability_method": "polyfit",
        "deg": 3,
        "uncertainty_level": 90,
    }
    plot_opts_base = {
        "title": "",
        "xlim": (year_range[0], year_range[1]),
        "xticks": np.arange(year_range[0], year_range[1] + 1, 10),
        "xlabel": "",
        "legend_position": "top_left",
    }
    plot_opts_clim = {
        **plot_opts_base,
        "ylim": (9, 18),
        "yticks": np.arange(9, 19, 3),
        "ylabel": "Annual mean temperature (°C)",
    }
    plot_opts_epi = {
        **plot_opts_base,
        "ylim": (0, 210),
        "yticks": np.arange(211, step=30),
        "ylabel": "Days suitable for transmission",
    }
    data_base_dir = pathlib.Path(__file__).parents[1] / "data"
    figure_dir = pathlib.Path(__file__).parents[1] / "figures"
    opts = {
        "year_range": year_range,
        "location_default": location_default,
        "locations_additional": locations_additional,
        "epi_model_name_default": epi_model_name_default,
        "epi_model_species_default": epi_model_species_default,
        "epi_model_name_additional": epi_model_name_additional,
        "epi_model_species_additional": epi_model_species_additional,
        "ensemble_stats_kwargs": ensemble_stats_kwargs,
        "plot_opts_clim": plot_opts_clim,
        "plot_opts_epi": plot_opts_epi,
        "data_base_dir": data_base_dir,
        "figure_dir": figure_dir,
    }
    return opts
