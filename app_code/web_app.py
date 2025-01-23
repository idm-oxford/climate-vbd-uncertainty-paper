"""Python script defining the web app."""

import pathlib

from climepi.app import run_app

if __name__ == "__main__":
    clim_dataset_example_base_dir = pathlib.Path(__file__).parents[1] / "data"
    clim_dataset_example_names = ["isimip_cities_daily"]
    epi_model_example_names = [
        "mordecai_ae_aegypti_niche",
        "mordecai_ae_albopictus_niche",
    ]
    run_app(
        clim_dataset_example_base_dir=clim_dataset_example_base_dir,
        clim_dataset_example_names=clim_dataset_example_names,
        enable_custom_epi_model=False,
        dask_distributed=True,
        address="0.0.0.0",
        port=7860,
        allow_websocket_origin=["*"],
    )
