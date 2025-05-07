fig_names = [
    "main",
    "summary_location_species",
    "sensitivity_splines",
]


rule figures:
    input:
        expand(
            "figures/{fig_name}.{extension}",
            fig_name=fig_names,
            extension=["svg", "pdf", "png"],
        ),


rule figures_svg:
    input:
        expand(
            "figures/{fig_name}.svg",
            fig_name=fig_names,
        ),


rule main_fig_svg:
    input:
        "paper_code/opts.py",
        "paper_code/figure_export.py",
        "figures/main_template.svg",
        "paper_code/main.ipynb",
    output:
        "figures/main.svg",
    notebook:
        "paper_code/main.ipynb"


rule weather_data:
    input:
        "paper_code/weather_data.ipynb",
    output:
        "data/weather_2020.nc",
    notebook:
        "paper_code/weather_data.ipynb"


rule summary_fig_svg:
    input:
        "data/weather_2020.nc",
        "paper_code/opts.py",
        "paper_code/figure_export.py",
        "paper_code/summary_plot.py",
        "paper_code/summary_location_species.ipynb",
    output:
        "figures/summary_location_species.svg",
    notebook:
        "paper_code/summary_location_species.ipynb"


rule splines_fig_svg:
    input:
        "paper_code/opts.py",
        "paper_code/figure_export.py",
        "figures/main_template.svg",
        "paper_code/sensitivity_splines.ipynb",
    output:
        "figures/sensitivity_splines.svg",
    notebook:
        "paper_code/sensitivity_splines.ipynb"


rule fig_svg_to_pdf:
    input:
        "figures/{fig_name}.svg",
    output:
        "figures/{fig_name}.pdf",
    shell:
        "inkscape --export-type=pdf --export-filename={output} {input}"


rule fig_svg_to_png:
    input:
        "figures/{fig_name}.svg",
    output:
        "figures/{fig_name}.png",
    shell:
        "inkscape --export-type=png --export-filename={output} {input}"
