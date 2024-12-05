fig_names = [
    "main",
    "sensitivity_location_model",
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
        "src/opts.py",
        "src/figure_export.py",
        "figures/main_template.svg",
        "src/main.ipynb",
    output:
        "figures/main.svg",
    notebook:
        "src/main.ipynb"


rule sensitivity_fig_svg:
    input:
        "src/opts.py",
        "src/figure_export.py",
        "src/summary_plot.py",
        "src/sensitivity_location_model.ipynb",
    output:
        "figures/sensitivity_location_model.svg",
    notebook:
        "src/sensitivity_location_model.ipynb"


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
