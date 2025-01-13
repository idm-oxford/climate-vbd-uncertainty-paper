fig_names = [
    "main",
    "summary_location_species",
    "sensitivity_niche",
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


rule summary_fig_svg:
    input:
        "src/opts.py",
        "src/figure_export.py",
        "src/summary_plot.py",
        "src/summary_location_species.ipynb",
    output:
        "figures/summary_location_species.svg",
    notebook:
        "src/summary_location_species.ipynb"


rule sensitivity_fig_svg:
    input:
        "src/opts.py",
        "src/figure_export.py",
        "src/sensitivity_niche.ipynb",
    output:
        "figures/sensitivity_niche.svg",
    notebook:
        "src/sensitivity_niche.ipynb"


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
