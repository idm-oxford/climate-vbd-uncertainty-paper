fig_names = [
    "main",
    "variance_plots",
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


rule fig_svg:
    input:
        "src/opts.py",
        "src/utils.py",
        "figures/main_template.svg",
        "src/{fig_name}.ipynb",
    output:
        "figures/{fig_name}.svg",
    notebook:
        "src/{wildcards.fig_name}.ipynb"


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
