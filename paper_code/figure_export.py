"""Module defining methods for exporting bokeh/holoviews figures to SVG."""

import bokeh.layouts as bl
import bokeh.models as bm
import bokeh.plotting as bp
import holoviews as hv
import opts
import svgutils.transform as svgt
from bokeh.io import export_svg
from lxml import etree
from selenium.webdriver import Chrome, ChromeOptions
from selenium.webdriver.chrome.service import Service as ChromeService
from webdriver_manager.chrome import ChromeDriverManager

WEBDRIVER_SERVICE = ChromeService(ChromeDriverManager().install())
WEBDRIVER_OPTIONS = ChromeOptions()
WEBDRIVER_OPTIONS.add_argument("--headless")


def export_main_figure(panel_list, file_name="main.svg"):
    """
    Export the main figure to SVG.

    Combines the panels, formats them, adds panel labels, and appends panels to a
    template, before exporting the figure.
    """
    panel_list = _render_panels(panel_list)
    panel_label_list = [
        "Internal climate variability",
        "Climate model uncertainty",
        "Scenario uncertainty",
        "Climate projection uncertainty",
        "Epidemiological projection uncertainty",
    ]
    attrs_list = [
        {
            "title.text": ll,
            "title.align": "left",
            "title.offset": -50,
            "title.standoff": 15,
        }
        for ll in panel_label_list
    ]
    for i in range(3):
        attrs_list[i] = {
            **attrs_list[i],
            "frame_height": 182,
            "title.text_font_size": "14pt",
            "yaxis.axis_label": "Annual mean temp. (°C)",
        }
    attrs_list[4]["legend.visible"] = False
    reverse_legend_entries_list = [True, False, False, True, None]

    for p, reverse_legend_entries, attrs in zip(
        panel_list,
        reverse_legend_entries_list,
        attrs_list,
        strict=True,
    ):
        _format_bokeh_panel(
            p,
            reverse_legend_entries=reverse_legend_entries,
            attrs=attrs,
        )
    col1 = bl.column(panel_list[0:3])
    col2 = bl.column(panel_list[3:])
    layout = bl.row([col1, bl.Spacer(width=50), col2])
    # Export combined panels to SVG
    figure_dir = opts.get_opts()["figure_dir"]
    panels_save_path = figure_dir / ("panels_" + file_name)
    _export_figure(layout, save_path=panels_save_path)
    # Add panels to figure template and save again
    template_path = figure_dir / "main_template.svg"
    final_save_path = figure_dir / file_name
    fig = svgt.fromfile(template_path)
    fig.append(svgt.fromfile(panels_save_path).getroot())
    fig.save(final_save_path)


def export_summary_figure(panel_list):
    """
    Export the summary figure to SVG.

    Combines the panels, formats them, and adds panel labels, before exporting the
    figure.
    """
    panel_list = _render_panels(panel_list)
    panel_label_list = ["A.", "B.", "C."]
    show_legend_list = [True, False, False]
    for p, show_legend, panel_label in zip(
        panel_list,
        show_legend_list,
        panel_label_list,
        strict=True,
    ):
        _format_bokeh_panel(
            p,
            attrs={
                "frame_width": 720,
                "frame_height": 180,
                "legend.visible": show_legend,
                # "legend.spacing": 1,
                # "legend.margin": 5,
                # "legend.padding": 3,
                # "legend.location": (0, -10),
                "title.text": panel_label,
                "title.offset": 777.5,
                "x_range.group_padding": 0.25,
                "x_range.range_padding": 0.1,
                "x_range.range_padding_units": "absolute",
            },
        )
    # Hack to have legend outside of the plot area for the first panel
    p = panel_list[0]
    p.legend.visible = False
    p.add_layout(
        bm.Legend(
            items=p.legend.items,
            location=(-710, 77.5),
            label_text_font_size="10pt",
            margin=5,
            padding=5,
        ),
        "right",
    )
    # Export combined panels to SVG
    layout = bl.column(panel_list)
    figure_dir = opts.get_opts()["figure_dir"]
    save_path = figure_dir / "summary_location_species.svg"
    _export_figure(layout, save_path=save_path)
    _italicize_mosquito_species(save_path)
    # Resize to remove whitespace from putting legend in right side panel
    _resize_svg(save_path, width=810)


def export_sensitivity_figure(panel_list):
    """
    Export the sensitivity figure to SVG.

    Combines the panels, formats them, and adds panel labels, before exporting the
    figure.
    """
    panel_list = _render_panels(panel_list)
    panel_label_list = ["A.", "B."]
    for p, panel_label in zip(
        panel_list,
        panel_label_list,
        strict=True,
    ):
        _format_bokeh_panel(
            p,
            reverse_legend_entries=True,
            attrs={
                "title.text": panel_label,
                "title.standoff": 2.5,
            },
        )
    layout = bl.row(panel_list)
    figure_dir = opts.get_opts()["figure_dir"]
    save_path = figure_dir / "sensitivity_niche.svg"
    _export_figure(layout, save_path=save_path)


def _render_panels(panel_list):
    # Render HoloViews objects to Bokeh figures
    panel_list_rendered = []
    for p in panel_list:
        if isinstance(p, bp.figure):
            panel_list_rendered.append(p)
        else:
            panel_list_rendered.append(hv.render(p))
    return panel_list_rendered


def _format_bokeh_panel(
    bokeh_fig,
    reverse_legend_entries=False,
    attrs=None,
):
    attrs = {
        "output_backend": "svg",
        "width": None,
        "height": None,
        "frame_width": 300,
        "frame_height": 300,
        "background_fill_color": None,
        "border_fill_color": None,
        "xaxis.axis_label_text_font_style": "normal",
        "yaxis.axis_label_text_font_style": "normal",
        "xaxis.axis_label_text_font_size": "12pt",
        "yaxis.axis_label_text_font_size": "12pt",
        "xaxis.major_label_text_font_size": "10pt",
        "yaxis.major_label_text_font_size": "10pt",
        "xaxis.group_text_font_size": "10pt",
        "legend.label_text_font_size": "10pt",
        "legend.margin": 5,
        "legend.padding": 5,
        "legend.title": None,
        "title.text": "",
        "title.text_font_size": "16pt",
        "title.text_font_style": "normal",
        "title.align": "right",
        "title.offset": 335,
        "title.standoff": -5,
        "toolbar_location": None,
        "min_border_left": 57,  # account for offset panel labels
        "min_border_right": 20,  # enough to avoid x-axis labels being squeezed
        **attrs,
    }
    for k, v in attrs.items():
        try:
            if "." in k:
                k_split = k.split(".")
                setattr(getattr(bokeh_fig, k_split[0]), k_split[1], v)
            else:
                setattr(bokeh_fig, k, v)
        except AttributeError:
            pass
    if reverse_legend_entries:
        bokeh_fig.legend.items.reverse()
    return bokeh_fig


def _export_figure(bokeh_fig, save_path):
    with Chrome(options=WEBDRIVER_OPTIONS, service=WEBDRIVER_SERVICE) as driver:
        export_svg(bokeh_fig, filename=save_path, webdriver=driver)


def _italicize_mosquito_species(svg_path):
    # Italicise mosquito species names in an SVG file (needed as Bokeh does not support
    # mixed font styles in text labels)
    fig = svgt.fromfile(svg_path)
    fig.set_size((fig.width, fig.height))
    for element in fig.root.findall(".//{http://www.w3.org/2000/svg}text"):
        text = element.text
        if text in ["(Ae. albopictus)", "(Ae. aegypti)"]:
            species_name = text[1:-1]
            # Clear the current text
            element.text = ""
            # Add a non-italic <tspan> for the opening bracket
            tspan_prefix = etree.Element("tspan")
            tspan_prefix.text = "("
            element.append(tspan_prefix)
            # Add an italic <tspan> for the mosquito name
            tspan_name = etree.Element("tspan")
            tspan_name.text = species_name
            tspan_name.set("font-style", "italic")
            element.append(tspan_name)
            # Add a non-italic <tspan> for the closing bracket
            tspan_suffix = etree.Element("tspan")
            tspan_suffix.text = ")"
            if species_name == "Ae. aegypti":
                # Adjust the x-offset for the closing bracket to avoid overlap
                tspan_suffix.set("dx", "1.5")
            element.append(tspan_suffix)
    fig.save(svg_path)


def _resize_svg(svg_path, width=None, height=None):
    # Resize an SVG file
    fig = svgt.fromfile(svg_path)
    width = str(width or fig.width)
    height = str(height or fig.height)
    fig.set_size((width, height))
    fig.save(svg_path)
