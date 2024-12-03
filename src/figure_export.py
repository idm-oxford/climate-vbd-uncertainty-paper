import bokeh.layouts as bl
import bokeh.plotting as bp
import holoviews as hv
import svgutils.transform as svgt
from bokeh.io import export_svg
from selenium import webdriver
from selenium.webdriver.chrome.service import Service as ChromeService
from webdriver_manager.chrome import ChromeDriverManager

import opts

WEBDRIVER_SERVICE = ChromeService(ChromeDriverManager().install())
WEBDRIVER_OPTIONS = webdriver.ChromeOptions()
WEBDRIVER_OPTIONS.add_argument("--headless")


def export_main_figure(panel_list):
    """
    Export the main figure to SVG by combining the panels, formatting them, adding
    panel labels, and then appending them to a template.
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
    ):
        _format_bokeh_panel(
            p,
            reverse_legend_entries=reverse_legend_entries,
            attrs=attrs,
        )
    col1 = bl.column(panel_list[0:3])
    col2 = bl.column(panel_list[3:])
    layout = bl.row([col1, bl.Spacer(width=50), col2])
    # Export combined panels to SV
    figure_dir = opts.get_opts()["figure_dir"]
    panels_save_path = figure_dir / "main_panels.svg"
    _export_figure(layout, save_path=panels_save_path)
    # Add panels to figure template and save again
    template_path = figure_dir / "main_template.svg"
    final_save_path = figure_dir / "main.svg"
    fig = svgt.fromfile(template_path)
    fig.append(svgt.fromfile(panels_save_path).getroot())
    fig.save(final_save_path)


def export_sensitivity_figure(panel_list):
    """
    Export the sensitivity figure to SVG by combining the panels, formatting them,
    and adding panel labels.
    """
    panel_list = _render_panels(panel_list)
    panel_label_list = ["A.", "B.", "C."]
    show_legend_list = [True, False, False]
    for p, show_legend, panel_label in zip(
        panel_list,
        show_legend_list,
        panel_label_list,
    ):
        _format_bokeh_panel(
            p,
            attrs={
                "frame_width": 720,
                "frame_height": 180,
                "legend.visible": show_legend,
                "title.text": panel_label,
                "title.offset": 772.5,
            },
        )
    layout = bl.column(panel_list)
    figure_dir = opts.get_opts()["figure_dir"]
    panels_save_path = figure_dir / "sensitivity_location_model.svg"
    _export_figure(layout, save_path=panels_save_path)


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
        "x_range.range_padding": 0.075,
        "x_range.group_padding": 0.5,
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
    with webdriver.Chrome(
        options=WEBDRIVER_OPTIONS, service=WEBDRIVER_SERVICE
    ) as driver:
        export_svg(bokeh_fig, filename=save_path, webdriver=driver)
