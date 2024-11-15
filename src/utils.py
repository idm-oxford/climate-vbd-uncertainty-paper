import logging
import pathlib
from math import ceil

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


def notebook_config():
    # Supress default INFO logging
    logging.getLogger().setLevel(logging.CRITICAL)


def get_data_base_path():
    return pathlib.Path(__file__).parents[1] / "data"


def export_figure(
    save_path,
    panel_list,
    tiling=None,
    show_legend_list=None,
    reverse_legend_entries_list=None,
    legend_location_list=None,
):
    panel_list = _render_panels(panel_list)
    if show_legend_list is None:
        show_legend_list = [True] * len(panel_list)
    if reverse_legend_entries_list is None:
        reverse_legend_entries_list = [False] * len(panel_list)
    if legend_location_list is None:
        legend_location_list = [None] * len(panel_list)
    no_panels = len(panel_list)
    if tiling is not None:
        rows, cols = tiling
    else:
        rows = 1 + (no_panels - 1) // 3
        cols = ceil(no_panels // rows)
    panel_label_list = [chr(65 + i) + "." for i in range(no_panels)]
    for p, show_legend, reverse_legend_entries, legend_location, panel_label in zip(
        panel_list,
        show_legend_list,
        reverse_legend_entries_list,
        legend_location_list,
        panel_label_list,
    ):
        _format_bokeh_panel(
            p,
            show_legend=show_legend,
            reverse_legend_entries=reverse_legend_entries,
            legend_location=legend_location,
            panel_label=panel_label,
        )
    grid = bl.gridplot(
        panel_list,
        ncols=cols,
        merge_tools=False,
        toolbar_location=None,
    )
    grid.width = 375 * cols
    grid.height = 345 * rows
    _export_figure(grid, save_path)


def export_main_figure(panel_list):
    panel_list = _render_panels(panel_list)
    panel_label_list = [
        "Internal climate variability",
        "Climate model uncertainty",
        "Scenario uncertainty",
        "Climate projection uncertainty",
        "Epidemiological projection uncertainty",
    ]
    show_legend_list = [True, True, True, True, False]
    reverse_legend_entries_list = [True, False, False, True, None]
    legend_location_list = [None, None, None, (7, 184), None]
    for p, show_legend, reverse_legend_entries, legend_location, panel_label in zip(
        panel_list,
        show_legend_list,
        reverse_legend_entries_list,
        legend_location_list,
        panel_label_list,
    ):
        _format_bokeh_panel(
            p,
            show_legend=show_legend,
            reverse_legend_entries=reverse_legend_entries,
            legend_location=legend_location,
            panel_label=panel_label,
        )
        p.title.align = "left"
        p.title.offset = -50
        p.title.standoff = 15
    for p in panel_list[0:3]:
        p.frame_height = 182
        p.title.text_font_size = "14pt"
        p.yaxis.axis_label = "Annual mean temp. (°C)"
        p.legend.margin = 5
        p.legend.padding = 5
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
    panel_list = _render_panels(panel_list)
    panel_label_list = ["A.", "B.", "C."]
    show_legend_list = [True] * len(panel_list)
    legend_location_list = [None, None, None]
    for p, show_legend, legend_location, panel_label in zip(
        panel_list,
        show_legend_list,
        legend_location_list,
        panel_label_list,
    ):
        _format_bokeh_panel(
            p,
            show_legend=show_legend,
            legend_location=legend_location,
            panel_label=panel_label,
        )
        p.title.offset = 755
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
    show_legend=True,
    reverse_legend_entries=False,
    legend_location=None,
    panel_label="",
):
    bokeh_fig.output_backend = "svg"
    bokeh_fig.background_fill_color = None
    bokeh_fig.border_fill_color = None
    bokeh_fig.xaxis.axis_label_text_font_style = "normal"
    bokeh_fig.yaxis.axis_label_text_font_style = "normal"
    if show_legend:
        bokeh_fig.legend.title = None
        if reverse_legend_entries:
            bokeh_fig.legend.items.reverse()
        if legend_location is not None:
            bokeh_fig.legend.location = legend_location
    else:
        bokeh_fig.legend.visible = False
    bokeh_fig.title.text = panel_label
    bokeh_fig.title.text_font_size = "16pt"
    bokeh_fig.title.text_font_style = "normal"
    bokeh_fig.title.align = "right"
    bokeh_fig.title.offset = 335
    bokeh_fig.title.standoff = -5
    bokeh_fig.toolbar_location = None
    bokeh_fig.min_border_left = 57  # account for offset panel labels
    bokeh_fig.min_border_right = 20  # enough to avoid x-axis labels being squeezed
    return bokeh_fig


def _export_figure(bokeh_fig, save_path):
    with webdriver.Chrome(
        options=WEBDRIVER_OPTIONS, service=WEBDRIVER_SERVICE
    ) as driver:
        export_svg(bokeh_fig, filename=save_path, webdriver=driver)
