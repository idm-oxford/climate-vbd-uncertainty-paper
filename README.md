---
sdk: docker
license: gpl-3.0
---

Python code and data accompanying manuscript “Internal climate variability amplifies the
need for vector-borne disease outbreak preparedness” by Hart *et al*.


Note that the author-developed `climepi` Python package
(distributed via `conda-forge`, with source code available at
https://github.com/idm-oxford/climate-epidemics) is used extensively in the manuscript
code.

This repository contains:
- Code for reproducing the figures in the manuscript
    ('paper_code' directory).
- [Climate projection](
    https://www.isimip.org/gettingstarted/terms-of-use/terms-use-publicly-available-isimip-data-after-embargo-period/)
    and [recorded temperature](https://dev.meteostat.net/terms.html) data used in the
    analyses ('data' directory).
- Results of the literature search reported in the Supporting Information of the
    manuscript ('literature_search' directory).
- Source code for a web app available at
    https://idm-oxford.github.io/climate-vbd-uncertainty-paper/, which is hosted on
    [Hugging Face spaces](
        https://huggingface.co/spaces/will-s-hart/climate-vbd-uncertainty) ('app_code'
    directory).

To reproduce the figures, download the source code and
create a ``conda`` virtual environment named 'climate-vbd-uncertainty' with the required
dependencies:
```
climate-vbd-uncertainty-paper $ conda env create -f environment.yml
```
Or alternatively, using `mamba` (may be faster):
```
climate-vbd-uncertainty-paper $ mamba env create -f environment.yml
```

The figures (in SVG format) can then be reproduced by executing the notebooks in the
``src`` directory using the 'climate-vbd-uncertainty' environment created in the above
step.

Alternatively, all figures can be reproduced using ``snakemake``:
```
climate-vbd-uncertainty-paper $ conda activate climate-vbd-uncertainty
(climate-vbd-uncertainty) climate-vbd-uncertainty-paper $ snakemake --cores 1 --forceall figures_svg
```

Note that the Google Chrome browser must be installed in order to export the figures to
SVG.

If [Inkscape](https://inkscape.org) is installed, PNG and PDF figures can be created by
running
```
(climate-vbd-uncertainty) climate-vbd-uncertainty-paper $ snakemake --cores 1
```
