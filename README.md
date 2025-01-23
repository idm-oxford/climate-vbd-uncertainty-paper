---
title: climate-vbd-uncertainty-paper
emoji: 📈
colorFrom: blue
colorTo: yellow
sdk: docker
pinned: false
license: gpl-3.0
short_description: Impact of climate uncertainty on VBD suitability projections
---

Code accompanying manuscript.

Note that the author-developed `climepi` Python package
(distributed via `conda-forge`, with source code available at
https://github.com/idm-oxford/climate-epidemics) is used extensively in the manuscript
code.

In addition to code for reproducing the figures in the manuscript, this repository
provides source code for a web app available at
https://idm-oxford.github.io/climate-vbd-uncertainty-paper/, which is hosted on
[Hugging Face spaces](
https://huggingface.co/spaces/will-s-hart/climate-vbd-uncertainty).

To reproduce the figures, download the source code and
create a `conda` virtual environment with the required dependencies:
```
climate-vbd-uncertainty-paper $ conda env create -f environment.yml
```
Or alternatively, using `mamba` (may be faster):
```
climate-vbd-uncertainty-paper $ mamba env create -f environment.yml
```

The figures (in SVG format) can then be reproduced by executing the notebooks in the
``src`` directory using the ``climate-vbd-uncertainty`` environment created in the above
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
