# climate-vbd-uncertainty-paper
Code accompanying manuscript. To reproduce the figures, download the source code and
create a `conda` virtual environment with the required dependencies:
```
climate-disease-uncertainty $ conda env create -f environment.yml
```
Or alternatively, using `mamba` (may be faster):
```
climate-disease-uncertainty $ mamba env create -f environment.yml
```

The figures (in SVG format) can then be reproduced by executing the notebooks in the
``src`` directory using the ``climate-disease-uncertainty`` environment created in the above
step.

Alternatively, all figures can be reproduced using ``snakemake``:
```
climate-disease-uncertainty $ conda activate climate-disease-uncertainty
(climate-disease-uncertainty) climate-disease-uncertainty $ snakemake --cores 1 --forceall figures_svg
```

Note that the Google Chrome browser must be installed in order to export the figures to
SVG.

If [Inkscape](https://inkscape.org) is installed, PNG and PDF figures can be created by
running
```
(climate-disease-uncertainty) climate-disease-uncertainty $ snakemake --cores 1
```
