# climate-disease-uncertainty
Code underlying manuscript. To reproduce the figures, download the source code and
create a `conda` virtual environment with the required packages:
```
climate-disease-uncertainty $ conda env create -f environment.yml
```
Or alternatively, using `mamba` (which is likely to lead to much faster environment
creation):
```
climate-disease-uncertainty $ mamba env create -f environment.yml
```

The figures (in svg format) can then be reproduced by executing the notebooks in the
`src` directory using the `climate-disease-uncertainty` environment created in the above
step.

Alternatively, all figures can be reproduced using `snakemake`:
```
climate-disease-uncertainty $ conda activate climate-disease-uncertainty
(climate-disease-uncertainty) climate-disease-uncertainty $ snakemake --cores 1 --forceall figures_svg

```