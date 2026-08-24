# Software requirements

The workflows use the following primary software and language environments.
Exact versions reported in the accepted manuscript should take precedence over
this operational list.

## Command-line software

- STACKS (`ref_map.pl`, `gstacks`, and `populations`)
- A reference-read aligner and SAM/BAM utilities used by the STACKS workflow
- VCFtools and/or bcftools
- PLINK
- IQ-TREE 3
- fastStructure
- PopCluster
- BPP
- BEAST and associated tree-log utilities where applicable
- `runeems_snps` and `rEEMSplots`

## R

Scripts use base R plus analysis-specific packages that include `adegenet`,
`ape`, `ggplot2`, `patchwork`, `data.table`, `vcfR`, `hierfstat`, `vegan`,
`sf`, `terra`, `rnaturalearth`, `coda`, and packages loaded explicitly at the
start of individual scripts.

## Python

Python scripts use packages that include `numpy`, `pandas`, `matplotlib`,
`scipy`, `scikit-learn`, `geopandas`, `shapely`, `pyproj`, and `rasterio`, as
loaded explicitly by individual scripts.

## Inputs not included

The archive does not duplicate raw FASTQ files, the reference assembly, large
VCF/PLINK matrices, external geographic rasters, third-party source code, or
MCMC output. Input paths are visible near the start of each script and should
be updated for the deposit or execution environment.

