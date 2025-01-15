## Figure 4

This code reproduces Figure 4.

How to generate the figure:

* First, obtain the following files that were generated earlier: the allele frequency data (`orch2021_Downsampled_META_Filtered.RData`) and inbred line SNP tables for each chromosome (`Orchard2021/inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`). Place the first file under `data/raw/` (`data/raw/orch2021_Downsampled_META_Filtered.RData`). Place the SNP tables under `data/snptables/Orchard2021/` (`inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`).

* Run all numbered R and Python scripts in this directory in order. These scripts will generate small tables in the `plot_data` folder that are used for plotting.

* Run `plot.Rmd` to generate the figure file.

### Programming environment

The R scripts require the following packages: `c("tidyverse", "reticulate", "broom", "egg", "ggnewscale")`.

The Python programming environment should be created with `conda` using the `env.yml` file in this directory.