## Figure 4

This code reproduces Figure 4.

How to generate the figure:

* First, obtain the following files that were generated earlier: the allele frequency data ()`orch2021_Downsampled_META_Filtered.RData`) and  inbred line SNP tables for each chromosomem (`inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`). Place these files into this folder under `data/raw`.

* Run all numbered R and Python scripts in this directory in order. These scripts will generate small tables in the `plot_data` folder that are used for plotting.

* Run `plot.Rmd` to generate the figure file.