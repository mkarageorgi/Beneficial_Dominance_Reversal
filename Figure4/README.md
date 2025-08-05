## Figure 4

This code reproduces Figure 4.

How to generate the figure:

### First, obtain the following files that were generated earlier:
- orch2021_Downsampled_META_Filtered.RData - the allele frequency data from the orchard 2021 experiment (available at https://doi.org/10.5061/dryad.w0vt4b937) 
- dm3.fa.out - repeat masker (available at https://doi.org/10.5061/dryad.w0vt4b937) 
- Orchard2021/inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric  - inbred line SNP tables for each chromosome (available at https://doi.org/10.5061/dryad.w0vt4b937)

### Part 1. Quasi-binomial logistic regression to estimate allele frequency differences in malathion-treated and untreated populations during and after the malathion treatment

* Run the code in the `glm` subfolder to generate GLM effect size estimates. For this part of the analysis, obtain and put under `data/raw/` the `orch2021_Downsampled_META_Filtered.RData` and `dm3.fa.out` files as `data/raw/orch2021_Downsampled_META_Filtered.RData` and `data/raw/dm3.fa.out`, respectively.
First, run `01_glm_malation.R`, passing each chromosome arm as an argument:
```
Rscript --vanilla 01_glm_malation.R 2L
Rscript --vanilla 01_glm_malation.R 2R
Rscript --vanilla 01_glm_malation.R 3L
Rscript --vanilla 01_glm_malation.R 3R
Rscript --vanilla 01_glm_malation.R X
```

These scripts are highly parallelized and computationally intensive, so a cluster environment is likely neccessary. Then, run `02_summarize_glm.R`. This script generates a file `data/raw/sigsite_malation.csv`(available at https://doi.org/10.5061/dryad.w0vt4b937) that would be used in downstream computations.

### Part 2. Linked SNP selection procedure to identify SNPs linked (Ace-linked) and not linked (control) to the sweeping R2 and R3 Ace alleles across the genome. Use these sets of SNPs to quantify the extent of both the forward Ace sweep during malathion treatment and the reverse Ace sweep post-treatment

- Place the SNP tables under `data/snptables/Orchard2021/` (`inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`).

* Run all numbered R and Python scripts in this directory in order. These scripts will generate small tables in the `plot_data` folder that are used for plotting.

* Run `plot.Rmd` to generate the figure panels.

### Programming environment

The R scripts require the following packages: `c("tidyverse", "reticulate", "broom", "egg", "ggnewscale", "purrr", "furrr", "biomartr")`.

The Python programming environment should be created with `conda` using the `env.yml` file in this directory.

Code and analysis by Marianna Karageorgi (glm) and Egor Lappo (linked SNPs)
